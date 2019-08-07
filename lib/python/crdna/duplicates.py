

import itertools
import tenkit.lane as tk_lane
import numpy as np
import math
import random
from collections import namedtuple
from crdna.read_filter import stringent_read_filter
from crdna.constants import SELF_FIVE_PRIME_POS_TAG, MATE_FIVE_PRIME_POS_TAG, DUPLICATE_COUNT_TAG
import crdna.bio_io as crdna_io
import martian

N_BINS = 5000
DEFAULT_FRAC_MAX_DISTANCE_DIFFUSION = 0.01
# For optical duplicate detection
OPTICAL_DUPLICATE_DISTANCE=100

# For diffusion duplicate detection, max distance over which diffusion is expected
MAX_DIFFUSION_DUP_DISTANCE=25e3

class Range:
    def __init__(self):
        self.min = None
        self.max = None
    def update(self, val):
        if self.min is None:
            self.min = val
        if self.max is None:
            self.max = val

        self.min = min(self.min, val)
        self.max = max(self.max, val)
    
    def overlap(self, other):
        if self.min > other.max or self.max < other.min:
            return 0.0
        l0 = min(self.min, other.min)
        r0 = max(self.max, other.max)
        l1 = max(self.min, other.min)
        r1 = min(self.max, other.max)
        return float(r1-l1)/float(r0-l0)

class XYrange:
    def __init__(self):
        self.xrange = Range()
        self.yrange = Range()

    def update(self, coords):
        self.xrange.update(coords[0])
        self.yrange.update(coords[1])

    def overlap(self, other):
        return self.xrange.overlap(other.xrange) * self.yrange.overlap(other.yrange)

    def __str__(self):
        return "Xrange : [" + str(self.xrange.min) + ", " + str(self.xrange.max) + "]\t" \
                "Yrange : [" + str(self.yrange.min) + ", " + str(self.yrange.max) + "]"

def estimate_flowcell_geometry(bam_in, lane_coordinate_system):

    bam_in.reset()

    lane_extents = {}
    flowcells = set()
    result = XYrange()

    for read in itertools.islice(bam_in, 100000):
        read_loc = tk_lane.extract_read_position(read)
        loc_coords = lane_coordinate_system.convert_to_lane_coords(read_loc)
        flowcells.add(read_loc.flowcell)
        key = read_loc.flowcell + '_' + read_loc.lane
        if key not in lane_extents:
            lane_extents[key] = XYrange()
        lane_extents[key].update(loc_coords)
        result.update(loc_coords)
    bam_in.reset()

    min_overlap = 1.0
    lanes = lane_extents.keys()
    for i in xrange(len(lanes)):
        for j in xrange(i+1, len(lanes)):
            min_overlap = min(lane_extents[lanes[i]].overlap(lane_extents[lanes[j]]), min_overlap)

    print "Flowcells : ", flowcells
    for k, v in lane_extents.items():
        print k, v

    print "Flowcell overlap = ", min_overlap

    if len(flowcells) > 1:
        print "Multiple flowcells found"

    if min_overlap < 0.9:
        print "Flowcells do not have identical geometry."
        return None

    return {"x": (result.xrange.min, result.xrange.max), "y": (result.yrange.min, result.yrange.max)}

def find_max_distance(flowcell_geometry):
    x = flowcell_geometry["x"]
    y = flowcell_geometry["y"]
    return np.linalg.norm((x[1]-x[0], y[1]-y[0]))

def default_diffusion_threshold(flowcell_geometry):
    if flowcell_geometry is None:
        return MAX_DIFFUSION_DUP_DISTANCE
    return find_max_distance(flowcell_geometry) * DEFAULT_FRAC_MAX_DISTANCE_DIFFUSION

class Bins:
    def __init__(self, flowcell_geometry, nbins=N_BINS):
        max_dist = find_max_distance(flowcell_geometry)
        self.bin_size = round(max_dist / nbins)
        self.nbins = nbins

    def get_bin(self, val):
        assert(val >= 0)
        n = max(min(math.ceil(float(val) / self.bin_size), self.nbins), 1)
        return n * self.bin_size - self.bin_size / 2.0

class BinnedCounts:
    def __init__(self, bin):
        self.distribution = dict()
        self.bin = bin
        for i in xrange(bin.nbins):
            val = i*bin.bin_size + bin.bin_size / 2.0
            self.distribution[bin.get_bin(val)] = 0

    def increment(self, val):
        self.distribution[self.bin.get_bin(val)] += 1


def compute_null_distribution(flowcell_geometry, nsamples=1000, seed=None):
    if seed is not None:
        np.random.seed(np.uint32(seed))
    
    x = flowcell_geometry["x"]
    y = flowcell_geometry["y"]

    def interp(p, factor):
        return p[0] + (p[1] - p[0]) * factor

    cc = np.random.rand(nsamples, 2)

    distribution = BinnedCounts(Bins(flowcell_geometry))

    for i in xrange(nsamples):
        xi = interp(x, cc[i][0])
        yi = interp(y, cc[i][1])
        for j in range(i+1, nsamples):
            xj = interp(x, cc[j][0])
            yj = interp(y, cc[j][1])
            dist = np.linalg.norm((xi-xj, yi-yj))
            distribution.increment(dist)

    return distribution.distribution

# Generator utilities -- move to tenkit?
def consumer(func):
    ''' decorator for initializing a generator consumer function '''
    def start(*args,**kwargs):
        c = func(*args,**kwargs)
        c.next()
        return c
    return start

def broadcast(source, consumers):
    ''' send each item in the source generator to each consumer in the list '''
    for item in source:
        for c in consumers:
            c.send(item)

    for c in consumers:
        c.close()

class DupMode:
    ESTIMATE = 0,
    MARK = 1

class DupSummary:
    def __init__(self, perfect_read_filter, sample_rate, split_bcs, description, lane_coordinate_system, output_bam=None, write_to_stdout=False,\
        mode=DupMode.MARK, flowcell_geometry=None, threshold=MAX_DIFFUSION_DUP_DISTANCE, tag_counts=False):
        ''' Summarize dups at a given subsampling rate, and barcode
            splitting policy.  If an open output_bam pysam.Samfile
            is passed, dups will be marked and reads will be written
            to output_bam '''

        self.perfect_read_filter = perfect_read_filter
        self.require_barcode = perfect_read_filter and split_bcs
        self.sample_rate = sample_rate
        self.split_bcs = split_bcs
        self.description = description
        self.output_bam = output_bam
        self.write_to_stdout = write_to_stdout
        self.lane_coordinate_system = lane_coordinate_system
        self.custom_diffusion_threshold = threshold
        self.tag_counts = tag_counts

        self.mode = mode
        if self.mode == DupMode.MARK:
            # default return value -- will be used if the sample_rate > 1.0
            self.result = (None, None, None, None)
            self.observed_distribution = None

        elif self.mode == DupMode.ESTIMATE:
            self.result = None
            assert flowcell_geometry is not None
            self.observed_distribution = BinnedCounts(Bins(flowcell_geometry))

    @classmethod
    def diffusion_estimator(cls, lane_coordinate_system, flowcell_geometry):
        return cls(True, 1.0, True, "full_use_bcs", lane_coordinate_system, mode=DupMode.ESTIMATE, flowcell_geometry=flowcell_geometry)

    def count_dups_by_distance(self, reads):
        '''Count number of nearby duplicates in a set of reads.  A pair is counted as 1'''
        # Get (flowcell, lane, surface, swath, tile, x, y) tuples for each read
        read_locs = []
        for (key, read, idx) in reads:
            read_loc = tk_lane.extract_read_position(read)
            if read_loc is not None:
                read_locs.append((read_loc, read))

        # Sort by flowcell_lane
        def flowcell_lane(read_loc):
            return "%s_%s" % (read_loc[0].flowcell, read_loc[0].lane)
        read_locs.sort(key=flowcell_lane)
        lane_groups = itertools.groupby(read_locs, flowcell_lane)

        opt_dups_found = 0   # really close dupes
        diff_dups_found = 0  # somewhat close dupes
        custom_diff_dups_found = 0 # Based on the custom threshold

        # Measure distances between all pairs in a lane
        for (lane, lane_reads) in lane_groups:
            lane_reads = list(lane_reads)

            layout = self.lane_coordinate_system.get_layout_for_read_loc(lane_reads[0][0])
            test_dups = layout.has_diffusion_duplicates(MAX_DIFFUSION_DUP_DISTANCE)

            if len(lane_reads) > 100:
                martian.log_info("Got dup cluster of size: %d" % len(lane_reads))
                first_read = lane_reads[0][1]
                martian.log_info("tid: %d, pos: %d, mapq: %d, seq: %s" % (first_read.reference_id, first_read.reference_start, first_read.mapping_quality, first_read.query_sequence))

            opt_dups = set()
            diff_dups = set()
            custom_diff_dups = set()
            dump = []
            cmp_reads = min(200, len(lane_reads))
            lane_loc_coords = [self.lane_coordinate_system.convert_to_lane_coords(loc) for (loc, _) in lane_reads]
            for i in range(cmp_reads):
                loc1, read1 = lane_reads[i]
                lane_loc1 = lane_loc_coords[i]

                for j in range(i+1, len(lane_reads)):
                    loc2, read2 = lane_reads[j]
                    lane_loc2 = lane_loc_coords[j]

                    dist = math.sqrt((lane_loc1[0]-lane_loc2[0])**2 + (lane_loc1[1] - lane_loc2[1])**2)
                    if test_dups and dist < MAX_DIFFUSION_DUP_DISTANCE:
                        diff_dups.add(j)
                        if self.write_to_stdout and j not in diff_dups:
                            dump.append(("%d\t"+("%d\t"*14)) % (dist,
                                                                loc1.surface, loc1.swath, loc1.tile, loc1.x, loc1.y, lane_loc1[0], lane_loc1[1],
                                                                loc2.surface, loc2.swath, loc2.tile, loc2.x, loc2.y, lane_loc2[0], lane_loc2[1]))

                    if dist < OPTICAL_DUPLICATE_DISTANCE:
                        opt_dups.add(j)
                    
                    if dist < self.custom_diffusion_threshold:
                        custom_diff_dups.add(j)

            if self.write_to_stdout and len(diff_dups) >= 2:
                for x in dump:
                    print ("%d\t%s" % (len(diff_dups), x))

            diff_dups_found += len(diff_dups)
            opt_dups_found += len(opt_dups)
            custom_diff_dups_found += len(custom_diff_dups)

        return (opt_dups_found, diff_dups_found, custom_diff_dups_found)


    def update_observed_distribution(self, reads):

        # Get (flowcell, lane, surface, swath, tile, x, y) tuples for each read
        read_locs = []
        for (key, read, idx) in reads:
            read_loc = tk_lane.extract_read_position(read)
            if read_loc is not None:
                read_locs.append((read_loc, read))

        # Sort by flowcell
        def flowcell(read_loc):
            return "%s" % (read_loc[0].flowcell)
        read_locs.sort(key=flowcell)
        lane_groups = itertools.groupby(read_locs, flowcell)

        # Measure distances between all pairs
        for (lane, lane_reads) in lane_groups:
            lane_reads = list(lane_reads)
            cmp_reads = min(200, len(lane_reads))
            lane_loc_coords = [self.lane_coordinate_system.convert_to_lane_coords(loc) for (loc, _) in lane_reads]
            for i in range(cmp_reads):
                loc1, read1 = lane_reads[i]
                lane_loc1 = lane_loc_coords[i]

                for j in range(i+1, len(lane_reads)):
                    loc2, read2 = lane_reads[j]
                    lane_loc2 = lane_loc_coords[j]

                    dist = math.sqrt((lane_loc1[0]-lane_loc2[0])**2 + (lane_loc1[1] - lane_loc2[1])**2)
                    self.observed_distribution.increment(dist)

    @consumer
    def read_consumer(self):
        read_footprint = namedtuple("read_footprint", ["barcode", "tid", "self_five_prime_pos", "mtid", "mate_five_prime_pos"])
        read_tuple = namedtuple("read_tuple", ["footprint", "read", "num"])

        def process_read_block(reads, count_hist, optical_dup_hist, diffusion_dup_hist, custom_diffusion_dup_hist):
            ''' dedup a block of reads, then write to BAM in original order '''
            read_tuples = []
            for read in reads:
                # Don't consider unmapped reads
                if read.is_unmapped or read.mate_is_unmapped or read.is_secondary:
                    continue

                bc_sequence = None

                if self.split_bcs:
                    bc_sequence = crdna_io.get_read_barcode(read)

                footprint = read_footprint(bc_sequence, read.reference_id, read.get_tag(SELF_FIVE_PRIME_POS_TAG), read.next_reference_id, read.get_tag(MATE_FIVE_PRIME_POS_TAG))

                # Include the read and the original index so that we can get back to the original order
                read_tuples.append(read_tuple(footprint, read, len(read_tuples)))

            # Sort by the read tuple -- need to do this for the groupby to get all the common items together
            read_tuples.sort(key=lambda x: x.footprint)

            # Group the reads by the read_tuple
            dup_groups = itertools.groupby(read_tuples, lambda x: x.footprint)

            for (key, dup_group) in dup_groups:
                # Note how many dups we have
                dup_group = list(dup_group)
                n_dups = len(dup_group)

                if self.mode == DupMode.ESTIMATE:
                    self.update_observed_distribution(dup_group)
                    continue

                if n_dups > 1:
                    optical_dups, diffusion_dups, custom_dups = self.count_dups_by_distance(dup_group)
                else:
                    optical_dups = 0
                    diffusion_dups = 0
                    custom_dups = 0

                non_proximal_dups = n_dups - custom_dups

                # If we are splitting on bcs, then only counts stats for read groups with BCs
                group_bc = key[0]
                if not self.split_bcs or group_bc is not None:
                    count_hist[non_proximal_dups] = count_hist.setdefault(non_proximal_dups,0) + 1
                    if optical_dups > 0:
                        optical_dup_hist['count'] = optical_dup_hist.setdefault('count',0) + optical_dups
                    diffusion_dup_hist['count'] = diffusion_dup_hist.setdefault('count',0) + diffusion_dups
                    custom_diffusion_dup_hist['count'] = custom_diffusion_dup_hist.setdefault('count',0) + custom_dups

                # We pick the read with the smallest query name as non-dup.
                # Its pair is guaranteed to be picked as non-dup in it's own dup group, in our
                # scheme of using five prime positions of read and mate for sorting
                min_qname_read = min(dup_group, key=lambda x: x.read.query_name)

                if self.tag_counts:
                    min_qname_read.read.set_tag(DUPLICATE_COUNT_TAG, non_proximal_dups-1) # -1 because 1 is unique

                # Mark dups
                if self.output_bam:
                    for i in range(0, n_dups):
                        dup_group[i].read.is_duplicate = True
                    min_qname_read.read.is_duplicate = False

            if self.output_bam:
                for read in reads:
                    self.output_bam.write(read)

            # done processing block of reads
            return


        # bam is sorted by SELF_FIVE_PRIME_POS tag, chrom and pos.
        current_bam_key = (-1, -1)
        current_reads = []

        # we will track the histogram of the reads/molecule
        molecule_count_hist = {}
        optical_dup_hist = {'count':0}
        diffusion_dup_hist = {'count':0}
        custom_diffusion_dup_hist = {'count':0}

        try:
            while True:
                # accept the next read
                read = (yield)

                # Sample to the requested subsample rate
                skip = (self.sample_rate < 1.0) and (random.random() > self.sample_rate)

                # Apply the perfect read filter - the perfect read filter requires barcodes if they are
                # being split on
                skip = skip or (self.perfect_read_filter and not stringent_read_filter(read, self.require_barcode))

                new_bam_key = (read.reference_id, read.get_tag(SELF_FIVE_PRIME_POS_TAG))

                # If the dup group gets extremely large we can run out of memory.
                # Process things in groups of 500K to prevent memory blow-up
                # May cause us to miss a few dups, but it doesn't really matter in these crazy regions
                if new_bam_key != current_bam_key or len(current_reads) > 500000:
                    # Otherwise process our block of reads
                    process_reads = current_reads
                    current_reads = []
                    current_bam_key = new_bam_key

                    if len(process_reads) > 0:
                        process_read_block(process_reads, molecule_count_hist, optical_dup_hist, diffusion_dup_hist, custom_diffusion_dup_hist)

                if not skip:
                    # accumulate block of reads with same start position
                    current_reads.append(read)

                    # Short circuit the bolus of reads
                    # that have no mapping info --
                    # they pass straight through
                    if current_bam_key == (-1, -1):
                        current_reads = []
                        if self.output_bam:
                            self.output_bam.write(read)


        except GeneratorExit:
            # Finish up final batch
            process_read_block(current_reads, molecule_count_hist, optical_dup_hist, diffusion_dup_hist, custom_diffusion_dup_hist)

            # save the result for later use
            if self.mode == DupMode.MARK:
                self.result = (molecule_count_hist, optical_dup_hist, diffusion_dup_hist, custom_diffusion_dup_hist)
            else:
                self.result = self.observed_distribution.distribution
            return

