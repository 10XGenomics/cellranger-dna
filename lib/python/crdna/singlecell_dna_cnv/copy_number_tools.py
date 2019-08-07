#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#import tenkit.constants
from tenkit.bio_io import get_bed_iterator
import numpy
import copy
import math

#...............................................................................
def deneighbor(positions, cutoff):
    pos = positions[:]
    pos.sort()
    n_reads = len(positions)
    if n_reads > 1:
        delta = [x - pos[i - 1] for i, x in enumerate(pos)][1:]
        for i in range(n_reads - 1, 0, -1):
            if delta[i - 1] < cutoff:
                del(pos[i])
            # if delta
        # for i
    # if n_reads
    return pos
# deneighbor

#...............................................................................
def make_bins(contig_length, bin_size):
    bins = range(0, int(contig_length), int(bin_size))
    if bins[-1] < contig_length:
        bins.append(contig_length)
    # if bins
    return bins
# make_bins

#...............................................................................
def generate_raw_profile(bam_in, reads, contigs, bins, deneighbor_cutoff=1e3):
    raw_profile = {}
    #read_positions = dict(zip(contigs.keys(), [] * len(contigs.keys())))
    read_positions = {}
    for contig_name in contigs.keys():
        read_positions[contig_name] = []
    for read in reads:
        fail = (read.is_unmapped 
                or (not read.is_read1) 
                or read.is_read2 
                or read.is_qcfail 
                or read.is_duplicate 
                or read.is_qcfail 
                or read.is_secondary 
                or read.is_supplementary 
                # or (read.mapping_quality < tenkit.constants.HIGH_CONF_MAPQ)
                )
        if not fail:
            try:
                #read_positions[bam_in.get_reference_name(read.tid)].append(read.pos)
                read_positions[read.reference_name].append(read.pos)
            except KeyError:
                print("KeyError: chr=%s, %s, tid=%d" % (bam_in.get_reference_name(read.tid), read.reference_name, read.tid))
                pass
            # try except
        # if !unmapped
    # for read
    #print(read_positions)
    for contig_name in contigs.keys():
        #print("BINS ", contig_name, ": ", bins[contig_name])
        #
        # Dedupe - this does not change anything!
        # read_positions[contig_name] = list(set(read_positions[contig_name]))
        #
        # Deneighbor
        deneighbored_read_positions = deneighbor(read_positions[contig_name], cutoff=deneighbor_cutoff)
        #
        raw_profile[contig_name] = numpy.histogram(
            deneighbored_read_positions,
            bins=bins[contig_name])[0]
    # for contig
    return raw_profile
# generate_raw_profile

#...............................................................................
def mask_profile(profile, bin_size, mask_bed_file, cutoff=0, keep_gt_cutoff=True):
    # Auxiliary function
    def mask_bin(masked_profile, chrom, start, end, bin_size):
        #bin_size = end - start + 1
        bin_index = int(math.floor(start / bin_size))
        try:
            masked_profile.get(chrom)[bin_index] = float("NaN")
        except Exception as e:
            #print(traceback.format_exc())
            print(chrom + ":" + str(start) + "-" + str(end) + " not found")
            print(e)
        # try except
    # mask_bin
    #
    masked_profile = copy.deepcopy(profile)
    # TODO: Implement
    # TODO: Test - the input dictionaries profile and mask must have identical keys
    # TODO: Test - for each key (chrName), the arrays profile[chrName] and mask[chrName]
    #       must have identical lengths
    mask_iterator = get_bed_iterator(bed_file_name=mask_bed_file, locus=None)
    for (chrom, start, end, value) in mask_iterator:
        #print("%s\t%d\t%d\t%f\n" % (chrom, start, end, value))
        if value > cutoff:
            if keep_gt_cutoff:
                pass # keep bins above cutoff
            else:
                # Mask bins above cutoff
                mask_bin(masked_profile, chrom, start, end, bin_size)
            # if keep_gt_cutoff
        else:
            if keep_gt_cutoff:
                # Mask bins below cutoff
                mask_bin(masked_profile, chrom, start, end, bin_size)
            else:
                pass # keep bins below cutoff
            # if keep_gt_cutoff
        # if value else
    # for chrom
    return masked_profile
# mask_profile

#...............................................................................
def count_reads_and_bins(raw_profile):
    n_reads = 0
    n_bins = 0
    for chromosome in raw_profile:
        n_reads += numpy.nansum(raw_profile[chromosome])
        # TODO: Exclude centromeres and N-base bins from the total bin count:
        # TODO: n_bins += numpy.nansum(mask[chromosome]) 
        n_bins += len(raw_profile[chromosome])
    # for chromosome
    return (n_reads, n_bins)
# count_reads_and_bins

#...............................................................................
def scale_raw_profile_all_chr(raw_profile):
    scaled_profile = copy.deepcopy(raw_profile)
    (n_reads, n_bins) = count_reads_and_bins(raw_profile)
    print "n_reads=%d, n_bins=%d" % (n_reads, n_bins)
    if (n_reads > 0) and (n_bins > 0):
        scale = float(n_reads) / float(n_bins) # Average number of reads per bin
        for chromosome in scaled_profile:
            chr_profile = [
                float(x) / scale for x in scaled_profile[chromosome]]
            scaled_profile[chromosome] = chr_profile
        # for chromosome
    # if n_reads
    return scaled_profile
# scale_raw_profile_all_chr

#...............................................................................
def scale_raw_profile(raw_profile, method="all"):
    if method == "all":
        return(scale_raw_profile_all_chr(raw_profile))
    elif method == "autosomes":
        # TODO: Implement
        pass
    else:
        # TODO: Error message
        pass
    # if method else
    return raw_profile
# scale_raw_profile

#...............................................................................
def normalize(scaled_profile):
    # TODO: Implement
    normalized_profile = scaled_profile
    return normalized_profile
# normalize

#...............................................................................
def aggregate_chr_bins(chr_profile, n_merge):
    if n_merge <= 0:
        # TODO: warning!
        return chr_profile
    # if n_merge
    n_bins = len(chr_profile)
    aggregated_profile = []
    for i in range(0, n_bins, n_merge):
        upper_bound = i + n_merge
        aggregated_profile.append(numpy.nansum(chr_profile[i:upper_bound]))
    # for i
    rescaled_profile = [x / n_merge for x in aggregated_profile]
    return rescaled_profile
# aggregate_chr_bins

#...............................................................................
def aggregate_bins(profile, n_merge):
    aggregated_profile = {}
    for chromosome in profile:
        aggregated_profile[chromosome] = aggregate_chr_bins(profile[chromosome], n_merge)
    # for chromosome
    return aggregated_profile
# aggregate_bins

#...............................................................................
def get_all_chromosomes(profile):
    chr_set = list(set([x.keys() for x in profile.keys()]))
    return chr_set
# get_all_chromosomes

#...............................................................................
def aggregate_barcodes(profiles, barcodes):
    merged_profile = {}
    initialize = True
    count = 0
    for bc in barcodes:
        if profiles.get("Profiles").has_key(bc):
            count = count + 1
            if initialize:
                for chromosome in profiles.get("Profiles")[bc]:
                    merged_profile[chromosome] = profiles.get("Profiles")[bc][chromosome]
                # for chromosome
                initialize = False
            else:
                for chromosome in merged_profile:
                    if profiles.get("Profiles")[bc].has_key(chromosome):
                        merged_profile[chromosome] = [
                            merged_profile[chromosome][i] + profiles.get("Profiles")[bc][chromosome][i]
                            for i in range(0,len(merged_profile[chromosome]))]
                    # if profiles
                # for chromosome
            # if initialize else
        # if profiles
    # for bc
    if count > 0:
        rescaled_profile = {}
        for chromosome in merged_profile:
            rescaled_profile[chromosome] = [x / count for x in merged_profile[chromosome]]
        # for chromosome
    else:
        # TODO: warning!
        rescaled_profile = merged_profile
    # if count else
    return rescaled_profile
# aggregate_barcodes

## compute diff variance ratio for a given ncell x nbin matrix
## returns a numpy vector of length ncell
def compute_diff_variance_ratio( M ):
    Mnorm = M/M.mean(axis=1)[:, numpy.newaxis]
    ## dev from mean > 3 sig is allowed if it occurs in less than 
    ## 50 % of cells (very conservative)
    spike_mask = (numpy.abs( M - M.mean(axis=1)[:, numpy.newaxis] )
                  > 3*numpy.std(M, axis=1)[:, numpy.newaxis]).sum(axis=0) < 0.5*M.shape[0]
    Mnorm = Mnorm[:, spike_mask]
    dMnorm = numpy.sort(numpy.diff( Mnorm, axis=1 ), axis=1)
    if dMnorm.shape[1] > 4:
        dMnorm = dMnorm[:, 2:-2]

    S = []
    for i in xrange(Mnorm.shape[0]):
        S.append(2*numpy.std(Mnorm[i])**2/numpy.std(dMnorm[i])**2 - 1)
    return numpy.array(S)

## aggregates counts over a given window
## input can be a list of numbers or a numpy array
## output is a numpy array of size approx len(input)/window
def aggregate_counts( vec, window ):
    assert window > 1
    assert len(vec) >= window
    n = (len(vec)-1)/window+1
    agg_vec = numpy.zeros(n, vec.dtype)
    ai = 0
    for i in xrange(0, len(vec), window):
        high = min(i + window, len(vec))
        agg_vec[ai] = vec[i:high].sum()
        ai += 1
    return agg_vec

## aggregates counts that are in a matrix across columns
## input is a numpy matrix of shape (ncell, nbin)
## output is a numpy matrix of shape approx (ncell, nbin/window)
def aggregate_matrix( M, window ):
    assert window > 1
    n = (M.shape[1]-1)/window+1
    agg_mat = numpy.zeros((M.shape[0], n), M.dtype)
    ai = 0
    for i in xrange(0, M.shape[1], window):
        high = min(i + window, M.shape[1])
        agg_mat[:, ai] = M[:, i:high].sum(axis=1)
        ai += 1
    return agg_mat

## from a ncell x nbin matrix per chromosome
## trim away outliers using when counts are outside
## using MAD score of mad_score
def get_dpcvs( M, mad_score = 6 ):
    assert len(M.shape) == 2
    ncells = M.shape[0]
    
    ## filter out outliers using (gaussian) MAD
    median = numpy.median(M, axis=1)[:, numpy.newaxis]
    mads = 1.4826*numpy.median(numpy.abs(M - median), axis=1)
    outlier_mask = numpy.array([numpy.abs(M[i] - median[i]) < mad_score*mads[i] 
                             for i in xrange(M.shape[0])])

    ## are there repeat outliers across > 25 % of cells?
    ## if so, mask those bins
    bad_frac = numpy.sum((~outlier_mask).astype(int), axis=0).astype(float)/M.shape[0]
    goodbins = (outlier_mask) & (bad_frac < 0.25)

    m = numpy.zeros(ncells)
    v = numpy.zeros(ncells)
    dpcvs = numpy.zeros(ncells)
    for i in xrange(ncells):
        m[i] = M[i][goodbins[i]].mean()
        v[i] = M[i][goodbins[i]].std()**2

    fail_mask = numpy.logical_or(m < 1e-6, v < m)
    dpcvs[fail_mask] = numpy.nan
    dpcvs[~fail_mask] = numpy.sqrt( (v-m)[~fail_mask] ) / m[~fail_mask]
    return dpcvs


