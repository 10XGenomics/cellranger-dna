#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Put the sample index barcode and 10X barcode into read tags in a BAM file.
# Screen each 10X barcode against the given barcode whitelist
#
import gzip
import itertools
import random
import numpy as np
import array
import json
import tenkit.seq as tk_seq
import tenkit.bam as tk_bam
import tenkit.fasta as tk_fasta
import crdna.bio_io as crdna_io
from crdna.read_filter import stringent_read_filter
from crdna.constants import PROCESSED_BARCODE_TAG, RAW_BARCODE_TAG, RAW_BARCODE_QUAL_TAG, SELF_FIVE_PRIME_POS_TAG, MATE_FIVE_PRIME_POS_TAG
from tenkit.constants import SAMPLE_INDEX_TAG, SAMPLE_INDEX_QUAL_TAG, ILLUMINA_QUAL_OFFSET, TRIM_TAG, TRIM_QUAL_TAG
import barcodes.utils as bc_utils


import martian

__MRO__ = """
 stage ATTACH_BCS(
     in  string barcode_whitelist,
     in  bam[]  align,
     in  map[]  chunks,
     in  bool   paired_end,
     in  bool   exclude_non_bc_reads,
     in  float  bc_confidence_threshold,
     in  json   bc_counts,
     out bam[]  outputs,
     out int    perfect_read_count,
     src py     "stages/reads/attach_bcs",
 ) split using (
     in  bam    align_chunk,
     in  map    chunk,
     in  int    chunk_index,
     out bam    output,
 )
"""

class GlobalFivePrimePosTagger:
    """
    Computes the location of the 5' most end of each read in chromosome-stiched global
    reference coordinates including soft-clipped regions. This will be useful in
    marking duplicates accurately.
    Attributes:
    """

    def __init__(self, bam_in):
        """
        Args:
            bam_in(pysam.AlignmentFile): Pysam bam reader object
        """
        chrom_lengths = bam_in.lengths # Lengths of chromosomes in the same order as pysam.AlignmentFile.references
        self.offsets = [0] * bam_in.nreferences
        for i in xrange(1, bam_in.nreferences):
            self.offsets[i] = self.offsets[i-1] + chrom_lengths[i-1]

    def compute_global_coords(self, read):
        """
        For unmapped read the global coordinate is -1
        For mapped reads, compute the 5' position in chromosome coordinates and add the offset
        Args:
            read(pysam.AlignedSegment): Pysam read object
        """
        MAX_I32 = 2147483647 # BAM tags are unfortunately of type int32, so wrap around
        # Overlaps are okay as samtools sort will use read position as a second index
        # when we ask it to sort by some tag
        if read.is_unmapped:
            return -1
        else:
            return (self.compute_five_prime_coords(read) + self.offsets[read.reference_id]) % MAX_I32

    def compute_five_prime_coords(self, read):
        """
        Computes the 5' position in chromosome coordinates
        Assumes that the read is mapped and the CIGAR exists.
        Args:
            read(pysam.AlignedSegment): Pysam read object
        """
        SOFT_CLIP = 4
        cigar = read.cigartuples
        if read.is_reverse:
            # Add the suffix clip to the alignment end position
            suffix_clip = sum([x[1] for x in itertools.takewhile(lambda x: x[0]==SOFT_CLIP, reversed(cigar))])
            return read.reference_end + suffix_clip
        else:
            # Subtract the prefix clip from the alignment position
            prefix_clip = sum([x[1] for x in itertools.takewhile(lambda x: x[0]==SOFT_CLIP, cigar)])
            return read.reference_start - prefix_clip

    def tag_reads(self, reads):
        """
        Adds an integer tag to the read, denoting the global position of the most 5' base
        Args:
            reads(pysam.AlignedSegment): List of Pysam read objects sharing the same qname
        """
        coords = {True: -1, False: -1}
        for read in reads:
            assert read.query_name == reads[0].query_name
            if not read.is_secondary:
                coords[read.is_read1] = self.compute_global_coords(read)

        for read in reads:
            read.set_tag(SELF_FIVE_PRIME_POS_TAG, self.compute_global_coords(read))
            read.set_tag(MATE_FIVE_PRIME_POS_TAG, coords[not read.is_read1])


def split(args):
    ''' Attach BCS to each chunk of the input files '''

    chunk_defs = [
            {'chunk': fastq_chunk, 'align_chunk': aln_chunk, 'chunk_index': i}
            for i, (fastq_chunk, aln_chunk) in enumerate(zip(args.chunks, args.align))]
    return {'chunks': chunk_defs}



def join(args, outs, chunk_defs, chunk_outs):
    ''' Pass through the BAM files for sorting '''
    outs.outputs = [chunk.output for chunk in chunk_outs]
    outs.perfect_read_count = sum([chunk.perfect_read_count for chunk in chunk_outs])


def main(args, outs):
    """ Attaches barcodes. Attaches raw barcode to RAW_BC tag and filters those to form set of PROCESSES_BARCODES """
    # this silences a weird non-failure in --strict=error mode
    # TODO(lhepler): remove this when martian upstream handles this itself
    outs.outputs = []

    chunk = args.chunk

    bam_in = tk_bam.create_bam_infile(args.align_chunk)
    bc_spec = "{}:{}".format(RAW_BARCODE_TAG, RAW_BARCODE_QUAL_TAG)

    # only comment the first chunk, otherwise later merge will duplicate the comments and could lead to:
    # samtools merge ... : '[finish_merged_header] Output header text too long'
    if args.chunk_index > 0:
        COs = None
    elif chunk['trim']:
        COs = ['10x_bam_to_fastq:R1({},TR:TQ,SEQ:QUAL)'.format(bc_spec), '10x_bam_to_fastq:R2(SEQ:QUAL)', '10x_bam_to_fastq:I1(BC:QT)']
    else:
        COs = ['10x_bam_to_fastq:R1({},SEQ:QUAL)'.format(bc_spec), '10x_bam_to_fastq:R2(SEQ:QUAL)', '10x_bam_to_fastq:I1(BC:QT)']

    bam_out, tids = tk_bam.create_bam_outfile(outs.output, None, None, template=bam_in, pgs=[tk_bam.make_pg_header(martian.get_pipelines_version(), "attach_bcs")], cos = COs)

    gp_tagger = GlobalFivePrimePosTagger(bam_in)

    if args.barcode_whitelist is None or args.bc_counts is None:
        # If there's no whitelist or counts then all high quality BC reads get allowed.
        barcode_whitelist = None
        wl_idxs = None
        bc_dist = None
    else:
        barcode_whitelist = bc_utils.load_barcode_whitelist(args.barcode_whitelist)

        # Load the bc counts for this GEM group
        counts = json.load(open(args.bc_counts, 'r'))
        counts = counts[str(chunk['gem_group'])]['bc_counts']

        # Prior distribution over barcodes, with pseudo-count
        bc_dist = np.array(counts, dtype=np.float) + 1.0
        bc_dist = bc_dist / bc_dist.sum()
        wl_idxs = { bc:idx for (idx,bc) in enumerate(sorted(list(barcode_whitelist))) }

    # set random seed to get deterministic subsampling
    random.seed(0)

    def open_maybe_gzip(fn):
        if fn[-2:] == "gz":
            return gzip.open(fn)
        else:
            return open(fn)

    if chunk['barcode']:
        processed_barcode_iter = get_raw_processed_barcodes(open_maybe_gzip(chunk['barcode']), barcode_whitelist, args.bc_confidence_threshold, chunk['gem_group'], chunk['barcode_reverse_complement'], wl_idxs, bc_dist)
        require_barcode_for_stringent = True
    else:
        processed_barcode_iter = itertools.repeat(None)
        require_barcode_for_stringent = False

    if chunk['trim']:
        trim_iter = tk_fasta.read_generator_fastq(open_maybe_gzip(chunk['trim']), paired_end=True)
    else:
        trim_iter = itertools.repeat(None)

    if chunk['sample_index']:
        sample_index_iter = tk_fasta.read_generator_fastq(open_maybe_gzip(chunk['sample_index']))
    else:
        sample_index_iter = itertools.repeat(None)

    iters = itertools.izip(processed_barcode_iter, sample_index_iter, trim_iter)

    # First read
    read = bam_in.next()

    # Number of perfect reads -- used to compute down-sampling rates in mark_duplicates
    perfect_read_count = 0

    # Due to secondary alignments, we must apply the tags to all
    # reads with the same cluster name.
    for (barcode_info, sample_index_info, trim_info) in iters:
        tags = []
        read_name = None

        if read is None:
            break

        if barcode_info:
            (bc_read_name, raw_bc_seq, processed_bc_seq, raw_bc_qual) = barcode_info
            tags.append((RAW_BARCODE_TAG, raw_bc_seq))
            tags.append((RAW_BARCODE_QUAL_TAG, raw_bc_qual))
            if processed_bc_seq is not None:
                tags.append((PROCESSED_BARCODE_TAG, processed_bc_seq))
            read_name = bc_read_name.split()[0]


        if sample_index_info:
            (si_read_name, seq, qual) = sample_index_info
            tags.append((SAMPLE_INDEX_TAG, seq))
            tags.append((SAMPLE_INDEX_QUAL_TAG, qual))

            if read_name != None:
                if si_read_name.split()[0] != read_name:
                    martian.log_info("mismatch: si_read_name: %s, bam_read_name: %s" % (si_read_name, read_name))
                assert(si_read_name.split()[0] == read_name)
            else:
                read_name = si_read_name.split()[0]

        r1_tags = tags
        r2_tags = list(tags)

        if trim_info:
            (trim1_read_name, trim1_seq, trim1_qual, trim2_read_name, trim2_seq, trim2_qual) = trim_info
            if len(trim1_seq) > 0:
                r1_tags.append((TRIM_TAG, trim1_seq))
                r1_tags.append((TRIM_QUAL_TAG, trim1_qual))
            if len(trim2_seq) > 0:
                r2_tags.append((TRIM_TAG, trim2_seq))
                r2_tags.append((TRIM_QUAL_TAG, trim2_qual))

        reads_attached = 0
        reads_to_attach = []

        while read.qname == read_name or read_name == None:
            tags = r1_tags if read.is_read1 else r2_tags
            if len(tags) > 0:
                existing_tags = read.tags
                existing_tags.extend(tags)
                read.tags = existing_tags

            if not (read_name is None):
                assert(read.qname == read_name)

            if reads_to_attach and (read.query_name != reads_to_attach[0].query_name or reads_to_attach[0].query_name is None):
                gp_tagger.tag_reads(reads_to_attach)
                reads_attached += len(reads_to_attach)
                for r in reads_to_attach:
                    if stringent_read_filter(r, require_barcode_for_stringent):
                        perfect_read_count += 1

                    if args.exclude_non_bc_reads:
                        if not(crdna_io.get_read_barcode(r) is None):
                            bam_out.write(r)
                    else:
                        bam_out.write(r)
                reads_to_attach = []

            reads_to_attach.append(read)

            try:
                read = bam_in.next()

            except StopIteration:
                read = None
                break

        gp_tagger.tag_reads(reads_to_attach)
        reads_attached += len(reads_to_attach)
        for r in reads_to_attach:
            if stringent_read_filter(r, require_barcode_for_stringent):
                perfect_read_count += 1

            if args.exclude_non_bc_reads:
                if not(crdna_io.get_read_barcode(r) is None):
                    bam_out.write(r)
            else:
                bam_out.write(r)

        # We may have more than 2 reads is there was a
        # secondary alignment, but less than 2 means
        # something went wrong
        assert(reads_attached >= 2)


    outs.perfect_read_count = perfect_read_count
    bam_out.close()


def get_raw_processed_barcodes(barcode_file, barcode_whitelist, bc_confidence_threshold, gem_group, barcodes_reverse_complement, wl_idxs, wl_dist):
    """ Stream the barcodes and the 'processed' barcode """
    bc_iterator = tk_fasta.read_generator_fastq(barcode_file)

    gem_group_str = "-" + str(gem_group)

    for (name, seq, qual) in bc_iterator:
        if barcodes_reverse_complement:
            seq = tk_seq.get_rev_comp(seq)
            qual = qual[::-1] #reverse qual string
        # Check for valid bc sequences
        if barcode_whitelist is None:
            # No whitelist case -- attach BC if there are no Ns
            if not ('N' in seq):
                processed_bc = seq + gem_group_str
                yield (name, seq, processed_bc, qual)
            else:
                yield (name, seq, None, qual)
        else:
            # whitelist case -- attach bc if posterior probability of best
            # BC sequence exceeds the confidence threshold

            bc_seq = handle_10x_barcode(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist)
            if bc_seq is None:
                yield (name, seq, None, qual)
            else:
                processed_bc = bc_seq + gem_group_str
                yield (name, seq, processed_bc, qual)

def handle_10x_barcode(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist):
    '''Compute the MAP whitelist BC sequence using the read QV * bc distribution
       model. If the MAP probability exceeds the bc_confidence_threshold,
       return that sequence, otherwise return None.  In the case the BC is on
       the whitelist, consider Hamming distance=2 neighbors, if not, consider
       Hamming distance=1 neighbors'''
    if seq in wl_idxs:
        return check_correct_bc(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist)
    else:
        return correct_bc_error(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist)


def correct_bc_error(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist):
    '''Attempt to correct an incorrect BC sequence by computing
       the probability that a Hamming distance=1 BC generated
       the observed sequence, accounting for the prior distribution
       of the whitelist barcodes (wl_dist), and the QV of the base
       that must have been incorrect'''

    # QV values
    qvs = np.fromstring(qual, dtype=np.byte) - ILLUMINA_QUAL_OFFSET

    # Char array of read
    a = array.array('c', seq)

    # Likelihood of candidates
    wl_cand = []
    likelihoods = []

    # Enumerate Hamming distance 1 sequences - if a sequence
    # is on the whitelist, compute it's likelihood.
    for pos in range(len(a)):
        existing = a[pos]
        for c in ['A', 'C', 'G', 'T']:
            if c == existing:
                continue
            a[pos] = c
            test_str = a.tostring()

            idx = wl_idxs.get(test_str)
            if idx is not None:
                # prior probability of this BC
                p_bc = wl_dist[idx]

                # probability of the base error
                edit_qv = min(40.0, float(qvs[pos]))
                p_edit = 10.0**(-edit_qv/10.0)
                wl_cand.append(test_str)
                likelihoods.append(p_bc * p_edit)

        a[pos] = existing

    if len(likelihoods) > 0:
        posterior = np.array(likelihoods)
        posterior /= posterior.sum()
        pmax = posterior.max()
        if pmax > bc_confidence_threshold:
            return wl_cand[np.argmax(posterior)]

    return None


def check_correct_bc(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist):
    '''Attempt to correct an incorrect BC sequence by computing
       the probability that a Hamming distance=1 BC generated
       the observed sequence, accounting for the prior distribution
       of the whitelist barcodes (wl_dist), and the QV of the base
       that must have been incorrect'''

    # QV values
    qvs = np.fromstring(qual, dtype=np.byte) - ILLUMINA_QUAL_OFFSET

    # Only examine the barcode if there is a QV <= 24
    if (qvs > 24).all():
        return seq

    # Char array of read
    a = array.array('c', seq)

    # Likelihood of candidates
    wl_cand = []
    likelihoods = []

    # include the un-edited case here -- this assumes the pErrs are small
    original_bc = a.tostring()
    wl_cand.append(original_bc)
    likelihoods.append(wl_dist[wl_idxs[original_bc]])


    # Enumerate Hamming distance 1 sequences - if a sequence
    # is on the whitelist, compute it's likelihood.
    for pos in range(len(a)):
        existing = a[pos]
        for c in ['A', 'C', 'G', 'T']:
            if c == existing:
                continue
            a[pos] = c

            for pos2 in range(pos+1, len(a)):
                existing2 = a[pos2]
                for c2 in ['A', 'C', 'G', 'T']:
                    if c2 == existing2:
                        continue
                    a[pos2] = c2
                    test_str = a.tostring()

                    idx = wl_idxs.get(test_str)
                    if idx is not None:
                        # prior probability of this BC
                        p_bc = wl_dist[idx]

                        # probability of the base errors
                        edit_qv1 = min(33.0, max(3.0, float(qvs[pos])-1.0))
                        edit_qv2 = min(33.0, max(3.0, float(qvs[pos2])-1.0))
                        p_edit = 10.0**(-edit_qv1/10.0) * 10.0**(-edit_qv2/10.0)
                        likelihoods.append(p_bc * p_edit)

                        wl_cand.append(test_str)

                a[pos2] = existing2
        a[pos] = existing

    if len(likelihoods) > 0:
        posterior = np.array(likelihoods)
        posterior /= posterior.sum()
        pmax = posterior.max()
        if pmax > bc_confidence_threshold:
            return wl_cand[np.argmax(posterior)]

    return None
