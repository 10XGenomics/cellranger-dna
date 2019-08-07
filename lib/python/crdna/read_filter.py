#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
__doc__="""Overrides tenkit.read_filter functions and provides other helper functions."""
import crdna.constants as crdna_constants
import tenkit.constants

from crdna.bio_io import get_read_barcode

def stringent_read_filter(bam_read, require_barcode):
    ''' Test for a very high-quality read. Only reads satisfying this predicate are
        used when computing summary dup rates to avoid spurious duplicates.
        Reads must have a high MAPQ, a perfect cigar, with a fragment size somewhat
        longer than the read length '''

    return bam_read.mapq >= tenkit.constants.HIGH_CONF_MAPQ and \
            not bam_read.is_secondary and \
           ((abs(bam_read.pos - bam_read.mpos) >= tenkit.constants.MIN_MATE_OFFSET_DUP_FILTER and \
                   bam_read.tid == bam_read.rnext) or bam_read.tid != bam_read.rnext) and \
           len(bam_read.cigar) == 1 and (not require_barcode or get_read_barcode(bam_read) is not None)

def inserts_per_alignment(rec):
    """
    Returns how much a given aligned read contributes towards the insert profile
    based on the configuration of the read and its mate. If a read pair maps
    perfectly, then the contribution is 0.5 per read. The contribution is 0 for
    a wasted read.
    """
    if rec.is_secondary or rec.is_supplementary:
        ## we are ignoring alternate alignments
        return 0.0
    if (rec.is_unmapped or
        rec.mapping_quality < crdna_constants.PROFILE_MAPQ_THRESHOLD or
        rec.is_duplicate):
        ## if read is unmapped, poor mapping quality or dup
        return 0.0
    if rec.mate_is_unmapped:
        ## if mate is unmapped
        return 1.0
    if (rec.reference_id == rec.next_reference_id):
        if rec.template_length > crdna_constants.PROFILE_INSERT_THRESHOLD:
            ## both reads on the same chromosome, very far away, and in ANY orientation
            return 1.0
        else:
            ## both reads on the same chromosome, close by, and in ANY orientation
            return 0.5
    else:
        ## reads are on different chromosomes
        return 1.0
