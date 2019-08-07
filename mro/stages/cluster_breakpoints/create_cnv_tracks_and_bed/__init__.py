#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import os
import pandas as pd
import numpy as np
import itertools as it

import crdna.constants
from longranger.cnv import contig_manager
import martian
from crdna.utils import load_h5

__MRO__ = '''
stage CREATE_CNV_TRACKS_AND_BED
(
    in  string reference_path,
    in  h5     cluster_data,
    in  h5     tracks,
    in  string sex,
    out h5     cnv_tracks,
    out bed    cnv_calls,
    out bed    unmerged_cnv_calls,
    src py     "stages/cluster_breakpoints/create_cnv_tracks_and_bed",
) split using ()   
'''

#
# This is a special value in the ploidy matrix signifying that no ploidy has been assigned
#
MISSING_VALUE = np.iinfo(np.int8).min
#
# This is a special value in the ploidy matrix signifying that a zero ploidy state has been imputed
#
IMPUTED_ZERO = MISSING_VALUE + 1
#
# Only impute ploidy values for segments less than or equal to MAX_IMPUTED_SEGMENT_SIZE
#
MAX_IMPUTED_SEGMENT_SIZE = 500000

def split(args):
    constants = load_h5(args.cluster_data, "constants").to_dict()
    matsize_gb = float(constants["ncells"]*constants["genomebins"])/1e9
    return {'chunks': [], 'join': {'__mem_gb' : int(np.ceil(4*matsize_gb + 1))}}

def join(args, outs, chunk_defs, chunk_outs): 
    args.coerce_strings()
    outs.coerce_strings()

    ref = contig_manager.contig_manager( args.reference_path )
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)

    store = pd.HDFStore( args.cluster_data, "r" )
    windows = store["windows"]
    Q = store["quantized"].values
    #
    # due to the int8 conversion and the use of -127 as a special value,
    # unpredictable bad things will happen if Q>126
    # in practice we don't expect such a thing to ever occur
    #
    martian.log_info( "Found %d bins with Q>126" % (Q>126).sum() )
    Q[ Q>126 ] = 126
    constants = store["constants"]
    store.close( )

    ## cnv track
    ncells = Q.shape[0]

    store = pd.HDFStore( args.tracks, "r" )

    window_size = store["constants"]["window_size"]

    nbins = 0
    for chrom in chroms:
        nbins += store["/map/"+chrom].shape[0]

    C = MISSING_VALUE * np.ones((ncells, nbins), dtype="int8")

    #
    # The C array is filled out with the following convention:
    # * unmasked bins have positive ploidies
    # * masked bins with imputed ploidies are recorded with negative ploidies
    #
    chrom_start = 0
    masked_chrom_start = 0
    chrom_bdy = {}
    for chrom in chroms:
        ctrack = store["/map/"+chrom].values
        cmask = ctrack > crdna.constants.MAPPABILITY_THRESHOLD
        chrom_end          = chrom_start + len(cmask)
        masked_chrom_end   = masked_chrom_start + cmask.sum( ) 
        chrom_bdy[chrom] = (chrom_start, chrom_end)
        C[:, chrom_start:chrom_end][:, cmask] = Q[:, masked_chrom_start:masked_chrom_end]  
        impute_ploidies_for_chromosome_nocall_boundaries(C,chrom_start,chrom_end,window_size)
        chrom_start        = chrom_end
        masked_chrom_start = masked_chrom_end

    store.close( )
    in_store = pd.HDFStore( args.cluster_data, "r" )
    out_store = pd.HDFStore( outs.cnv_tracks, "w" )
    out_store["/cnv_tracks"] = pd.DataFrame(C)
    out_store["/windows"] = windows
    out_store["constants"] = constants
    out_store["/ploidy_conf"] = in_store["/ploidy_conf"]
    out_store["/reads_per_bin"] = in_store["/reads_per_bin"]
    out_store["/scale_factor"] = in_store["/scale_factor"]
    out_store.close( )
    in_store.close( )

    ## break up profile into segments and write to BED

    with open(outs.cnv_calls,"w") as out_bed, open(outs.unmerged_cnv_calls,"w") as out_unmerged_bed:
        for cell in xrange(ncells):
            for chrom in chroms:
                chrom_start, chrom_end = chrom_bdy[chrom]
                ## chrom piece of CNV
                chrom_piece = C[cell,chrom_start:chrom_end]

                for b in get_event_blocks_v2(cell,chrom,chrom_piece,window_size,ref):
                    out_bed.write("\t".join(map(str, b)) + os.linesep)
                for b in get_event_blocks_v2(cell,chrom,chrom_piece,window_size,ref,merge_imputed_blocks=False):
                    out_unmerged_bed.write("\t".join(map(str, b)) + os.linesep)

def get_event_blocks_v2(cell,chrom,chrom_piece,window_size,ref,merge_imputed_blocks=True):
    """Traverse chrom_piece, emitting BED records for
    every contiguous piece with the same ploidy.

    If merge_imputed_blocks is True then identical ploidy blocks will be merged in the output,
    including blocks with imputed ploidy.
    
    This method produces identical results to _old but is significantly faster (2x)
    I also tried speeding everything up by processing all cells at once with the np.where/np.diff
    operations, but that didn't have an appreciable effect.
    """
    if merge_imputed_blocks:
        # convert imputed zeroes back to real zeros
        # don't mutate the original data
        cp2 = chrom_piece.copy()
        cp2[ cp2==IMPUTED_ZERO ] = 0
        breaks = np.where(np.diff(np.abs(cp2))!=0)[0]+1
    else:
        breaks = np.where(np.diff(chrom_piece)!=0)[0]+1
    i = 0
    if len(breaks)==0 or breaks[-1]<len(chrom_piece):
        breaks = it.chain(breaks,[len(chrom_piece)])
    for j in breaks:
        if chrom_piece[i]==MISSING_VALUE:
            i = j
            continue
        #
        # n.b. converting to int here is necessary to convert -128 into 128
        # else chrom_piece[i] is an int8 type and abs(-128)=-128 in int8 land
        #
        ploidy = abs(int(chrom_piece[i]) if chrom_piece[i]!=IMPUTED_ZERO else 0)
        if not merge_imputed_blocks and chrom_piece[i]<0:
            conf = 0.0
        else:
            conf = 1.0
        #
        # BED file interval semantics are zero-based, half-open
        #
        b = [chrom, i*window_size,
             min(j*window_size, ref.contig_lengths[chrom]),
             cell, ploidy, conf]
        yield b
        i = j

def impute_ploidies_for_chromosome_nocall_boundaries(C,chrom_start,chrom_end,window_size):
    """This imputation strategy imputes ploidies across masked regions
    by using the called ploidies on the boundaries of the masked segment
    only if the ploidies agree.  If they disagree then the region is nocalled
    (left as MISSING_VALUE)."""
    if C.shape[0] == 0: return
    # half-open interval for current masked segment
    start, end = -1,-1
    in_missing_value = False
    max_bins_for_imputation = MAX_IMPUTED_SEGMENT_SIZE / window_size
    for i in xrange(chrom_start,chrom_end):
        if not in_missing_value and C[0,i]==MISSING_VALUE:
            start = i
            in_missing_value = True
        if in_missing_value and C[0,i]!=MISSING_VALUE:
            in_missing_value = False
            end = i
            # process current masked segment
            if start>chrom_start:
                left_ploidies = C[:,start-1]
            else:
                left_ploidies = None
            right_ploidies = C[:,end]
            # if we're not at the left edge, impute for internal masked segments
            if left_ploidies is not None and end-start<=max_bins_for_imputation:
                impute_indices = (left_ploidies==right_ploidies)
                C[impute_indices,start:end] = np.where(left_ploidies[impute_indices]>0,
                    -left_ploidies[impute_indices], IMPUTED_ZERO)[:,np.newaxis]
    # don't need to process final missing segment if it exists
    # it's left as missing value

def impute_ploidies_for_chromosome_halfway(C,chrom_start,chrom_end):
    """This imputation strategy imputes ploidies across masked regions
    by using the called ploidies on the boundaries of the masked segment
    and imputing calls halfway in from each side."""
    #
    # impute CNVs over masked regions by lifting over the ploidy call to the left
    #
    first_bin = None
    for index in xrange(chrom_start+1,chrom_end):
        if first_bin is None and C[0, index] != MISSING_VALUE:
            first_bin = index
        if C[0, index] == MISSING_VALUE:
            if C[0, index-1] >= 0:
                C[:, index] = np.where( C[:,index-1]>0, -C[:, index-1], IMPUTED_ZERO )
            else:
                #
                # lift over imputed ploidy
                #
                C[:, index] = C[:, index-1]
    #
    # impute CNVs at the beginning by lifting over the ploidy call to the right
    #
    if first_bin > 0:
        C[:, chrom_start:first_bin] = np.where( C[:,first_bin]>0, -C[:, first_bin], IMPUTED_ZERO )[:, np.newaxis]
    #
    # now improve upon the imputed CNVs by working from the right and filling in each
    # call halfway
    #
    i = chrom_end - 1
    while i>=0:
        if C[0,i]<0 and i<chrom_end-1:
            right_ploidy = C[:,i+1]
            #
            # find the extent of this masked segment
            #
            j = 1
            while i-j>=0 and C[0,i-j]<0:
                j += 1
            #
            # now [i-j+1,i] is the extent of the imputed segment
            # fill [i-j/2+1,i] with right_ploidy
            #
            half = j/2
            if half>0:
                C[:, i-half+1:i+1 ] = right_ploidy[:,np.newaxis]
            i = i - j
        i = i - 1
