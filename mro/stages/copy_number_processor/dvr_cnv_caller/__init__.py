#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
#import json
#from longranger.cnv import contig_manager
from longranger.cnv import coverage_matrix
#import subprocess
import pandas as pd
#import shutil
import numpy as np
from crdna.singlecell_dna_cnv import differential_variance_ratio as dvr
import crdna.singlecell_dna_cnv.single_cell_normalizer as scn
#import os

__MRO__ = """
stage DETECT_CNVS(
    in  h5      raw_profiles,
    in  h5      normalized_profiles,
    in  string  reference_path,
    in  string  sex,
    in  map     params,
    #
    out bed     cnvs,
    #
    src py      "stages/copy_number_processor/dvr_cnv_caller",
)        
"""

################################################################################
## Unimplemented
#def split(args):
#    pass
## split

################################################################################
## Compute a CNV calls from the profile for a specific choice of cluster
## and contig.
def main(args, outs):
    # Read in heuristics
    dvr_segment_length = int(args.params["dvr_segment_length"])
    dvr_f_cutoff       = args.params["dvr_f_cutoff"]
    dvr_trim_level     = args.params["dvr_trim_level"]
    dvr_min_length     = int(args.params["dvr_min_length"])
    #
    print('-' * 80)
    print('Entering __init__.main()')
    #
    default_n_merge = 50
    default_bin_size = 2e4
    target_bin_count = 200.0
    confident_genome_fraction = 1.0
    raw_profiles, _ = coverage_matrix.load_matrix(
        args.raw_profiles, args.reference_path)
    #
    # Iterate over all chrom_name 
    chromosomes = coverage_matrix.list_primary_contigs(
        args.normalized_profiles, args.reference_path, 
        allow_sex_chromosomes=False) # TODO: allow sex chromosomes
    normalized_profiles_h5 = pd.HDFStore(args.normalized_profiles, 'r')
    try:
        original_bin_size = normalized_profiles_h5["constants"]["window_size"]
    except Exception as error:
        print('__init__.main() caught an error %s' % repr(error))
        original_bin_size = default_bin_size
    # try/except
    #
    segments = []
    profiles_list = []
    for chrom_name in chromosomes:
        tmp = []
        chr_profiles = normalized_profiles_h5[
            '/contigs/' + chrom_name].astype(np.float64).values
        mask = coverage_matrix.contig_to_mask(
            normalized_profiles_h5, chrom_name)
        #
        #print(mask[:10])
        assert(mask.shape[0] == chr_profiles.shape[1])
        n_nodes = chr_profiles.shape[0]
        for node_id in xrange(n_nodes):
            # Dynamically determine bin resolution
            raw_cell_counts = scn.get_single_cell_counts(raw_profiles, node_id)
            n_bins = raw_cell_counts.shape[0]
            n_merge = scn.estimate_n_merge(raw_cell_counts, target_bin_count, confident_genome_fraction)
            if (np.isnan(n_merge) | 
                ~np.isfinite(n_merge) | 
                (n_merge > n_bins)):
                n_merge = default_n_merge
            # if NaN
            bin_size = original_bin_size * n_merge
            #print('n_merge=%d' % n_merge)
            #print('bin_size=%d' % bin_size)
            #
            weights = mask.copy()
            weights = weights.astype(float)
            weights = scn.merge_bins_single_cell_single_chrom(
                weights, n_merge, average=True)
            weights[weights <= 0.2] = np.nan
            #print('weights:')
            #print(weights)
            #
            cnv_profile = chr_profiles[node_id, :]
            #print('cnv_profile before masking:')
            #print(cnv_profile[:20])
            cnv_profile[~mask] = 0.0 
            #print('cnv_profile after masking:')
            #print(cnv_profile[:20])
            #
            # Merge bins:
            cnv_profile = scn.merge_bins_single_cell_single_chrom(
                cnv_profile, n_merge, average=True)
            cnv_profile /= weights
            #
            #print('type(cnv_profile):')
            #print(type(cnv_profile))
            #print('cnv_profile.shape:')
            #print(cnv_profile.shape)
            #print('cnv_profile:')
            #print(cnv_profile)
            #print(cnv_profile.tolist())
            #
            # Call CNVs on this profile. The result is an array with the copy number
            # per bin at each node_id.
            #
            block = dvr.partition_profile(
                cnv_profile, segment_length=dvr_segment_length, 
                f_cutoff=dvr_f_cutoff, trim_level=dvr_trim_level, 
                min_length=dvr_min_length)
            block['NodeID'] = node_id
            block['Chr'] = chrom_name
            block['Start'] = block['Start'] * bin_size + 1
            block['End'] = (block['End'] + 1) * bin_size
            segments.append(block)
            #
            # Debugging:
            tmp.append(cnv_profile)
            #
        # for node_id
        #
        # Debugging:
        profiles_list.append(tmp)
        #
    # for chrom_name
    normalized_profiles_h5.close()
    #
    #print('block:')
    #print(block)
    export_segments(outs.cnvs, segments)
    # 
    # Debugging:
    np.save(outs.cnvs + '_profiles.npy', profiles_list)
    print('Leaving __init__.main()')
    print('.' * 80)
# main

################################################################################
## Unimplemented
#def join(args, outs, chunk_defs, chunk_outs):
#    pass
## join

################################################################################
##
def export_segments(file_name, segments):
    output_file = open(file_name, 'a')
    confidence = 1.0 # TODO: replace this stub with actual confidence score
    #n_blocks = len(segments)
    #diploid = 2
    for block in segments:
        n_segments = block.shape[0]
        for segment in range(n_segments):
            ploidy = block.loc[segment, 'CopyNumber']
            chrom_name = block.loc[segment, 'Chr']
            node_id = str(block.loc[segment, 'NodeID'])
            if True: #ploidy != diploid:
                start = block.loc[segment, 'Start']
                end = block.loc[segment, 'End']
                try:
                    out_text = '%s\t%i\t%i\t%i\t%f\t%s\n' % (
                        chrom_name, start, end, ploidy, confidence, node_id)
                    output_file.write(out_text)
                except Exception as error:
                    print('export_segments() caught an error: %s' % repr(error))
                    print('chrom_name, start, end, ploidy, confidence, node_id:')
                    print(chrom_name, start, end, ploidy, confidence, node_id)
                # try/except
            # if ploidy
        # for segment
    # for block
    output_file.close()
# export_segments

