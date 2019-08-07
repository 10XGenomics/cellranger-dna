#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import json
from longranger.cnv import contig_manager
from longranger.cnv import coverage_matrix
import subprocess
import pandas as pd
import shutil
import numpy as np
from crdna.singlecell_dna_cnv import differential_variance_ratio as dvr
import os

__MRO__ = """
# Inputs: normalized, scaled, clustered profiles at high resolution - 20kb
stage CALL_CNVS(
    in  string  reference_path,
    in  h5      coverage_profile,
    in  map     params,
    in  json    clusters,
    #   
    out bed     cnvs,
    #   
    src py      "stages/copy_number_processor/cnv_caller",
) split using (
    in  int[]   cluster_indices,
    in  string  chroms,
)
"""

################################################################################
## Split on the cartesian product of contigs and clusters
def split(args):
    args.coerce_strings( )
    #
    ctg_mgr = contig_manager.contig_manager(args.reference_path)
    chroms  = ctg_mgr.primary_contigs(allow_sex_chromosomes=False)
    #
    # Handle case when clusters = None
    if args.clusters is None:
        ncells = coverage_matrix.get_num_cells(args.coverage_profile, args.reference_path)
        clusters = [[x] for x in xrange(ncells)]
    else:
        f = open(args.clusters)
        clusters = json.load(f)
    cart_prod = []
    for chrom in chroms:
        for ci, cluster in enumerate(clusters):
            chunk_def = {'chrom': chrom, 'cluster_index': ci}
            cart_prod.append(chunk_def)
        # for cluster
    # for chrom
    #
    # Split these pieces into at most MAX_CHUNKS chunks
    MAX_CHUNKS = 100
    npieces = len(cart_prod)
    pieces_per_chunk = npieces / MAX_CHUNKS + int(npieces % MAX_CHUNKS != 0)
    chunks = []
    start = 0
    while start < npieces:
        chunk_def = {"chroms": [], "cluster_indices": []}
        end = min(start+pieces_per_chunk, npieces)
        for i in xrange(start, end):
            chunk_def["chroms"].append(cart_prod[i]["chrom"])
            chunk_def["cluster_indices"].append(cart_prod[i]["cluster_index"])
        chunks.append(chunk_def)
        start += pieces_per_chunk
    assert len(chunks) <= MAX_CHUNKS
    return {'chunks': chunks}
# split

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
    # Iterate over all chrom,cluster_index pairs in chunk
    for chrom, cluster_index in zip(args.chroms, args.cluster_indices):
        scaled_cluster_profiles = pd.HDFStore(args.coverage_profile, 'r')
        cnv_profile = scaled_cluster_profiles['/contigs/' + chrom].astype(np.float64).values[cluster_index, :]
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
        # per bin at each cell.
        #
        cnv_segments = dvr.partition_profile(
            cnv_profile, segment_length=dvr_segment_length, 
            f_cutoff=dvr_f_cutoff, trim_level=dvr_trim_level, 
            min_length=dvr_min_length)
        #
        # TODO: recover original bin indices (preceding removal of NAs)
        #

        # try/except
        #
        # Now write that array out to a BED file; coalescing adjecent bins with the
        # same copy number
        # Fetch bin_size from profile.h5
        # bin_size = all_profiles["constants"]["window_size"]
        bin_size = 20000
        # export_segments(outs.cnvs, cnv_segments, bin_size, chrom, barcodes) # TODO: export cluster ID 
        export_segments(outs.cnvs, cnv_segments, bin_size, chrom, cluster_index)
    # for chrom, cluster_index
        print('Leaving __init__.main()')
    print('.' * 80)
# main

################################################################################
## Join step.  Each chunk returns a bed file. This just concatenates
## the bed file from each chunk and sorts them.
def join(args, outs, chunk_defs, chunk_outs):
    tmpfile = 'tmp.bed'
    print tmpfile
    outf = open(tmpfile, 'w')
    #
    for chunk in chunk_outs:
        if not os.path.exists(chunk.cnvs):
            continue
        # if !exists
        inf = open(chunk.cnvs)
        shutil.copyfileobj(inf, outf, 1024 * 1024)
        inf.close()
    # for chunk
    outf.close()
    #
    final_output = open(outs.cnvs, 'w')
    subprocess.check_call(
        ['bedtools', 'sort', '-i', tmpfile], stdout=final_output)
# join

################################################################################
##
def export_segments(file_name, segments, resolution, chrom, cluster_index):
    output_file = open(file_name, 'a')
    n_segments = segments.shape[0]
    confidence = 1.0 # TODO: replace this stub with actual confidence score
    #diploid = 2
    for segment in range(n_segments):
        ploidy = int(round(segments.loc[segment, 'CopyNumber'] * 2))
        if True: #ploidy != diploid:
            start = segments.loc[segment, 'Start'] * resolution + 1
            end = (segments.loc[segment, 'End'] + 1) * resolution
            out_text = '%s\t%i\t%i\t%i\t%f\t%s\n' % (
                chrom, start, end, ploidy, confidence, str(cluster_index))
            output_file.write(out_text)
        # if ploidy
    # for segment
    output_file.close()
# export_segments

