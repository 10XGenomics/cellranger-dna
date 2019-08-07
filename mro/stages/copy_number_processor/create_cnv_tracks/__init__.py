#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Convert a BED file into an h5 genome track
#

import pandas as pd
#from crdna.utils import create_chunks
import numpy as np
from longranger.cnv import contig_manager

__MRO__ = """
stage CREATE_CNV_TRACKS(
    in  bed    cnv_calls,
    in  string reference_path,
    in  string sex,
    in  int    window_size,
    out h5     cnv_tracks,
    src py     "stages/copy_number_processor/create_cnv_tracks",
)
"""

# TODO finish parallelization
#def split( args ):
#    indices = set([])
#    for line in open( args.cnv_calls, "r"):
#        index = int(line.strip( ).split( )[-1])
#        indices.add(index)
#
#    MAX_CHUNKS = 100
#    all_indices = range(len(indices))
#    chunked_indices = create_chunks(all_indices, MAX_CHUNKS)
# 
#    chunks = []
#    for indices in chunked_indices:
#        chunk_def = {"indices": indices, "__mem_gb" : 6}
#        chunks.append(chunk_def)

def main(args, outs):
    ref = contig_manager.contig_manager( args.reference_path )
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)

    ## read in calls as dataframe
    calls = pd.read_csv( args.cnv_calls, sep="\t", names = ["chrom", "start", "end",
        "ploidy", "confidence", "cluster_index"] )

    ## figure out dimensions of cnv_tracks and gather chrom data
    nclusters = len(np.unique(calls["cluster_index"].values))
    window_size = args.window_size
    contig_sizes = ref.get_contig_lengths( )
    chrom_bin_sizes = {}
    nbins = 0
    chrom_offset = {}
    for chrom in chroms:
        chrom_offset[chrom] = nbins
        csize = contig_sizes[chrom]
        cbins = csize / window_size + int(csize % window_size != 0)
        chrom_bin_sizes[chrom] = cbins
        nbins += cbins

    cnv_tracks = np.zeros((nclusters, nbins), dtype="int32")
    for chrom in chroms:
        p = ref.expected_ploidy(chrom, args.sex)
        csize = chrom_bin_sizes[chrom]
        cnv_tracks[:, chrom_offset[chrom]:chrom_offset[chrom]+csize] = p
        
    nclusters = calls["cluster_index"].unique( ).shape[0]
    for ci in xrange(nclusters):
        print ci
        cluster_calls = calls[calls["cluster_index"] == ci] 
        for _, row in cluster_calls.iterrows( ):
            offset = chrom_offset[row["chrom"]]
            assert (row["start"]-1) % window_size == 0
            #assert row["end"] % window_size == 0
            sbin = offset + (row["start"]-1)/window_size
            ebin = offset + row["end"]/window_size
            cnv_tracks[ci, sbin:ebin] = row["ploidy"]

    out_store = pd.HDFStore( outs.cnv_tracks, "w" )
    out_store["cnv_tracks"] = pd.DataFrame(cnv_tracks)
    out_store.close( )

