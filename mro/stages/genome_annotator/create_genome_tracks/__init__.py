#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Summarize genome coverage from BAM file
#

import numpy as np
import pandas as pd
from longranger.cnv import contig_manager
import tenkit.reference as tk_ref

__MRO__ = """
stage CREATE_GENOME_TRACKS(
    in  string   reference_path,
    in  int      window_size,
    out h5       tracks,
    src py       "stages/reporter/create_genome_tracks",
) split using (
    in  string[] chroms,
)
"""

## split up bam by chromosomes
def split(args):
    ctg_mgr = contig_manager.contig_manager(args.reference_path)

    ## every primary chromosome gets its own chunk
    ## all the secondary pieces are in one chunk
    chrom_chunks = []
    non_primary_chunk = []
    for chrom in ctg_mgr.list_all_contigs():
        if ctg_mgr.is_primary_contig(chrom, allow_sex_chromosomes=True):
            chrom_chunks.append([chrom])
        else:
            non_primary_chunk.append(chrom)
    if len(non_primary_chunk) > 0:
        chrom_chunks.append( non_primary_chunk )
    
    chunk_defs = [{'chroms': chroms, '__mem_gb': 12} for chroms in chrom_chunks]

    return {'chunks': chunk_defs, 'join': {'__mem_gb': 12}}

## process each chunk
def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    genome = tk_ref.open_reference(args.reference_path)

    store = pd.HDFStore( outs.tracks, "w" )
    for chrom in args.chroms:
        gc_count = []
        N_count  = []
        
        seq = genome[chrom]
        pos = 0
        while pos < len(seq):
            gc = 0
            N  = 0
            total_bases = min(pos+args.window_size, len(seq)) - pos
            for i in xrange(total_bases):
                base = seq[pos + i].upper( )
                gc += int(base == "G" or base == "C")
                N  += int(base == "N")
            gc_denominator = float(max(total_bases-N, 1))
            gc_count.append( float(gc)/gc_denominator )
            N_count.append( float(N)/float(total_bases) )
            pos += args.window_size

        store["/GC/%s"%chrom] = pd.Series(gc_count)
        store["/N/%s"%chrom]  = pd.Series(N_count)
    store.close( )

## join chunks
def join(args, outs, chunk_defs, chunk_outs):
    store = pd.HDFStore( outs.tracks, "w" )
    print "# GC content per bin"
    print "contig, nbins, mean, min, p25, p50, p75, max"
    for chunk_out in chunk_outs:
        in_track = pd.HDFStore( chunk_out.tracks, "r" )
        for key in in_track.keys( ):
            store[key] = in_track[key]
            # logging
            if key[:4] == "/GC/":
                chrom = key.lstrip("/GC/")
                nbins = len(in_track[key])
                mean = np.nanmean(in_track[key])
                vals = np.nanpercentile(in_track[key], [0, 25, 50, 75, 100])
                print ("{}, {}, " + ", ".join(["{:.3f}"]*6)).format(chrom, nbins, mean, *vals)
        in_track.close( )
    store["constants"] = pd.Series({"window_size": args.window_size})
    store.close( )
