#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Summarize genome coverage from BAM file
#
import martian
import os
import tenkit.bam as tk_bam
import tenkit.stats as tk_stats
import pandas as pd
import numpy as np
import crdna.read_filter
from crdna.constants import PROCESSED_BARCODE_TAG, MAPPABILITY_THRESHOLD
from longranger.cnv import contig_manager

__MRO__ = """
stage MAKE_PROFILES(
    in  bam      possorted,
    in  string   reference_path,
    in  h5       tracks,
    in  int      window_size,
    in  map      cell_barcodes,
    out h5       profiles,
    src py       "stages/reporter/make_profiles",
) split using (
    in  string[] chroms,
)
"""

#...............................................................................
## split up bam by chromosomes
def split(args):
    ref = contig_manager.contig_manager(args.reference_path)
    ## every primary chromosome gets its own chunk
    ## all the secondary pieces are in one chunk
    chrom_chunks = []
    non_primary_chunk = []
    for chrom in ref.list_all_contigs():
        if ref.is_primary_contig(chrom, allow_sex_chromosomes=True):
            chrom_chunks.append([chrom])
        else:
            non_primary_chunk.append(chrom)
    
    if len(non_primary_chunk) > 0:
        chrom_chunks.append( non_primary_chunk )
    
    chrom_sizes = ref.contig_lengths
    max_size = 0
    for chroms in chrom_chunks:
        chunk_size = sum([chrom_sizes[chrom] for chrom in chroms])
        max_size = max(max_size, chunk_size)
    
    nbcs = 0
    for v in args.cell_barcodes.itervalues():
        nbcs += len(v)

    max_mat_size = 4*nbcs*max_size/args.window_size
    chunk_mem_gb = int(np.ceil((1.0*max_mat_size/1e9) + 1))
    join_mem_gb = int(np.ceil(1.0*max_mat_size/1e9 +
                              1.0*sum(chrom_sizes.values())/args.window_size/1e9 + 1))

    chunk_defs = [{'chroms': chroms, '__mem_gb': chunk_mem_gb} for chroms in chrom_chunks]
    return {'chunks': chunk_defs, 'join': {'__mem_gb': join_mem_gb}}
# split

#...............................................................................
## process each chunk
def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    ref = contig_manager.contig_manager(args.reference_path)

    bam_in = tk_bam.create_bam_infile(args.possorted)
    assert bam_in.has_index()
    
    chrom_size = dict(zip(bam_in.references, bam_in.lengths))
    
    ## put cell barcodes in canonical order and index them
    bc_list = []
    for v in args.cell_barcodes.itervalues():
        bc_list.extend(v.keys())
    # for v
    bc_list = list(set(bc_list))
    bc_list.sort( )
    
    nbcs = len(bc_list)
    bc_index = {}
    for bci in xrange(nbcs):
        bc_index[bc_list[bci]] = bci
    # for bci
    
    genomebins = 0
    in_track = pd.HDFStore(args.tracks, "r")
    store = pd.HDFStore(outs.profiles, "w")
    for chrom in args.chroms:
        chrom_length = chrom_size[chrom]
        nbins = chrom_length/args.window_size + int(chrom_length % args.window_size!=0)
        genomebins += nbins
        counts = np.zeros( (nbcs, nbins), dtype="float32" )
        for rec in bam_in.fetch( str(chrom) ):
            if not rec.has_tag(PROCESSED_BARCODE_TAG):
                continue
            # if rec
            
            inserts = crdna.read_filter.inserts_per_alignment(rec)
            if inserts < 1e-6:
                continue
            # if inserts

            ## if this is a non cell barcode, then skip
            bc = rec.get_tag(PROCESSED_BARCODE_TAG)
            bci = bc_index.get( bc, None)
            if bci is None:
                continue    
            # if bci
            
            w = rec.reference_start/args.window_size
            counts[bci, w] += inserts
        # for rec

        store["/contigs/" + chrom] = pd.DataFrame( counts )
        ## mask is defined based on the mappability information alone
        mask = pd.Series(in_track["/map/"+chrom] > MAPPABILITY_THRESHOLD)
        if ref.is_primary_contig(chrom) and tk_stats.robust_divide(sum(mask), len(mask)) == 0:
            martian.exit('Chromosome {} has no mappable bins.'.format(chrom))
        store["/masks/" + chrom] = mask

    # for chrom

    ## store the window size in the h5
    store["constants"] = pd.Series({"window_size": args.window_size,
        "ncells": nbcs, "genomebins": genomebins})
    in_track.close( )
    store.close( )
# main

#...............................................................................
## join chunks
def join(args, outs, chunk_defs, chunk_outs):
    store = pd.HDFStore( outs.profiles, "w" )
    
    ## put cell barcodes in canonical order
    bc_list = []
    for v in args.cell_barcodes.itervalues():
        bc_list.extend(v.keys())
    # for v
    bc_list = list(set(bc_list))
    bc_list.sort( )
    
    store["barcodes"] = pd.Series( bc_list )

    # load ref to determine primary-ness

    all_chroms = []
    masks = []
    genomebins = 0
    ncells = None
    for chunk_out, chunk_def in zip(chunk_outs, chunk_defs):
        chroms = chunk_def.chroms
        all_chroms.extend(chroms)
        profile_chunk = pd.HDFStore( chunk_out.profiles, "r" )
        for chrom in chroms:
            mask = profile_chunk["/masks/" + chrom]
            masks.extend(mask)
            store["/contigs/" + chrom] = profile_chunk["/contigs/" + chrom]
            store["/masks/" + chrom] = mask
        genomebins += profile_chunk["constants"]["genomebins"]
        if ncells is None:
            ncells = profile_chunk["constants"]["ncells"]

        # for chrom
        profile_chunk.close( )
    # for chunk_out, chunk_def
    ## store the window size in the h5
    store["constants"] = pd.Series({"window_size": args.window_size,
        "ncells" : ncells, "genomebins" : genomebins})

    ref = contig_manager.contig_manager(args.reference_path)
    write_mask_bed(outs.mappable_regions,store,all_chroms,args.window_size,ref,
        args)
    
    store.close( )
# join

def write_mask_bed(bedfile,store,chroms,window_size,ref, args):
    """Write a BED file corresponding to the mask=True regions
    for our profiles."""
    chroms = sorted(chroms)
    with open(bedfile,'w') as outfile:
        version = martian.get_pipelines_version()
        outfile.write("#cellranger-dna {}\n".format(version))
        outfile.write("#reference genome: {}\n".format(args.reference_path))
        outfile.write("#chrom\tstart\tend\n")

        for chrom in chroms:
            chrom_length = ref.contig_lengths[chrom]
            mask = store['/masks/'+chrom]
            start,end = None,None
            in_mask = False
            for i in xrange(len(mask)):
                if not in_mask and mask[i]:
                    start = i
                    in_mask = True
                if in_mask and not mask[i]:
                    end = i
                    in_mask = False
                    outfile.write( '\t'.join(str(s) for s in [\
                        chrom,
                        start*window_size,
                        min(end*window_size,chrom_length),
                        ]) + os.linesep
                    )
            if in_mask and \
               (start is not None) and \
               (end is None or end*window_size<chrom_length):
                outfile.write( '\t'.join(str(s) for s in [\
                    chrom,
                    start*window_size,
                    chrom_length,
                    ]) + os.linesep
                )
