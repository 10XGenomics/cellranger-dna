#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
import longranger.cnv.coverage_matrix as coverage_matrix
from crdna.singlecell_dna_cnv import vesna
import martian
import numpy as np
import pandas as pd

__MRO__ = """
stage NORMALIZE_GC_BIAS(
    in h5       raw_singlecell_profiles,
    in h5       tracks,
    in string   reference_path,
    in float    linear,
    in float    quadratic,
    # 
    out h5      normalized_singlecell_profiles, 
    #
    src py      "stages/copy_number_processor/normalize_gc_bias",
)       
"""
# 
ordered_chromosomes = ['hg19_chr' + str(i) for i in range(1, 23) + ['X', 'Y']]
#
#............................................................................
def get_chrom_index(chrom_name, chromosomes):
    chrom_index = chromosomes.index(chrom_name)
    return(chrom_index)
# get_chrom_index
#
#............................................................................
def get_mappability(bin_parameters, chrom_name, chromosomes):
    chrom_index = get_chrom_index(chrom_name, chromosomes)
    mappability = bin_parameters[chrom_index][5, :]
    return(np.array(mappability))
# get_mappability
#
#............................................................................
def get_gc(bin_parameters, gc0, chrom_name, chromosomes):
    chrom_index = get_chrom_index(chrom_name, chromosomes)
    gc = bin_parameters[chrom_index][3, :]
    gc_gc0 = np.array(gc) - gc0
    return(gc_gc0)
# get_gc
#
#...............................................................................
def split(args):
    raise Exception("Split is unimplemented")
# split
#
#...............................................................................
def main(args, outs):
    normalized_profiles = []
    raw_profiles, mask = coverage_matrix.load_matrix(
        args.raw_singlecell_profiles, args.reference_path)
    print('len(mask)=%d' % len(mask))
    print('len(raw_profiles)=%d' % len(raw_profiles))
    
    chromosomes = coverage_matrix.list_primary_contigs(
        args.raw_singlecell_profiles, args.reference_path)
    print('chromosomes:')
    print(chromosomes)
    n_chrom = len(chromosomes)
    #
    # Get mappability, GC content:
    bin_parameters = []
    vesna.load_track_parameters(args.tracks, bin_parameters)
    n_cells = raw_profiles[0].shape[0]
    linear = args.linear
    quadratic = args.quadratic
    gc0 = 0.45 # TODO: Replace this with mean of GC in good bins across entire genome
    #
    remove = []
    for chrom_index, chrom_name in enumerate(chromosomes):
        try:
            mappability =  get_mappability(bin_parameters, chrom_name, ordered_chromosomes)
            gc_gc0 = get_gc(bin_parameters, gc0, chrom_name, ordered_chromosomes)
            print('len(mappability)=%d' % len(mappability))
            print('len(gc_gc0)=%d' % len(gc_gc0))
            print('raw_profiles[chrom_index].shape:')
            print(raw_profiles[chrom_index].shape)
            expectation = mappability * (1.0 + linear * gc_gc0 + quadratic * gc_gc0 * gc_gc0)
            #print('expectation')
            #print(expectation.tolist())
            tmp = np.zeros(raw_profiles[chrom_index].shape, dtype='float')
            for cell in range(n_cells):
                #print('tmp[cell, :] before:')
                #print(tmp[cell, :].tolist())
                tmp[cell, :] = raw_profiles[chrom_index][cell, :] / expectation
                tmp[cell, tmp[cell, :] < 0.0] = 0.0
                #print('tmp[cell, :] after:')
                #print(tmp[cell, :].tolist())
            # for cell
            normalized_profiles.append(tmp)
        except Exception as error:
            martian.alarm("stages/copy_number_processor/normalize_gc_bias/__init__ encountered an exception. Error: %s" % repr(error))
            print("stages/copy_number_processor/normalize_gc_bias/__init__ encountered an exception. Error: %s" % repr(error))
            print('Removing chrom_name=%s, chrom_index=%d (absent from input raw profiles)' % (chrom_name, chrom_index))
            remove.append(chrom_name)
        # try/except
    # for chrom
    for chrom_name in remove:
        if chrom_name in chromosomes:
            chromosomes.remove(chrom_name)
        # if chrom_name
    # for chrom_name
    #
    # Export normalized cell profiles
    bin_size = 20000 # TODO: Fetch this value from input raw_profiles h5 file
    tracks = pd.HDFStore(args.tracks, 'r')
    coverage_matrix.store_matrix(
        file_name=outs.normalized_singlecell_profiles,
        chroms=chromosomes,
        profiles=normalized_profiles,
        tracks=tracks,
        window_size=bin_size)
    tracks.close()
# main
#
#...............................................................................
def join(args, outs, chunk_defs, chunk_outs):
    raise Exception("Join is unimplemented")
# join

