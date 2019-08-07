#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
# 
# This is a no op module that has the same interface as the normalization module
# except that it does nothing
#

import longranger.cnv.coverage_matrix as coverage_matrix
import pandas as pd

__MRO__ = """
stage DO_NOTHING_ABOUT_GC_BIAS(
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

def main(args, outs):
    raw_profiles, mask = coverage_matrix.load_matrix(
        args.raw_singlecell_profiles, args.reference_path)
    chromosomes = coverage_matrix.list_primary_contigs(
        args.raw_singlecell_profiles, args.reference_path)
    print(chromosomes)

    bin_size = coverage_matrix.get_bin_size(args.raw_singlecell_profiles)
    tracks = pd.HDFStore(args.tracks, 'r')
    coverage_matrix.store_matrix(
        file_name=outs.normalized_singlecell_profiles,
        chroms=chromosomes,
        profiles=raw_profiles,
        tracks=tracks,
        window_size=bin_size)
    tracks.close()

