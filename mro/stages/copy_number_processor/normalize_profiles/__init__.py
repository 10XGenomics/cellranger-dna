#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
from longranger.cnv import coverage_matrix, contig_manager
import json
from crdna.gc_normalizer import gc_normalize
from crdna.utils import load_h5

import pandas as pd
import numpy as np

__MRO__ = """
#
Normalize GC bias using GC bias parameters from the previous stage
stage NORMALIZE_PROFILES(
    in h5       raw_profiles,
    in h5       tracks,
    in string   reference_path,
    in json     gc_norm_params,
    #
    out h5      normalized_profiles,
    #
    src py      "stages/copy_number_processor/normalize_profiles",
) split using ()
"""


def split(args):
    mat_size_gb = coverage_matrix.get_genome_matrix_size_gb(args.raw_profiles)
    mem_gb = int(np.ceil(3*mat_size_gb + 1))
    return {'chunks': [], 'join': {'__mem_gb' : mem_gb}}


def load_data(profile_name, tracks_name, chroms):
    profile_store = pd.HDFStore(profile_name, "r")
    profiles = []
    mask = []
    for chrom in chroms:
        profiles.append(profile_store["/contigs/" + chrom].values)
        mask.append(profile_store["/masks/" + chrom].values) 
    profile_store.close()

    tracks_store = pd.HDFStore(tracks_name, "r")
    gc = []
    for chrom in chroms:
        gc.append(tracks_store["/GC/" + chrom].values)
    tracks_store.close()

    return profiles, gc, mask


def join(args, outs, chunk_defs, chunk_outs):
    ref = contig_manager.contig_manager(args.reference_path)
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)

    profiles, gc, mask = load_data(args.raw_profiles, args.tracks, chroms)
    gc_norm_params = json.load(open(args.gc_norm_params, "r"))
    scale = gc_norm_params["scale"]
    linear = gc_norm_params["linear"]
    quadratic = gc_norm_params["quadratic"]

    norm_profiles = gc_normalize(profiles, gc, linear, quadratic, chroms)

    bin_size = coverage_matrix.get_bin_size(args.raw_profiles)

    coverage_matrix.store_matrix(file_name=outs.normalized_profiles,
        chroms=chroms, profiles=norm_profiles, tracks=None,
        window_size=bin_size, masks=mask, dtype="float32")

    store = pd.HDFStore(outs.normalized_profiles, "a")
    constants = load_h5(args.raw_profiles, "constants")
    store["constants"] = constants
    store.close()

    store = pd.HDFStore(outs.normalized_profiles, "a")
    store["/gc_params/scale"]     = pd.Series(scale)
    store["/gc_params/linear"]    = pd.Series(linear)
    store["/gc_params/quadratic"] = pd.Series(quadratic)
    store.close()
