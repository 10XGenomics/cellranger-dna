#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
from longranger.cnv import coverage_matrix, contig_manager
from crdna.utils import load_genome_data, load_gc_data
from crdna.gc_normalizer import estimate_gc_normalization
import json
import numpy as np

__MRO__ = """
# Estimate linear and quadratic GC bias coefficients and scale
stage ESTIMATE_NORMALIZATION_PARAMETERS(
    in h5       raw_profiles,
    in h5       tracks,
    in string   reference_path,
    #
    out json    gc_norm_params,
    #
    src py      "stages/copy_number_processor/estimate_normalization_parameters",
)
"""

def split(args):
    mat_size_gb = coverage_matrix.get_genome_matrix_size_gb(args.raw_profiles)
    mem_gb = int(np.ceil(3*mat_size_gb + 1))
    return {'chunks': [], 'join': {'__mem_gb' : mem_gb}}

def join(args, outs, chunk_defs, chunk_outs):
    ref = contig_manager.contig_manager(args.reference_path)
    ## only run normalization for single species samples
    species_list = ref.list_species()
    if len(species_list) == 1:
        chroms = ref.primary_contigs(allow_sex_chromosomes=True)
        profiles, mask, _ = load_genome_data(args.raw_profiles, args.tracks, chroms)
        gc = load_gc_data(args.tracks, chroms)
        scale, linear, quadratic = estimate_gc_normalization(profiles, gc, mask)
    else:
        ncells = coverage_matrix.get_num_cells(args.raw_profiles, 
            args.reference_path)
        scale     = [1.0]*ncells
        linear    = [0.0]*ncells
        quadratic = [0.0]*ncells
    with open(outs.gc_norm_params, "w") as out:
        gc_norm_data = {"scale": scale, "linear": linear, 
            "quadratic": quadratic}
        json.dump(gc_norm_data, out, indent=4)
