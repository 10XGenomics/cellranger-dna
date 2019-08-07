#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
import longranger.cnv.coverage_matrix as coverage_matrix
import martian
from crdna.singlecell_dna_cnv.gc_bias_estimator import estimate_gc_bias

__MRO__ = """
stage ESTIMATE_GC_BIAS_COEFFICIENTS(
    in h5       raw_singlecell_profiles,
    in h5       tracks,
    in string   reference_path,
    #
    out float   linear,
    out float   quadratic,
    #
    src py      "stages/copy_number_processor/estimate_gc_bias_coefficients",
)
"""

#...............................................................................
def split(args):
    raise Exception("Split is unimplemented")
# split

#...............................................................................
def main(args, outs):
    raw_profiles, mask = coverage_matrix.load_matrix(args.raw_singlecell_profiles, args.reference_path)
    # Get mappability, GC content
    ncells = raw_profiles[0].shape[0]
    # Default:
    (intercept, linear, quadratic) = (1.0, 0.0, 0.0)
    # Sum up all single-cell profiles
    try:
        print('DEBUG 0')
        result = estimate_gc_bias(args.raw_singlecell_profiles, args.tracks, args.reference_path)
        print('DEBUG result')
        print(result)
        (quadratic, linear, intercept) = result['Summary']['quadratic_coefficients']
        print('DEBUG intercept=%f, linear=%f, quadratic=%f' % (intercept, linear, quadratic))
    except Exception as error:
        martian.alarm("stages/copy_number_processor/estimate_gc_bias_coefficients/__init__ encountered an exception. Error: %s" % repr(error))
    # try/except
    #
    # Export scale factor and GC bias coefficients
    outs.linear = linear
    outs.quadratic = quadratic
# main

#...............................................................................
def join(args, outs, chunk_defs, chunk_outs):
    raise Exception("Join is unimplemented")
# join

