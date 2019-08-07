#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
import tenkit.safe_json
import evenness_ab
import json
import pandas as pd
import martian
import longranger.cnv.contig_manager as contig_manager
import longranger.cnv.coverage_matrix as coverage_matrix

__MRO__ = """
stage REPORT_AB(
    in  h5 copy_number_profiles "",
    in  json clusters "",
    out json evenness_ab "",
    src py     "stages/repoter_report_ab",
) 
"""

#...............................................................................
def split(args):
    raise Exception("Split is unimplemented")

#...............................................................................
def main(args, outs):
    mask = []
    raw_profiles = []

    raw_profiles, mask = coverage_matrix.load_matrix(args.coverage_profile, args.reference_path);

    clusters = json.load(open(args.clusters, 'r'))
    
    results = evenness_ab.evaluate_ab(raw_profiles, clusters[0], mask, n_trials=200, seed=0, n_merge=5)


    crmap = {}

    contigs = contig_manager.contig_manager(args.reference_path)
    # IS THIS RIGHT? Maybe look at hdf file???
    chroms = contigs.primary_contigs()
    for i in range(0,len(chroms)):
        crmap[chroms[i]] = results[i]

    out_file = open(outs.evenness_ab, 'w')
    out_file.write(tenkit.safe_json.safe_jsonify(crmap))
    out_file.close()
# main

#...............................................................................
def join(args, outs, chunk_defs, chunk_outs):
    raise Exception("Join is unimplemented")
# join

