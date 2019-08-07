#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
import tenkit.safe_json
import cluster_jedna
import longranger.cnv.coverage_matrix as coverage_matrix
import martian

__MRO__ = """
stage GENERATE_FINAL_CLUSTERS(
    in  string coverage_profile "",
    out string cell_clusters        "",
    src py     "stages/copy_number_processor/profile_clusterer/generate_final_clusters",
) split using (
    in  string chunk_start,
    in  string chunk_end,
)
"""

#...............................................................................
def split(args):
    raise Exception("Split is unimplemented")

#...............................................................................
def main(args, outs):
    mask = []
    raw_profiles = []
    
    raw_profiles, mask = coverage_matrix.load_matrix(args.coverage_profile, args.reference_path)

    ncells = raw_profiles[0].shape[0]
    try:
        if args.skip_clustering:
            print('Skipping clustering.')
            results = [range(ncells)]
        else:
            ## NOTE: this is a temporary short circuit of clustering when there are more than
            ## 500 cells. We will revisit this module and fix the issue later.
            if ncells < 500:
                results = cluster_jedna.cluster(raw_profiles, mask, n_merge=25, score_cutoff=10)
            else:
                martian.alarm("Too many cells for clustering. Putting all cells in one cluster.")
                results = [range(ncells)]
            # if ncells else
        # if skip_clustering else
    except:
        martian.alarm("Clustering encountered an exception. Putting all cells in one cluster.")
        results = [range(ncells)]
    # try/except

    out_file = open(outs.clusters, 'w')
    out_file.write(tenkit.safe_json.safe_jsonify(results))
    out_file.close()
# main

#...............................................................................
def join(args, outs, chunk_defs, chunk_outs):
    raise Exception("Join is unimplemented")
# join

