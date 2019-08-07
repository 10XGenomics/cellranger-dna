#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
import tenkit.safe_json
from crdna.singlecell_dna_cnv import cluster_jedna
import longranger.cnv.coverage_matrix as coverage_matrix
import martian

__MRO__ = """
stage GENERATE_CLUSTERS(
    in  h5      normalized_singlecell_profiles,
    in  string  reference_path, 
    in  bool    skip_clustering,
    #
    out json    clusters,
    # 
    src py      "stages/copy_number_processor/profile_clusterer/generate_clusters",
)   
"""

#...............................................................................
def split(args):
    raise Exception("Split is unimplemented")

#...............................................................................
def main(args, outs):
    normalized_singlecell_profiles, mask = coverage_matrix.load_matrix(
        args.normalized_singlecell_profiles, args.reference_path)

    print('DEBUG generate_final_clusters/__init__.main():')
    print('normalized_singlecell_profiles[0].shape')
    print(normalized_singlecell_profiles[0].shape)

    ncells = normalized_singlecell_profiles[0].shape[0]
    results = [range(ncells)]
    try:
        if args.skip_clustering:
            print('Skipping clustering.')
        else:
            ## NOTE: this is a temporary short circuit of clustering when there are more than
            ## 500 cells. We will revisit this module and fix the issue later.
            if True: # ncells < 500:
                # results = cluster_jedna.cluster(normalized_singlecell_profiles, mask, n_merge=25, score_cutoff=10)
                results = cluster_jedna.cluster(
                    normalized_singlecell_profiles, mask, n_merge=25, score_cutoff=5)
            else:
                martian.alarm("Too many cells for clustering. Putting all cells in one cluster.")
            # if ncells else
        # if skip_clustering else
    except Exception as error:
        martian.alarm("Clustering encountered an exception. Putting all cells in one cluster. Error: %s" % repr(error))
    # try/except
    #
    out_file = open(outs.clusters, 'w')
    out_file.write(tenkit.safe_json.safe_jsonify(results))
    out_file.close()
# main

#...............................................................................
def join(args, outs, chunk_defs, chunk_outs):
    raise Exception("Join is unimplemented")
# join

