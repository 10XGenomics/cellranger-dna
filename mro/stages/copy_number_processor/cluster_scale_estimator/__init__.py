#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import json
from longranger.cnv import coverage_matrix
import pandas as pd
import numpy as np
import crdna.singlecell_dna_cnv.scale_estimator as sp
from crdna.utils import create_chunks

__MRO__ = """
# Inputs: cluster indices from previous stage, normalized single cell profiles at 20kb
# Output: scaled 20 kb clustered profiles, normalized
stage ESTIMATE_CLUSTER_SCALE(
    in  string  reference_path,
    in  h5      normalized_singlecell_profiles,
    in  json    clusters,
    #   
    out float[] cluster_scaling_factors,
    out h5      scaled_cluster_profiles,
    # 
    src py      "stages/copy_number_processor/cluster_scale_estimator",
) split using (
    in  int[][] clusters,
)       
"""

################################################################################
## Split on clusters
def split(args):
    if args.clusters is None:
        ncells = coverage_matrix.get_num_cells( args.normalized_singlecell_profiles,
                                                args.reference_path )
        clusters = [[x] for x in xrange(ncells)]
    else:
        f = open(args.clusters)
        clusters = json.load(f)
    # create chunks imposing MAX_CHUNKS
    MAX_CHUNKS = 30
    cluster_chunks = create_chunks(clusters, MAX_CHUNKS)
    chunks = []
    for chunk in cluster_chunks:
        chunks.append({"clusters": chunk})
    #
    return {'chunks': chunks}
# split

################################################################################
## Load data and masks
def load_data(file_name, reference_path):
    profiles, mask = coverage_matrix.load_matrix(file_name, reference_path)
    return((profiles, mask))
# load_data

################################################################################
## Compute scale from the profile for a specific choice of cluster
## and contig.
def main(args, outs):
    #
    profiles, mask = load_data(args.normalized_singlecell_profiles, args.reference_path)
    #
    # Find cluster scaling factors for each cluster
    cluster_scaling_factors = []
    for cluster in args.clusters:
        print(cluster)
        cluster_scaling_factors.append(sp.get_scale(profiles, mask, cluster))
    outs.cluster_scaling_factors = cluster_scaling_factors
    #
    #
    # Write out scaled profiles for each chunk
    #
    #
    store = pd.HDFStore( outs.scaled_cluster_profiles, "w" )
    chroms = coverage_matrix.list_primary_contigs(args.normalized_singlecell_profiles,
        args.reference_path)
    for chrom_index, chrom in enumerate(chroms):
        nclusters = len(args.clusters)
        chrom_bins = profiles[chrom_index].shape[1]
        # construct the scaled matrix for each chrom
        S = np.zeros((nclusters, chrom_bins), dtype="float32")
        for i, c in enumerate(args.clusters):
            S[i] = (profiles[chrom_index][c, :].sum(axis=0).astype("float32") / 
                    cluster_scaling_factors[i])
        # and write to disk
        store["/contigs/"+chrom] = pd.DataFrame(S)
    store.close( )

# main

################################################################################
## Join step.  Each chunk returns a floating point number per cluster. 
## This just concatenates the scale factors into an array
## 
def join(args, outs, chunk_defs, chunk_outs):
    #
    clusters = []
    for chunk in chunk_defs:
        clusters.extend(chunk.clusters)
    #
    cluster_scaling_factors = []
    for chunk in chunk_outs:
        cluster_scaling_factors.extend(chunk.cluster_scaling_factors)
    outs.cluster_scaling_factors = cluster_scaling_factors
    #
    nclusters = len(clusters)
    chroms = coverage_matrix.list_primary_contigs(args.normalized_singlecell_profiles,
        args.reference_path)
    out_store = pd.HDFStore(outs.scaled_cluster_profiles, "w")
    for chrom_index, chrom in enumerate(chroms):
        # this is the final nclusters x chrom_bins matrix for the chrom
        X = None
        start = 0
        for chunk_out in chunk_outs:
            in_store = pd.HDFStore(chunk_out.scaled_cluster_profiles, "r")
            chunk = in_store["/contigs/"+chrom].values
            in_store.close( )
            clusters_per_chunk, chrom_bins = chunk.shape
            if X is None:
                # initialize X
                X = np.zeros((nclusters, chrom_bins), dtype="float32")
            # add piece to X
            X[start:start+clusters_per_chunk, :] = chunk
            start += clusters_per_chunk
        assert start == nclusters, "incorrectly implemented split"
        out_store["/contigs/"+chrom] = pd.DataFrame(X)
        # add masks and constants from original profiles
        in_store = pd.HDFStore(args.normalized_singlecell_profiles, "r")
        out_store["/masks/"+chrom] = in_store["/masks/"+chrom]
        out_store["constants"] = in_store["constants"]
        in_store.close( )

    out_store["/clusters"] = pd.Series(range(nclusters))
    out_store.close( )
# join

################################################################################
##

################################################################################
def export_scaled_cluster_profiles(profiles, mask, chroms, clusters, cluster_scaling_factors, out_file, window_size):
    S = []
    n_clusters = len(clusters)
    n_chrom = len(profiles)
    for chrom in range(n_chrom):
        N = profiles[chrom]
        n_bins = N.shape[1]
        chr_block = pd.DataFrame(index=range(n_clusters), columns=range(n_bins))
        for cluster_index in range(n_clusters):
            cluster = clusters[cluster_index]
            # Add up cell profiles belonging to this cluster:
            tmp = N[cluster].sum(axis=0)
            # Divide the aggregated cluster profile by the appropriate scale factor:
            tmp /= cluster_scaling_factors[cluster_index]
            # Store the result in chr_block:
            print('=' * 80)
            print('chrom=%d' % chrom)
            print('cluster_index=%d' % cluster_index)
            print('tmp.shape:')
            print(tmp.shape)
            print('chr_block.shape:')
            print(chr_block.shape)
            chr_block.loc[cluster_index] = tmp
        # for cluster_index
        S.append(chr_block)
    # for chrom
    #   
    ## Store scaled profiles
    store = pd.HDFStore(out_file, "w" )
    for i, chrom in enumerate(chroms):
        store["/contigs/"+chrom] = pd.DataFrame(S[i])
        store["/masks/"+chrom] = pd.Series(np.ones(S[i].shape[0], dtype=bool))
    store["/clusters"] = pd.Series(range(n_clusters))
    store["/constants"] = pd.Series({"window_size": window_size})
    store.close( )
# export_scaled_cluster_profiles


