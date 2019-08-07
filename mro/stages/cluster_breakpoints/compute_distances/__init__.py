#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Compute all pair-wise distances between a set of points
#

import pandas as pd
import numpy as np
from plot_utils import aggregate_matrix

__MRO__ = """
stage COMPUTE_DISTANCES(
    in  h5       features     "num_cells x n matrix",
    in  string   feature_key  "which h5 key to use",
    out h5       distances    "upper triangle of distance matrix (flattened)",
    src py       "stages/cluster_breakpoints/compute_distances",
) split using (
    in  map      pairs,
)
"""

def split(args):

    store = pd.HDFStore( args.features, "r" )
    ncells, nfeatures = store[args.feature_key].shape
    store.close( )

    npairs = ncells * (ncells - 1)/2
    MAX_CHUNKS = 30
    MIN_PAIRS_PER_CHUNK = 300000
    
    nchunks = npairs / MIN_PAIRS_PER_CHUNK + int(npairs % MIN_PAIRS_PER_CHUNK != 0)
    nchunks = np.clip(nchunks, 1, MAX_CHUNKS)
    pairs_per_chunk = npairs/nchunks + int(npairs % nchunks != 0)

    pi = 0
    chunk_defs = []
    for ci in xrange(nchunks):
        chunk = {"start" : pi, "end" : min(npairs, pi + pairs_per_chunk), "n": ncells}
        pi += pairs_per_chunk
        chunk_defs.append({"pairs"    : chunk,
                           "__mem_gb" : 24})

    return {'chunks': chunk_defs, 'join': {'__mem_gb': 12}}

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    start_pair = args.pairs["start"]
    end_pair   = args.pairs["end"]
    ncells     = args.pairs["n"]
    npairs     = end_pair - start_pair

    start_i, start_j = translate_pair_index(start_pair, ncells)
    
    ## load the window size for segmentation
    store = pd.HDFStore( args.features, "r" )
    segment_windows = store["constants"]["segment_windows"]
    store.close( )

    ## aggregate features at segment window size to improve
    ## break point position accuracy
    X = aggregate_matrix(
        load_features(args.features, args.feature_key).astype("int32"),
        segment_windows)/segment_windows
    
    i = start_i
    j = start_j
    distances = np.zeros(npairs, dtype="float32")

    di = 0
    while i < ncells:
        while j < ncells:
            distances[di] = dist_L1(X[i], X[j])
            j  += 1
            di += 1
            if di == npairs:
                break
        if di == npairs:
            break
        i += 1
        j =  i+1
    
    store = pd.HDFStore( outs.distances, "w" )
    store["/distances"] = pd.Series(distances)
    store.close( )

def join(args, outs, chunk_defs, chunk_outs):
    
    distances = []
    for chunk_out in chunk_outs:
        instore = pd.HDFStore( chunk_out.distances, "r" )
        distances.extend( instore["/distances"].values )
        instore.close( )

    store = pd.HDFStore( outs.distances, "w" )
    store["/distances"] = pd.Series(distances)
    store.close( )

## support functions

def dist_L1( x, y ):
    return np.abs(x - y).sum( )

def load_features( features_file, feature_key ):
    store = pd.HDFStore( features_file )
    features = store[feature_key].values
    store.close( )
    return features

def translate_pair_index(q, n):
    """translate the index of the flattened upper triangle into an (i, j)"""
    bins = np.arange( 1, n )[::-1]
    cbins = np.cumsum(bins)
    i = np.searchsorted(cbins, q, "right")
    if i == 0:
        j = q + i + 1
    else:
        j = q - cbins[i-1] + i + 1
    return (i, j)
