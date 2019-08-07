#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import pandas as pd
import numpy as np

from longranger.cnv import contig_manager
from crdna.clustering import compute_linkage
from crdna.constants import MAPPABILITY_THRESHOLD
from plot_utils import aggregate_matrix
from crdna.utils import load_h5

__MRO__ = '''
stage CLUSTER_CELLS(
    in  h5       cnv_tracks,
    in  h5       tracks,
    in  string   reference_path,
    out h5       data,
    src py       "stages/cluster_breakpoints/cluster_by_distances",
) 
'''

def split(args):
    constants = load_h5(args.cnv_tracks, "constants").to_dict()
    matsize_gb = float(constants["ncells"]*constants["genomebins"])/1e9
    return {'chunks': [], 'join': {'__mem_gb' : int(np.ceil(matsize_gb * 12 + 2))}}

def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()
    
    ref = contig_manager.contig_manager( args.reference_path )
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)

    ## Load data
    store = pd.HDFStore( args.cnv_tracks, "r" )
    Q = store["/cnv_tracks"].values
    sf = store["/scale_factor"]
    rpb = store["/reads_per_bin"]
    segment_windows = store["constants"]["segment_windows"]
    store.close( )

    if args.tracks is None:
        gmask = np.ones(Q.shape[1], dtype=bool)
    else:
        gmask = []
        maptrack = pd.HDFStore(args.tracks, "r")
        for chrom in chroms:
            x = maptrack["/map/"+chrom].values 
            ## TODO make this consistent across stages
            gmask.extend( x > MAPPABILITY_THRESHOLD )
        maptrack.close( )
        gmask = np.array(gmask)

    ## Aggregate all cells to the same resolution and compute L1 norm
    Q_agg = np.round(aggregate_matrix( Q[:, gmask].astype("float32"), 
        segment_windows)/segment_windows).astype("int32")

    distances, Z = compute_linkage( Q_agg )

    out_store = pd.HDFStore( outs.data, "w")
    out_store["/Z"] = pd.DataFrame(Z)
    out_store["distances"] = pd.Series(distances)
    out_store["constants"] = pd.Series({"segment_windows": segment_windows})
    out_store["scale_factor"] = sf
    out_store["reads_per_bin"] = rpb
    out_store.close( )
