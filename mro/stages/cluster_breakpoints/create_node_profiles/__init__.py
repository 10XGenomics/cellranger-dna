#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import pandas as pd
import numpy as np

from longranger.cnv import contig_manager
from crdna.utils import load_h5

__MRO__ = '''
stage CREATE_NODE_PROFILES(
    in  string reference_path,
    in  h5     profiles,
    in  h5     cluster_data,
    out h5     node_profiles,
    out h5     internal_node_profiles,
    src py     "stages/cluster_breakpoints/create_node_profiles",
) split using ()
'''

def split(args):
    ref = contig_manager.contig_manager( args.reference_path )
    chroms = ref.list_all_contigs( )
    max_chrom_size = max([ref.contig_lengths[chrom] for chrom in chroms])
    constants = load_h5(args.profiles, "constants").to_dict()
    ncells = constants["ncells"]
    window_size = constants["window_size"]
    max_mat_size_gb = float(2*ncells*max_chrom_size/window_size)/1e9*4
    mem_gb = int(np.ceil(max_mat_size_gb*4 + 1))
    return {'chunks' : [], 'join': {'__mem_gb': mem_gb}}

def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()

    ref = contig_manager.contig_manager( args.reference_path )
    chroms = ref.list_all_contigs( )

    # load tree
    store = pd.HDFStore( args.cluster_data, "r" )
    Z = store["/Z"].values
    cell_scale_factors = store["scale_factor"].values
    cell_cov = store["reads_per_bin"].values
    store.close( )

    print Z
    print cell_scale_factors
    print cell_cov

    ncells = Z.shape[0] + 1
    nnodes = 2*ncells - 1

    #desc = get_descendant_matrix( Z )    

    in_store = pd.HDFStore( args.profiles, "r" )
   
    node_store = pd.HDFStore( outs.node_profiles, "w" )
    internal_node_store = pd.HDFStore( outs.internal_node_profiles, "w" )
    node_store["/constants"] = in_store["/constants"]
    internal_node_store["/constants"] = in_store["/constants"]

    ## compute scale factors for internal nodes. Each node is assumed to 
    ## scale like its constituent cells. Store in internal node profiles.
    scale_guess = np.zeros(len(Z))
    scale_factors = -np.ones(nnodes)
    cov = -np.ones(nnodes)
    scale_factors[:ncells] = cell_scale_factors
    cov[:ncells] = cell_cov

    for i, z in enumerate(Z):
        node = ncells+i
        left, right = z.astype(int)[:2]
        sf_left, sf_right = scale_factors[[left, right]]
        assert sf_left >= 0 and sf_right >= 0
        lsize = int(Z[left-ncells, 3]) if left >= ncells else 1
        rsize = int(Z[right-ncells, 3]) if right >= ncells else 1
        nsize = z[-1]
        assert nsize == lsize+rsize, "tree structure invalid"
        cov_left, cov_right = cov[[left, right]]
        assert cov_left >= 0 and cov_right >= 0, "cov not initialized"
        cov_node = cov_left+cov_right
        cov[node] = cov_node
        if cov_node > 0:
            scale_guess[i] = nsize/cov_node*(sf_left*cov_left/lsize +
                sf_right*cov_right/rsize)
        else:
            scale_guess[i] = 0
        scale_factors[node] = scale_guess[i]
    internal_node_store["/scale_guess"] = pd.Series(scale_guess)

    for chrom in chroms:
        X = in_store["/contigs/"+chrom].values
        nbins = X.shape[1]

        ## first store all the cell profiles
        Y = np.zeros( (nnodes, nbins), dtype=X.dtype )
        Y[:ncells, :] = X

        ## next all the internal nodes
        for inode, row in enumerate(Z):
            left, right = row[:2].astype(int)
            Y[inode+ncells, :] = Y[left, :] + Y[right, :]
        #for i in xrange(nnodes):
        #    Y[i] = X[desc[i], :].sum(axis=0)
        node_store["/contigs/"+chrom] = pd.DataFrame(Y)
        node_store["/masks/" + chrom] = in_store["/masks/" + chrom]
        
        ## next store all the internal nodes
        internal_node_store["/contigs/"+chrom] = pd.DataFrame(Y[ncells:, :])
        internal_node_store["/masks/" + chrom] = in_store["/masks/" + chrom]

    node_store.close( )
    internal_node_store.close( )
    in_store.close( )


