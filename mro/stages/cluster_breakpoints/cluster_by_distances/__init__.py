#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import pandas as pd
import scipy.cluster

__MRO__ = '''
stage CLUSTER_BY_DISTANCES(
    in  h5       distances,
    in  h5       features,
    in  string   feature_key,
    out h5       data,
    src py       "stages/cluster_breakpoints/cluster_by_distances",
)
'''

def main( args, outs ):
    args.coerce_strings()
    outs.coerce_strings()
    
    ## These are L1-norm/manhattan distances
    store = pd.HDFStore( args.distances, "r" )
    distances = store["/distances"].values
    store.close( )

    ## These are L1-norm/manhattan distances
    store = pd.HDFStore( args.distances, "r" )
    distances = store["/distances"].values
    store.close( )
    
    #store = pd.HDFStore( args.features, "r" )
    #features = store[args.feature_key].values
    #store.close( )
    
    ## Greedy hierarchical clustering from scipy.cluster
    ## From the scipy docs:
    ## method="complete" assigns
    ##           d(u,v)=max(dist(u[i],v[j]))
    ## for all points i
    ## in cluster u and j in cluster v. This is also known as the 
    ## Farthest Point Algorithm or Voor Hees Algorithm.
    Z = scipy.cluster.hierarchy.linkage( distances, method="complete" )

    out_store = pd.HDFStore( outs.data, "w")
    out_store["/Z"] = pd.DataFrame(Z)
    out_store.close( )

