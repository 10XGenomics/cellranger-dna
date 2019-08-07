#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import numpy as np
import scipy.cluster
import scipy.stats
import time
from collections import Counter


def compute_linkage( Q ):
    distances = scipy.cluster.hierarchy.distance.pdist(Q, metric="cityblock")

    if Q.shape[0] < 2:
        return distances, np.zeros((0, 4), dtype=float)

    ## Greedy hierarchical clustering from scipy.cluster
    ## From the scipy docs:
    ## method="complete" assigns
    ##           d(u,v)=max(dist(u[i],v[j]))
    ## for all points i
    ## in cluster u and j in cluster v. This is also known as the 
    ## Farthest Point Algorithm or Voor Hees Algorithm.
    Z = scipy.cluster.hierarchy.linkage( distances, method="complete" )

    return distances, Z

def get_descendant_matrix( Z ):
    """ Given a hierarchical clustering Z determine the leaves that are 
        descendants of a given cluster node for all nodes in the tree.
        Returns a bool matrix of shape (2*n-1) x nleaves.
    """
    nleaves = Z.shape[0] + 1
    nnodes = 2*nleaves - 1
    desc = np.zeros((nnodes, nleaves), dtype=bool)
    desc[0:nleaves, :] = np.eye(nleaves)
    for i, z in enumerate(Z):
        desc[nleaves+i] = np.logical_or(desc[int(z[0])], desc[int(z[1])])
    return desc

def get_cluster_node_profile( Y, node_id, desc ):
    """ Given a node_id return the aggregated profile corresponding to the 
        descendant leaves.
    """
    return Y[desc[node_id], :].sum(axis=0)

def compute_heterogeneity(Q, Z, cmask, windows, levels=6):
    """ For each internal node of the tree compute the heterogeneity of
        CNV calls across all cells that are associated with the node.
        Returns a (ncells - 1, nbins) matrix. At a given genome bin, we 
        treat each ploidy level as an "allele", and heterogeneity is defined 
        as 1.0 - max allele fraction. Heterogeneity is assigned special values
        - in a masked region, it is -1
        - we ignore heterogeneity in a neighborhood (blur window) around a 
          breakpoint. Rather than zeroing them out, we simply add a - sign. The
          size of the blur window is the max of window sizes for each cell in the
          cluster.
    """
    t0 = time.time()

    ninodes = Z.shape[0]
    ncells = ninodes + 1
    nbins = Q.shape[1]

    desc = get_descendant_matrix(Z)
    nchunks = 100
    chunk_size = max(int(np.ceil(nbins*1.0/nchunks)), 1)
    start = 0
    het = -np.ones((ninodes, nbins), dtype="float32")
    node_size = ((desc.sum(axis=1)[ncells:])[:,
        np.newaxis, np.newaxis]).astype(float)

    ## break up the genome into chunks and compute heterogeneity over
    ## each chunk. This is done because we create a data structure of
    ## size (# nodes, # levels, # bins) which would be huge if we did not
    ## chunk.
    while start < nbins:
        end = min(start+chunk_size, nbins)
        mbins = end - start
        ## the ploidy distribution for a given internal node of the tree
        occupancy = -np.ones((ninodes, levels, mbins), dtype="float32")
        for inode in xrange(ninodes):
            lnode, rnode = map(int, Z[inode, :2])
            ## The ploidy distribution over the left node
            if lnode < ncells:
                M = np.clip(Q[lnode, start:end], None, levels-1)
                ## Note: M has the -128 calls aka no calls intact
                lo = np.zeros((levels, mbins))
                for level in xrange(levels):
                    lo[level] = (M == level)
            else:
                lo = occupancy[lnode - ncells, :, :]
            
            ## The ploidy distribution over the right node
            if rnode < ncells:
                M = np.clip(Q[rnode, start:end], None, levels-1)
                ro = np.zeros((levels, mbins))
                for level in xrange(levels):
                    ro[level] = (M == level)
            else:
                ro = occupancy[rnode - ncells, :, :]
            occupancy[inode, :, :] = lo + ro
        occupancy /= node_size
        
        het[:, start:end] = 1 - np.max(occupancy, axis=1)
        start = end
    
    assert (het < 0).sum() == 0, "het uninitialized"
    
    ## map from masked genome to full genome
    cmasksum = np.cumsum(cmask)
    nmaskedbins = cmasksum[-1]
    cmap = -np.ones(cmask.sum(), dtype=int)
    for i, v in enumerate(cmask):
        cmap[cmasksum[i]-1] = i
    
    ## we expect heterogeneity simply due to breakpoint positional error. 
    ## Make the het values negative here to indicate this.
    for inode in xrange(ninodes):
        cell_filter  = desc[inode+ncells]
        blur_window = int(np.max(windows[cell_filter])/2)
        mploidy = Q[ncells+inode, cmask]
        mbpos = np.where(np.diff(mploidy) != 0)[0]
        zeromask = np.zeros_like(cmask)
        
        for mb in mbpos:
            mbmin = max(mb - blur_window, 0)
            mbmax = min(mb+blur_window, nmaskedbins)
            bmin = cmap[mbmin]
            bmax = cmap[mbmax] if mbmax < len(cmap) else cmap[-1]
            zeromask[bmin:bmax] = True
        #zero_frac = zeromask.sum()*1.0/zeromask.shape[0]
        het[inode, zeromask] *= -1
    
    ## het values over masked genome are -1
    het[:, ~cmask] = -1

    print "Time elapsed %.2f s"%(time.time()- t0)
    return het

def tree_cut_het2( root, het_agg, hblocks, desc, clusters,
    pvalue_func1, pvalue_threshold, verbose=False ):
    """ Traverse the tree from the root towards the leaves. At each node in the
        traversal, decide whether to declare the node as a cluster by itself
        or split into the left and right child nodes. The decision is made by
        computing the maximum block of heterogeneity at the node. If this block
        is larger than the smallest observable value, then split. Otherwise,
        declare the node a cluster.
    """
    queue = [root]
    while len(queue) > 0:
        new_root = queue.pop()
        left = new_root.get_left()
        right = new_root.get_right()
        root_node = new_root.get_id()

        if left is None and right is None:
            clusters.append(root_node)
            continue
        ncells = len(hblocks)+1
        blocks = hblocks[root_node-ncells]
        max_het_block = max(blocks) if len(blocks) else 0
        pvalue = pvalue_func1(max_het_block)
        cells_per_node = desc[root_node].sum()
        if verbose:
            print "%4d: %3d cells, het block: max %d, total %d"%(root_node,
                desc[root_node].sum( ), np.max(blocks), np.sum(blocks)),
            print "pvalue", pvalue,
        pvalue_correction = cells_per_node*(cells_per_node-1)/2
        if pvalue >= pvalue_threshold/pvalue_correction:
            clusters.append(root_node)
            if verbose:
                print " <- CLUSTER"
        else:
            if verbose:
                print
            queue.append(right)
            queue.append(left)

def estimate_log_a(Z, hblocks):
    nmerges = 3
    hb_counter  =Counter()
    count_merges = 0
    for i in xrange(Z.shape[0]):
        if Z[i][3] > 2:
            continue
        count_merges += 1
        hb_counter += Counter(hblocks[i])
        if count_merges == nmerges:
            break
    if sum(hb_counter.values()) > 20 and count_merges == nmerges:
        xvs = np.arange(1, 10)
        yvs = np.array([hb_counter[x]*1.0 for x in xvs])
        yvs/= yvs.sum()
        log_y = np.log(yvs)
        log_a, _ = np.polyfit((xvs-1)[~np.isinf(log_y)],
            log_y[~np.isinf(log_y)], 1)
    else:
        log_a = -2.0
    return log_a

def tree_cut_het( root, het_agg, hblocks, desc, clusters, threshold,
    verbose=False ):
    """ Traverse the tree from the root towards the leaves. At each node in the
        traversal, decide whether to declare the node as a cluster by itself
        or split into the left and right child nodes. The decision is made by
        computing the maximum block of heterogeneity at the node. If this block
        is larger than the smallest observable value, then split. Otherwise,
        declare the node a cluster.
    """ 
    left = root.get_left()
    right = root.get_right()
    root_node = root.get_id()
    
    if left is None and right is None:
        clusters.append(root_node)
        return
 
    ncells = len(hblocks)+1
    blocks = hblocks[root_node-ncells]
    if verbose:
        print "%4d: %3d cells, het block: max %d, total %d"%(root_node, 
            desc[root_node].sum( ), np.max(blocks), np.sum(blocks)),

    if np.max(blocks) < threshold:
        clusters.append(root_node)
        if verbose:
            print " <- CLUSTER"
    else:
        if verbose:
            print
        tree_cut_het(left, het_agg, hblocks, desc, clusters, threshold,
            verbose=verbose)
        tree_cut_het(right, het_agg, hblocks, desc, clusters, threshold,
            verbose=verbose)

def compute_het_blocks(het, min_frac = 0.01):
    """ Compute blocks of heterogeneity (defined by min_frac) for all nodes.
    """
    hblocks = []
    for h in het:
        blocks = []
        i = 0
        while i < len(h):
            j = i+1
            while j < len(h):
                if (h[i] >= min_frac) != (h[j] >= min_frac):
                    break
                j += 1
            if h[i] > min_frac:
                blocks.append(j-i)
            i = j
        hblocks.append(blocks)
    return hblocks

def tree_cut( root, Z, D, desc, clusters,
    pvalue_threshold=1e-10, max_dist = 1000,
    verbose = True ):
    """ Traverse the tree from the root towards the leaves. At each node in the
        traversal, decide whether to declare the node as a cluster by itself
        or split into the left and right child nodes. The decision is made by
        computing
        -   the "inter" distances between cells in the left child and cells in
            the right child
        -   the "intra" distances between cells in the left child
        -   the "intra" distances between cells in the right child
        Use a Welch t test to decide whether the "inter" distances and the
        "intra" distances belong to the same distribution or not. If they are
        the same distribution, then keep the node intact, otherwise split.

        pvalue_threshold: determines the pvalue at which we split. Increasing
        this parameter makes us split more permissively resulting in smaller
        clusters.
        max_dist: when we cannot compute the t test because of insufficient data
        just use a distance threshold of max_dist to decide whether to split
    """
    def intra_dist( leaves, D, eps = 1e-3 ):
        intra = []
        for i in xrange(len(leaves)):
            for j in xrange(i+1, len(leaves)):
                intra.append( D[leaves[i],leaves[j]] + eps )
        return intra
    def median_intra( leaves, D ):
        intra = intra_dist( leaves, D )
        d_intra = np.median(intra) if len(intra) else 0
        return d_intra

    root_node = root.get_id( )
    left = root.get_left( )
    right = root.get_right( )
    if left is None and right is None:
        clusters.append(root_node)
        return

    left_leaves = np.where(desc[left.get_id( )])[0]
    right_leaves = np.where(desc[right.get_id( )])[0]

    inter_dist = []
    for l in left_leaves:
        for r in right_leaves:
            inter_dist.append( D[l, r] )

    d_inter = np.median(inter_dist) if len(inter_dist) else 0

    left_intra  = intra_dist( left_leaves, D )
    right_intra = intra_dist( right_leaves, D )

    ttest_left = scipy.stats.ttest_ind( left_intra, inter_dist,
        equal_var=False )[1]
    ttest_right = scipy.stats.ttest_ind( right_intra, inter_dist,
        equal_var=False )[1]

    if verbose:
        print "Root:", root_node
        print "# left:", len(left_leaves), "# right:", len(right_leaves)
        print "Left cells: ", left_leaves[:10], "..."
        print "Right cells: ", right_leaves[:10], "..."
        print len(left_intra), len(right_intra), len(inter_dist)
        if len(inter_dist):
            print "left-inter", ttest_left,
            print "right-inter", ttest_right

    if len(left_leaves) <= 2 and len(right_leaves) <= 2:
        split = d_inter > max_dist
    elif len(left_leaves) <= 2:
        assert not np.isnan(ttest_right)
        split  = ttest_right < pvalue_threshold
    elif len(right_leaves) <= 2:
        assert not np.isnan(ttest_left)
        split = ttest_left < pvalue_threshold
    else:
        assert not np.isnan(ttest_right)
        assert not np.isnan(ttest_left)
        split = min(ttest_left, ttest_right) < pvalue_threshold

    if verbose:
        if split:
            print "SPLIT"
        print "-"*40

    if not split:
        ## keep cluster intact
        clusters.append(root_node)
    else:
        ## split and recursively run on child nodes
        tree_cut( left, Z, D, desc, clusters, pvalue_threshold, max_dist,
            verbose  )
        tree_cut( right, Z, D, desc, clusters, pvalue_threshold, max_dist,
            verbose )


def traverse( root, Y, desc, clusters, min_corr ):
    """ Recursively traverse tree from root towards leaves. Y = profiles per cell.
        For each node, look at how the left and right profiles compare using
        Pearson correlation. If correlation < min_corr, keep cluster intact,
        otherwise split cluster.
    """
    root_node = root.get_id( )
    left = root.get_left( )
    right = root.get_right( )
    if left is None and right is None:
        clusters.append(root_node)
        return
    ## get left and right nodes
    lnode = left.get_id( )
    rnode = right.get_id( )

    ## construct read count profiles for each
    P_left = get_cluster_node_profile( Y, lnode, desc )
    P_right = get_cluster_node_profile( Y, rnode, desc )

    ## compute correlation
    corr = get_corr( P_left, P_right )

    if corr > min_corr:
        ## keep cluster intact
        clusters.append(root_node)
    else:
        ## split and recursively run on child nodes
        traverse( left, Y, desc, clusters, min_corr=min_corr )
        traverse( right, Y, desc, clusters, min_corr=min_corr )

def cut_max_diameter( root, Z, max_dist, clusters ):
    """ Recursively traverse tree from root towards leaves. If the child nodes
        of the root are more than distance max_dist apart, then split. Otherwise
        keep the root node intact. Returns node IDs of flat clusters.
    """
    root_node = root.get_id( )
    left = root.get_left( )
    right = root.get_right( )
    if left is None and right is None:
        clusters.append(root_node)
        return

    nleaves = Z.shape[0]+1
    assert root_node >= nleaves, "%d root node, nleaves %d"%(root_node, nleaves)

    child_dist = Z[root_node - nleaves][2]

    if child_dist < max_dist:
        ## keep cluster intact
        clusters.append(root_node)
    else:
        ## split and recursively run on child nodes
        cut_max_diameter( left, Z, max_dist, clusters )
        cut_max_diameter( right, Z, max_dist, clusters )

def get_cluster_node_size( node_id, desc ):
    """ Given a descendant matrix and node_id, compute the number of leaves
        that are descendants of node_id."""
    return desc[node_id].sum( )

def get_corr( x, y ):
    """Pearson correlation"""
    return (np.sum(x.astype(float)*y.astype(float)))/np.sqrt(np.power(x.astype(float), 2).sum( ) * np.power(y.astype(float), 2).sum( ))

