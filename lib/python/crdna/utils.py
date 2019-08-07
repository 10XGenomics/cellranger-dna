#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import pandas as pd
import numpy as np
from plot_utils import aggregate_counts
from crdna.constants import MAPPABILITY_THRESHOLD
import os

def get_version():
    fn = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        "/../../../.version")
    if os.path.exists(fn):
        return open(fn, "r").read().strip()
    return "UnknownVersion"
#
#...............................................................................
def remove_missing_chroms(chroms, profiles):
    # Clean up input chroms:
    keys = profiles.keys()
    #print('keys:')
    #print(keys)
    remove = []
    for chrom in chroms:
        if not (('/contigs/' + chrom) in keys):
            print('Removing %s' % chrom)
            remove.append(chrom)
        # if not in
    # for chrom
    for chrom in remove:
        chroms.remove(chrom)
    # for chrom
    return(chroms)
# remove_missing_chroms
#
#...............................................................................
def load_genome_data( profiles_name, tracks_name, chroms, integer = True,
    start_cell = None, end_cell = None,
    mappability_threshold = MAPPABILITY_THRESHOLD, rounding = False ):
    #print('Entering utils.load_genome_data()')
    """
    Given h5 profiles and tracks, create a genome-wide matrix of read counts,
    a compatible mappability mask, and a record of the chromosome boundaries.
    start_cell and end_cell optionally allow only loading a subset of cells.
    """
    profiles = pd.HDFStore( profiles_name, "r" )

    chroms = remove_missing_chroms(chroms, profiles)
    
    ## load mask and define chromosome boundaries
    #bdy = [0]
    #gmask = []
    #for chrom in chroms:
    #    x = profiles["/masks/" + chrom].values
    #    gmask.extend(x)
    #    z = bdy[-1] + len(x)
    #    bdy.append(z)
    ## for chrom
    #bdy = np.array(bdy)
    #gmask = np.array(gmask)

    ## restrict to highly mappable regions if tracks_name is not None
    ## otherwise all bins are mappable
    gmask = []
    bdy = [0]
    if tracks_name is None:
        for chrom in chroms:
            cbins = profiles["/contigs/"+chrom].shape[1]
            gmask.extend( np.ones(cbins, dtype=bool) )
            z = bdy[-1] + cbins
            bdy.append(z)
        # for chrom
    else:
        maptrack = pd.HDFStore(tracks_name, "r")
        for chrom in chroms:
            x = maptrack["/map/"+chrom].values 
            gmask.extend( x > mappability_threshold )
            z = bdy[-1] + len(x)
            bdy.append(z)
        # for chrom
        maptrack.close( )
    # if tracks else
    gmask = np.array(gmask)
    bdy = np.array(bdy)

    nbins  = bdy[-1]
    
    assert (start_cell is None) == (end_cell is None)

    ## find the number of cells
    total_cells = 0
    for chrom in chroms:
        total_cells, _ = profiles["/contigs/"+chrom].shape
        break
    # for chrom

    start = 0 if start_cell is None else start_cell
    end = end_cell
    ncells = total_cells if start_cell is None else end_cell - start_cell
    
    #print('ncells=%d, total_cells=%d, start_cell=%d, end_cell=%d' %
    #      (ncells, total_cells, start_cell, end_cell))
    assert ncells >= 0 and ncells <= total_cells
    dtype = "int32" if integer else "float32"
    X = np.zeros((ncells, nbins), dtype=dtype)
    for ci, chrom in enumerate(chroms):
        P = profiles["/contigs/"+chrom].values[start:end, :]
        if rounding:
            X[:, bdy[ci]:bdy[ci+1]] = np.round(P).astype(dtype)
        else:
            X[:, bdy[ci]:bdy[ci+1]] = P.astype(dtype)
    profiles.close( )

    return X, gmask, bdy
# load_genome_data

def load_track_data(tracks_name, chroms, track):
    store = pd.HDFStore(tracks_name, "r")
    keys = set(store.keys())
    stack = []
    for chrom in chroms:
        key = "/{}/{}".format(track, chrom)
        if key in keys:
            stack.append(store[key].values)
    store.close()
    return np.hstack(stack)


def load_gc_data(tracks_name, chroms):
    return load_track_data(tracks_name, chroms, "GC")


def load_mappability_data(tracks_name, chroms):
    return load_track_data(tracks_name, chroms, "map")
#
#...............................................................................
# TODO: Move this to view (crdna/singlecell_dna_view/)
def add_cell_plots_to_axes(Y, axes, bdy=None, plot_res=10, plot_args = {}):
    ncells, nbins = Y.shape
    for i in xrange(ncells):
        z = aggregate_counts(Y[i], plot_res)
        axes[i].plot(z/z.mean(), rasterized=True, **plot_args)
        axes[i].set_ylim(0, 5)
        axes[i].set_yticks(range(1, 5))
        if bdy is not None:
            for x in bdy:
                axes[i].axvline(x/plot_res, color="k")
            if i < ncells - 1:
                axes[i].set_xticks([])
            else:
                axes[i].set_xticks(bdy/plot_res)
                axes[i].set_xticklabels(axes[i].get_xticks( ).astype(int), rotation=45)
            # if i else
        # if bdy
        axes[i].set_xlim(0, nbins/plot_res)
    # for i
# add_cell_plots_to_axes

def create_chunks( input_chunks, max_chunks ):
    """
        Divide up the chunks into at most max_chunks chunks. Input is a list 
        and output is a list of lists.
    """
    assert max_chunks > 0
    ninput = len(input_chunks)
    items_per_chunk = ninput / max_chunks + int(ninput % max_chunks != 0)
    start = 0
    chunks = []
    while start < ninput:
        end = min(start+items_per_chunk, ninput)
        chunks.append(input_chunks[start:end])
        start += items_per_chunk
    return chunks

def load_h5( name, key ):
    """ Return pandas object associated to HDFStore name. """
    store = pd.HDFStore(name, "r")
    X = store[key]
    store.close( )
    return X

def load_mask_bdy(profile_name, chroms):
    """ Load genome mask and chromosome boundaries. """
    store = pd.HDFStore(profile_name, "r")
    cmask = []
    bdy = [0]
    for chrom in chroms:
        cmask.extend(store["/masks/"+chrom].values)
        bdy.append(len(cmask))
    cmask = np.array(cmask)
    bdy = np.array(bdy)
    store.close( )
    return cmask, bdy
