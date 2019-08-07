#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Compute GC effect on read counts
#

import pandas as pd
import numpy as np
import json
import longranger.cnv.contig_manager as contig_manager
from plot_utils import aggregate_counts
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from crdna.constants import MAPPABILITY_THRESHOLD, GC_RES, MIN_GC, MAX_GC, NUM_GC_BINS, MIN_POINTS_PER_BIN

#...............................................................................
def estimate_gc_bias(profiles, tracks, reference_path):
    ## load genome tracks and profiles skipping sex chromosomes
    
    ref = contig_manager.contig_manager(reference_path)
    chroms = ref.primary_contigs(allow_sex_chromosomes=False)

    maptrack = pd.HDFStore(tracks, "r")
    cmask = []
    gctrack = []
    bdy = [0]
    mtrack = []
    for chrom in chroms:
        x = maptrack["/map/"+chrom].values > MAPPABILITY_THRESHOLD
        cmask.extend(x)
        z = bdy[-1] + len(x)
        gctrack.extend(maptrack["/GC/"+chrom].values)
        mtrack.extend(maptrack["/map/"+chrom].values)
        bdy.append(z)
    cmask = np.array(cmask)
    maptrack.close( )
    gctrack = np.array(gctrack)
    mtrack = np.array(mtrack)
    bdy = np.array(bdy)

    nbins = bdy[-1]
    pstore = pd.HDFStore(profiles, "r")
    ncells = len(pstore["/barcodes"].values)
    X = np.zeros((ncells, nbins), dtype="int32")
    for ci, chrom in enumerate(chroms):
        X[:, bdy[ci]:bdy[ci+1]] = pstore["/contigs/"+chrom].values
    pstore.close( )

    ## genome wide profile of all cells @ GC_RES resolution
    ## restricted to mappable regions
    y = aggregate_counts(X.sum( axis=0 )[cmask], GC_RES).astype(float)
    y /= y.mean( )
    gc = aggregate_counts(gctrack[cmask], GC_RES)/GC_RES

    gcbins = np.linspace(MIN_GC, MAX_GC, NUM_GC_BINS+1)
    gc_vals = 0.5 * (gcbins[1:] + gcbins[:-1])
    gc_bin_index = np.searchsorted(gcbins, gc)
    gc0 = np.nanmean(gc_vals)

    ## group data points by GC bins and compute the median
    x_vals = []
    y_vals = []
    for bi in xrange(1, NUM_GC_BINS+1):
        bin_filter = gc_bin_index == bi
        num_data_points = bin_filter.sum( )
        if num_data_points < MIN_POINTS_PER_BIN:
            continue
        bin_gc = gc_vals[bi-1]
        bin_val = np.median(y[bin_filter])
        x_vals.append(bin_gc)
        y_vals.append(bin_val)
    # for bi
    x_vals = np.array(x_vals) - gc0
    
    ## fit to ax^2 + bx + c
    a, b, c = np.polyfit(x_vals, y_vals, 2)
    
    ## GC metric is mean absolute deviation away from 1.0
    gc_metric = np.abs(np.array(y_vals) - 1.0).sum( ) / len(y_vals)

    ## store gc data in summary
    summary = {}
    summary["GC_content"] = x_vals
    summary["scaled_read_counts"] = y_vals
    summary["quadratic_coefficients"] = [a, b, c]
    summary["gc_cells_only"] = gc_metric
    summary["gc0"] = gc0
   
    #with open(outs.summary, "w") as out:
    #    json.dump(summary, out, indent=4)
    #
    return( {'GCMetric': gc_metric, 'Summary': summary})
# estimate_gc_bias
