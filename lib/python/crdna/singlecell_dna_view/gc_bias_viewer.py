#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#   
# Compute GC effect on read counts
#

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from crdna.constants import MIN_GC, MAX_GC, GC_RES

#...............................................................................
def plot_gc_bias(gc, y, xvals, yvals, a, b, c, gc_plot_file):
    ## make a nice plot
    _, ax = plt.subplots(figsize=(4,4))
    plt.scatter(gc, y, marker=".", alpha = 0.05, rasterized=True, label="raw" )
    plt.xlim(0.2, 0.7)
    plt.ylim(0, np.percentile(y, 99)*1.5)
    plt.plot(xvals, yvals, "o", color="r", label="median" )

    xvs2 = np.linspace(MIN_GC, MAX_GC, 100)
    plt.plot(xvs2, a*np.power(xvs2,2) + b*xvs2 + c, color="orange", lw=2, label="fit")
    ax.grid(True)
    plt.xlabel("GC content per %d kb bin"%(GC_RES*20))
    plt.ylabel("read counts per %d kb bin"%(GC_RES*20))
    plt.title("gc_cells_only = %.3f"%gc_metric)
    plt.savefig(gc_plot_file, bbox_inches="tight")
# plot_gc_bias
