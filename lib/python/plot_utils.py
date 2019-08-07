import os
import pandas as pd
import numpy as np
import scipy
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.rcParams['figure.figsize'] = (12.0, 9.0)
matplotlib.rcParams.update({'font.size': 14})
import matplotlib.gridspec as gridspec

import sys
from statsmodels.nonparametric.smoothers_lowess import lowess
import time

import longranger.cnv.coverage_matrix as coverage_matrix
#import longranger.cnv.contig_manager as contig_manager

def aggregate_counts( vec, window ):
    n = (len(vec)-1)/window+1
    agg_vec = np.zeros(n, vec.dtype)
    ai = 0
    for i in xrange(0, len(vec), window):
        high = min(i + window, len(vec))
        agg_vec[ai] = vec[i:high].sum()
        ai += 1
    return agg_vec

def aggregate_matrix( M, window ):
    n = (M.shape[1]-1)/window+1
    agg_mat = np.zeros((M.shape[0], n), M.dtype)
    ai = 0
    for i in xrange(0, M.shape[1], window):
        high = min(i + window, M.shape[1])
        agg_mat[:, ai] = M[:, i:high].sum(axis=1)
        ai += 1
    return agg_mat
    
    
## load data and masks
def load_data( fn, reference_path):

    X, Xm= coverage_matrix.load_matrix(fn, reference_path)

    store = pd.HDFStore( fn, "r" )
    barcodes = np.array(store["barcodes"])
    store.close( )
    return X, Xm, barcodes

## from a ncell x nbin matrix per chromosome
## trim away outliers using when counts are outside
## using MAD score of mad_score
def get_dpcvs( M, mad_score = 6 ):
    assert len(M.shape) == 2
    ncells = M.shape[0]
    
    ## filter out outliers using (gaussian) MAD
    median = np.median(M, axis=1)[:, np.newaxis]
    mads = 1.4826*np.median(np.abs(M - median), axis=1)
    outlier_mask = np.array([np.abs(M[i] - median[i]) < mad_score*mads[i]
                             for i in xrange(M.shape[0])])
    
    ## are there repeat outliers across > 25 % of cells?
    ## if so, mask those bins
    bad_frac = np.sum((~outlier_mask).astype(int), axis=0).astype(float)/M.shape[0]
    goodbins = (outlier_mask) & (bad_frac < 0.25)
    
    m = np.zeros(ncells)
    v = np.zeros(ncells)
    dpcvs = np.zeros(ncells)
    for i in xrange(ncells):
        m[i] = M[i][goodbins[i]].mean()
        v[i] = M[i][goodbins[i]].std()**2

    fail_mask = np.logical_or(m < 1e-6, v < m)
    if fail_mask.sum() > 0:
        dpcvs[fail_mask] = np.nan
    dpcvs[~fail_mask] = np.sqrt( (v-m)[~fail_mask] ) / m[~fail_mask]
    return dpcvs

def get_genome_profile(profile_name, reference_path, conf_filter):

    X, Xm, barcodes = load_data( profile_name, reference_path)
    chroms = coverage_matrix.list_primary_contigs(profile_name, reference_path)

    nchroms = len(X)
    print nchroms
    print len(conf_filter) if conf_filter is not None else "No conf filter"
    bdy = [(None, None) for _ in xrange(nchroms)]
    ncells = X[0].shape[0]

    ## coarse grain to 100 kb and compute dpcv
    Y = [aggregate_matrix( X[chrom], 5 ) for chrom in xrange(nchroms)]
    Ym = [aggregate_counts(Xm[chrom].astype(int), 5) == 5 for chrom in xrange(nchroms)]

    nbins = sum([ym.sum() for ym in Ym])
    G = np.zeros( (ncells, nbins) )
    conf = np.ones( nbins, dtype=bool )
    if conf_filter is not None:
        si = 0
        for c in xrange(nchroms):
            cbins = Ym[c].sum()
            G[:, si:si+cbins] = Y[c][:, Ym[c]]
            print "---"
            print c
            print Ym[c]
            print si+cbins
            print conf_filter[c]
            print len(conf_filter[c])
            print len(Ym[c])
            print len(conf)
            print si
            print cbins
            conf[si:si+cbins] = conf_filter[c][Ym[c]]
            bdy[c] = (si, si+cbins)
            si += cbins
    ## compute gc track
    
    gcbed100kb = "/mnt/opt/meowmix_git/genometracks/hg19-scdna/gc_content_100kb.bed"
    hg19gc = defaultdict(list)
    for line in open(gcbed100kb):
        f = line.strip().split("\t")
        hg19gc["hg19_"+f[0]].append(float(f[-1]))
    gctrack = []
    ci = 0
    for chrom in chroms:
        if len(hg19gc[chrom]) == 0:
            continue
        print chrom, ci, len(hg19gc[chrom]), len(Ym), Ym[ci].shape
        gctrack.extend(np.array(hg19gc[chrom])[Ym[ci]])
        ci += 1
    gctrack = np.array(gctrack)
    return G, barcodes, conf, gctrack, bdy

def plot_cell_subset( G_plot, bdf_cells, subset, vlines ):
    nplots = len(subset)
    assert nplots > 0
    assert nplots < 11
    nbins = G_plot.shape[1]
    xmax = nbins + 700
    fig, ax = plt.subplots(nplots, 1, figsize = (12,9), 
                           sharex=True, sharey=True)
    if nplots == 1:
        ax = [ax]
    for ci in xrange(nplots):
        cell = subset[ci]
        rds = G_plot[cell].mean()
        y = G_plot[cell]/rds
        #ax[ci].fill_between(np.arange(nbins), y-err, y+err, alpha=0.2)
        ax[ci].plot(y)
        yticks = [0,1,2]
        ax[ci].set_yticks(yticks)
        for yt in [0.5, 1, 1.5]:
            ax[ci].axhline(yt, 0, float(nbins)/xmax, color="k", ls="dotted" )
        #ax[ci].grid(True, axis="y")
        for b in vlines:
            ax[ci].axvline(b, color="k")

        display_text = "cell %d,amp %.2f,\n"\
                       "@100kb:inserts/bin %d\n"\
                       "poisson %.2f,dpcv %.2f"%(
                        cell, bdf_cells["amp.rate"].values[cell],
                        rds/10.0, 1/np.sqrt(rds/10.0),
                        bdf_cells["dpcv"].values[cell])
        ax[ci].text(2150, 1.5, display_text, horizontalalignment="left", verticalalignment="center",
                   fontsize=12)
        ax[ci].set_xticks([])
        if ci == nplots - 1:
            ax[ci].set_xlabel("autosome confident regions")
    ax[0].set_title("Insert density depicted at 1 mb, stats at 100 kb", fontsize=14)
    plt.ylim(0, 3)
    plt.xlim(0, xmax)
    plt.subplots_adjust(hspace=0.1, wspace=0)

    return fig

def plot_all_cells(G_plot, bdf_cells, bdy, plot_name, 
        order=None, plots_per_page = 10):
    assert os.path.exists(os.path.dirname(plot_name))
    ncells = G_plot.shape[0]
    if order is None:
        order = np.arange(ncells)
    pp = PdfPages(plot_name)
    for start in xrange(0, ncells, plots_per_page):
        if start % 10 == 0:
            print start, "cells done"
        subset = order[start:start+plots_per_page]
        fig = plot_cell_subset( G_plot, bdf_cells, subset, bdy )
        pp.savefig(fig)
        plt.close()
    pp.close()

def plot_cell_gc_norm( G_auto_conf, cell, bdf_cells,
	gcnormprofile, normprofile):
    _ , ax = plt.subplots(3,1, sharey=True)

    for a in ax:
        a.grid(True, axis="y")
    ax[0].plot(aggregate_counts(G_auto_conf[cell], 10)
               , alpha=0.5, label="dpcv %.2f"%bdf_cells["dpcv"].values[cell] ) 
    ax[0].legend()
    ax[0].set_xticks([])
    ax[0].set_title("genome-wide profile - unnormalized")

    ax[1].plot(aggregate_counts(gcnormprofile[cell], 10)
               , alpha=0.5, label="dpcv %.2f"%bdf_cells["gc_dpcv"].values[cell])
    ax[1].set_title("genome-wide profile - GC normalized")
    ax[1].set_xticks([])
    ax[1].legend()

    ax[2].plot( aggregate_counts(normprofile[cell], 10)
               , alpha=0.5, label="dpcv %.2f"%bdf_cells["norm_dpcv"].values[cell] )
    ax[2].set_title("normalized by mean over all cells")
    ax[2].legend()
    plt.subplots_adjust(hspace=0.35)

def compute_normalization( G_auto_conf, gctrack, conf, bdy ):
    auto_bdy = bdy[21][1]
    print G_auto_conf.shape
    print "this takes time...."
    sys.stdout.flush()
    t0 = time.time()
    gceffect = np.zeros_like(G_auto_conf)
    ncells, nbins = G_auto_conf.shape
    for cell in xrange(ncells):
        if cell % 5 == 0:
            print cell, "cells", 
            sys.stdout.flush()
        gceffect[cell, :] = lowess(G_auto_conf[cell], gctrack[0:auto_bdy][conf[0:auto_bdy]], frac=0.05, return_sorted=False)
    print
    print (time.time()-t0)/60, "mins"
    gcnorm = gceffect/gceffect.mean(axis=1)[:, np.newaxis]
    gcnan = np.isnan(gcnorm)
    gcnorm[gcnan] = 1
    gcnormprofile = G_auto_conf/gcnorm

    avgprofile = (gcnormprofile/gcnormprofile.mean(axis=1)[:, np.newaxis]).mean(axis=0)

    normprofile = gcnormprofile/avgprofile
    nannorm = np.isnan(normprofile)
    normprofile[nannorm] = gcnormprofile[nannorm]

    return gceffect, gcnormprofile, normprofile

def detect_replicating_cells( G_plot, lam_low = 0.05, lam_high = 0.95 ):
    # distribution is gaussian at x, sig1, and gaussian at 2x, sig2
    # two gaussians are weighted by lambda
    # parameterized as z = (x, sig1, sig2, lam)
    def log_dist(point, z):
        x, s1, s2, lam = z
        return np.logaddexp(np.log(lam) + scipy.stats.norm.logpdf(point, x, s1),
                     np.log(1-lam) + scipy.stats.norm.logpdf(point, 2*x, s2))
    def neg_log_like_data(z, points):
        return -np.sum(log_dist(points, z))

    ncells = G_plot.shape[0]
    lams = np.zeros(ncells)
    for cell in xrange(ncells):
        points = G_plot[cell]
        x0, sig1, sig2, lams[cell] = scipy.optimize.fmin(neg_log_like_data, (600, 100, 100, 0.5), args=(points,),
                                                  disp=False)
    return np.where((lams > lam_low) & (lams < lam_high))[0]
        
#### OLD VERSION - DO NOT USE  #####
def plot_cell_subset_full( G_auto_conf, G_plot, gctrack, gceffect, conf,
        bdy, bdf_cells, subset ):
    nplots = len(subset)
    assert nplots > 0
    assert nplots < 11
    nbins = G_plot.shape[1]

    if conf.sum( ) == nbins:
        xlabel = "autosome ALL regions"
    else:
        xlabel = "autosome CONFIDENT regions"

    ## chrom boundaries in 1 Mb scale over conf regions
    chrom_ends = [0]
    for i in xrange(1, 22):
        prev = chrom_ends[-1]
        chrom_ends.append(prev + conf[bdy[i-1][0]:bdy[i][0]].sum()/10.0)
    auto_bdy = bdy[21][1]
    chrom_mps = []
    for i in xrange(len(chrom_ends)-1):
        chrom_mps.append( 0.5*(chrom_ends[i] + chrom_ends[i+1]) )

    xsize = 16.0
    ysize = 9.0*nplots/10.0
    fig = plt.figure(figsize = (xsize, ysize))

    gs = gridspec.GridSpec(nplots, 4,width_ratios=[5,1,1,1])
    #fig, ax = plt.subplots(nplots, 3, , sharey=True)
    #if nplots == 1:
    #    ax = [ax]
    ymax = 3
    for ci in xrange(nplots):
	cell = subset[ci]
	rds = G_plot[cell].mean()
	y = G_plot[cell]/rds
	#ax[ci].fill_between(np.arange(nbins), y-err, y+err, alpha=0.2)
	
	ax1 = plt.subplot(gs[ci,0])
	ax1.plot(y, rasterized=True)
	yticks = [0,1,2]
	ax1.set_yticks(yticks)
	for yt in [0.5, 1, 1.5]:
	    ax1.axhline(yt, color="k", ls="dotted" )
	#ax[ci].grid(True, axis="y")
	for b in chrom_ends:
	    ax1.axvline(b, color="k")
	ax1.set_xlim(0, nbins)
	ax1.set_ylim(0, ymax)
	if ci == nplots-1:
            ax1.set_xticks(chrom_mps)
            ax1.set_xticklabels(map(str, range(1,23)), rotation=45, fontsize=10)
        else:
            ax1.set_xticks([])

	ax2 = plt.subplot(gs[ci,1])
	ax2.scatter(gctrack[0:auto_bdy][conf[0:auto_bdy]], G_auto_conf[cell]/G_auto_conf[cell].mean(), 
			 marker="o", alpha=0.01, rasterized=True)
	order = np.argsort(gctrack[0:auto_bdy][conf[0:auto_bdy]])
	ax2.plot(gctrack[0:auto_bdy][conf[0:auto_bdy]][order], 
		      gceffect[cell][order]/gceffect[cell].mean(), color="k")
	ax2.text(0.5,ymax-0.1,"delta gc %.2f"%bdf_cells["delta_gc"].values[cell], horizontalalignment="center", 
		 verticalalignment="top", fontsize=12)
	ax2.set_ylim(0, ymax)
	ax2.set_yticks([])
	ax2.set_xticks([])
	    
	ax3 = plt.subplot(gs[ci,2])
	ax3.hist(y, normed=True, histtype="stepfilled", alpha=0.5, bins=np.linspace(0, ymax, 100))
	for x in [0.5, 1.0, 1.5]:
	    ax3.axvline(x, ls="dashed", color="k")
	ax3.set_xticks([])
	ax3.set_yticks([])
        ax3.set_xlim(-0.1, None)
	ax3.text(ymax-0.1, ax3.get_ylim()[1]-0.2, "reads/1 mb\n%d  "%rds, horizontalalignment="right", verticalalignment="top",
		fontsize=10)
	
	display_text = "cell %d,amp %.2f,\n"\
		       "poisson %.2f,\ndpcv %.2f"%(
			cell, bdf_cells["amp.rate"].values[cell],
			1/np.sqrt(rds/10.0),
			bdf_cells["dpcv"].values[cell])
	ax4 = plt.subplot(gs[ci,3])
	ax4.text(ax4.get_xlim()[1]/2.0, 0, display_text, horizontalalignment="center", 
		 verticalalignment="bottom", fontsize=12)
	
	ax4.set_ylim(0, ymax)
	ax4.set_xticks([])
	ax4.set_yticks([])
	
	if ci == nplots-1:
	    ax1.set_xlabel(xlabel)
	    ax2.set_xlabel("GC frac/100 kb")
	    ax2.set_xticks([0.3, 0.4, 0.5, 0.6])
            ax2.set_xticklabels([".3", ".4", ".5", ".6"], rotation=45, fontsize=10)
	    ax3.set_xticks(range(3))
            ax3.set_xticklabels(range(3), rotation=45, fontsize=10)
	    
	if ci == 0:
	    ax1.set_title("read density per cell @ 1 mb")
	    ax2.set_title("GC effect")
	    ax3.set_title("Reads/1mb")
	    ax4.set_title("Stats/100kb")
	
    #ax.set_xlabel("autosome confident regions")
    #ax[0,0].set_title("Insert density depicted at 1 mb, stats at 100 kb", fontsize=14)
    plt.ylim(0, ymax)
    plt.subplots_adjust(hspace=0.1, wspace=0)

    return fig


