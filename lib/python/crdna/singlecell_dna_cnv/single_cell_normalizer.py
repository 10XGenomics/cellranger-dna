#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
#import json
#import pandas as pd
#from scipy.optimize import curve_fit
#from scipy.cluster import vq
from scipy.stats import norm
import crdna.constants

#import longranger.cnv.coverage_matrix as coverage_matrix

#from crdna.singlecell_dna_cnv import cluster_jedna
#from crdna.singlecell_dna_cnv import scale_estimator as sp
#from crdna.singlecell_dna_cnv import vesna
#from crdna.singlecell_dna_cnv import differential_variance_ratio as dvr
from crdna.singlecell_dna_cnv import scale_estimator as se

#............................................................................
which = lambda lst:list(np.where(lst)[0])
#
#............................................................................
#def get_default_chromosome_names():
#    chromosomes = ['hg19_chr' + str(i) for i in range(1, 23)]
#    return(chromosomes)
# get_default_chromosome_names

#............................................................................
def get_chrom_mappability(bin_parameters, chrom_name, chrom_names):
    chrom_index = get_chrom_index(chrom_name, chrom_names)
    mappability = bin_parameters[chrom_index][5, :]
    return(np.array(mappability))
# get_chrom_mappability
#
#............................................................................
def get_chrom_gc(bin_parameters, chrom_name, chrom_names):
    chrom_index = get_chrom_index(chrom_name, chrom_names)
    gc = bin_parameters[chrom_index][3, :]
    return(gc)
# get_chrom_gc
#
#............................................................................
def get_position(bin_parameters, chrom_name, chrom_names):
    chrom_index = get_chrom_index(chrom_name, chrom_names)
    position = bin_parameters[chrom_index][1, :]
    return(position)
# get_position
#
#............................................................................
def get_chrom_index(chrom_name, chrom_names):
    chrom_index = chrom_names.index(chrom_name)
    return(chrom_index)
# get_chrom_index
#
#............................................................................
def get_gc(bin_parameters, chrom_names):
    gc = []
    for chrom_name in chrom_names:
        chrom_index = get_chrom_index(chrom_name, chrom_names)
        gc.append(bin_parameters[chrom_index][3, :])
    # for chrom_name
    return(gc)
# get_gc
#
#............................................................................
def get_single_cell_counts(raw_profiles, cell):
    cell_counts = []
    n_chrom = get_n_chrom(raw_profiles)
    #print('n_chrom=%d' % n_chrom)
    for chrom in xrange(n_chrom):
        tmp = raw_profiles[chrom][cell, :]
        cell_counts.append(tmp)
    # for chrom
    cell_counts = np.array(cell_counts)
    #print('cell_counts.shape=')
    #print(cell_counts.shape)
    return(cell_counts)
# get_single_cell_counts
# 
#............................................................................
def get_genome_mappability(bin_parameters, chrom_names):
    mappability = []
    for chrom_name in chrom_names:
        tmp = get_chrom_mappability(bin_parameters, chrom_name, chrom_names)
        mappability.append(tmp)
    # for chrom_name
    return(mappability)
# get_genome_mappability
#
#............................................................................
def get_genome_gc(bin_parameters, chrom_names):
    gc = []
    for chrom_name in chrom_names:
        tmp = get_chrom_gc(bin_parameters, chrom_name, chrom_names)
        gc.append(tmp)
    # for chrom_name
    return(gc)
# get_genome_gc
# 
#............................................................................
def get_n_chrom(profiles):
    return(len(profiles))
# get_n_chrom
# 
#............................................................................
def merge_bins_single_cell(profiles, cell, n_merge, average=False):
    tmp = []
    n_chrom = get_n_chrom(profiles)
    for chrom_index in xrange(n_chrom):
        tmp.append(profiles[chrom_index][cell, :])
    # for chrom_index
    return(merge_bins(tmp, n_merge, average))
# merge_bins_single_cell
#
#............................................................................
def merge_bins(profiles, n_merge, average=False):
    # TODO: handle exceptions
    if n_merge < 2:
        # TODO: warning: illegal n_merge
        return(profiles)
    # if n_merge
    n_chrom = get_n_chrom(profiles)
    if n_chrom < 1:
        # TODO: warning: no chromosomes!
        return(profiles)
    # if n_chrom
    n_cells = get_n_cells(profiles)
    if n_cells < 1:
        # TODO: warning: no cells!
        return(profiles)
    # if n_cells
    merged_profiles = [None] * n_chrom
    for chrom in range(n_chrom):
        n_bins = profiles[chrom].shape[1]
        if n_bins < n_merge:
            # TODO: warning
            pass
        # if n_bins
        #n_merged_bins = int(math.ceil(n_bins / float(n_merge)))
        tmp = np.cumsum(profiles[chrom], axis=1, dtype=float)
        tmp[:, n_merge:] = tmp[:, n_merge:] - tmp[:, :-n_merge]
        merged_profiles[chrom] = tmp[:, n_merge - 1::n_merge] # This throws away any leftover bins!
        if average:
            merged_profiles[chrom] /= n_merge
        # if average
    # for chrom
    return(merged_profiles)
# merge_bins
# 
#..........................................................................
def merge_bins_single_cell_single_chrom(profile, n_merge, average=False):
    #n_bins = profile.shape[0]
    #if n_bins < n_merge:
    #    # TODO: warning
    #    pass
    ## if n_bins
    #print('profile.shape: ', profile.shape)
    #print('profile[np.newaxis, :].shape: ', profile[np.newaxis, :].shape)
    merged_profile = merge_bins([profile[np.newaxis, :]], n_merge, average)[0].ravel()
    #print('merged_profile:')
    #print(merged_profile)
    return(merged_profile)
# merge_bins_single_cell_single_chrom

#..........................................................................
def gc_bias(gc, scale, lin, quad):
    gc0 = crdna.constants.GC_ORIGIN
    gc_gc0 = gc - gc0
    f = scale * (1.0 + lin * gc_gc0 + quad * gc_gc0 * gc_gc0)
    return(f)
# gc_bias

#..........................................................................
def extrapolated_gc_bias(gc, scale, linear, quadratic):
    trend = gc_bias(gc, scale, linear, quadratic)
    left_trend = gc_bias(crdna.constants.MIN_GC, scale, linear, quadratic)
    right_trend = gc_bias(crdna.constants.MAX_GC, scale, linear, quadratic)
    trend[gc < crdna.constants.MIN_GC] = left_trend
    trend[gc > crdna.constants.MAX_GC] = right_trend
    return(trend)
# extrapolated_gc_bias

#..........................................................................
def gauss(x, a, m0, sd):
    return(a * a * norm.pdf(x, loc=m0, scale=sd))
# gauss

#.................................................................................
def entropy(density, delta_x):
    epsilon = 1e-7
    e = -np.nansum(density * np.log(density + epsilon)) * delta_x
    return(e)
# entropy

#.................................................................................
def minimize_entropy(
        genome_profile,
        genome_merged_gc,
        cell,
        plot=True):
    #
    grid_points = np.array([
    [-9, 20],
    [-8.5, 18],
    [-8.5, 19],
    [-8.5, 20],
    [-8, 15],
    [-8, 16],
    [-8, 17],
    [-8, 18],
    [-8, 19],
    [-8, 20],
    [-7.5, 13],
    [-7.5, 14],
    [-7.5, 15],
    [-7.5, 16],
    [-7.5, 17],
    [-7.5, 18],
    [-7.5, 19],
    [-7.5, 20],
    [-7, 11],
    [-7, 12],
    [-7, 13],
    [-7, 14],
    [-7, 15],
    [-7, 16],
    [-7, 17],
    [-7, 18],
    [-7, 19],
    [-7, 20],
    [-6.5, 8],
    [-6.5, 9],
    [-6.5, 10],
    [-6.5, 11],
    [-6.5, 12],
    [-6.5, 13],
    [-6.5, 14],
    [-6.5, 15],
    [-6.5, 16],
    [-6.5, 17],
    [-6.5, 18],
    [-6.5, 19],
    [-6.5, 20],
    [-6, 6],
    [-6, 7],
    [-6, 8],
    [-6, 9],
    [-6, 10],
    [-6, 11],
    [-6, 12],
    [-6, 13],
    [-6, 14],
    [-6, 15],
    [-6, 16],
    [-6, 17],
    [-6, 18],
    [-6, 19],
    [-6, 20],
    [-5.5, 3],
    [-5.5, 4],
    [-5.5, 5],
    [-5.5, 6],
    [-5.5, 7],
    [-5.5, 8],
    [-5.5, 9],
    [-5.5, 10],
    [-5.5, 11],
    [-5.5, 12],
    [-5.5, 13],
    [-5.5, 14],
    [-5.5, 15],
    [-5.5, 16],
    [-5.5, 17],
    [-5.5, 18],
    [-5.5, 19],
    [-5.5, 20],
    [-5, 0],
    [-5, 1],
    [-5, 2],
    [-5, 3],
    [-5, 4],
    [-5, 5],
    [-5, 6],
    [-5, 7],
    [-5, 8],
    [-5, 9],
    [-5, 10],
    [-5, 11],
    [-5, 12],
    [-5, 13],
    [-5, 14],
    [-5, 15],
    [-5, 16],
    [-5, 17],
    [-5, 18],
    [-5, 19],
    [-5, 20],
    [-4.5, -2],
    [-4.5, -1],
    [-4.5, 0],
    [-4.5, 1],
    [-4.5, 2],
    [-4.5, 3],
    [-4.5, 4],
    [-4.5, 5],
    [-4.5, 6],
    [-4.5, 7],
    [-4.5, 8],
    [-4.5, 9],
    [-4.5, 10],
    [-4.5, 11],
    [-4.5, 12],
    [-4.5, 13],
    [-4.5, 14],
    [-4.5, 15],
    [-4.5, 16],
    [-4.5, 17],
    [-4.5, 18],
    [-4.5, 19],
    [-4.5, 20],
    [-4, -4],
    [-4, -3],
    [-4, -2],
    [-4, -1],
    [-4, 0],
    [-4, 1],
    [-4, 2],
    [-4, 3],
    [-4, 4],
    [-4, 5],
    [-4, 6],
    [-4, 7],
    [-4, 8],
    [-4, 9],
    [-4, 10],
    [-4, 11],
    [-4, 12],
    [-4, 13],
    [-4, 14],
    [-4, 15],
    [-4, 16],
    [-4, 17],
    [-4, 18],
    [-4, 19],
    [-4, 20],
    [-3.5, -7],
    [-3.5, -6],
    [-3.5, -5],
    [-3.5, -4],
    [-3.5, -3],
    [-3.5, -2],
    [-3.5, -1],
    [-3.5, 0],
    [-3.5, 1],
    [-3.5, 2],
    [-3.5, 3],
    [-3.5, 4],
    [-3.5, 5],
    [-3.5, 6],
    [-3.5, 7],
    [-3.5, 8],
    [-3.5, 9],
    [-3.5, 10],
    [-3.5, 11],
    [-3.5, 12],
    [-3.5, 13],
    [-3.5, 14],
    [-3.5, 15],
    [-3.5, 16],
    [-3.5, 17],
    [-3.5, 18],
    [-3.5, 19],
    [-3.5, 20],
    [-3, -9],
    [-3, -8],
    [-3, -7],
    [-3, -6],
    [-3, -5],
    [-3, -4],
    [-3, -3],
    [-3, -2],
    [-3, -1],
    [-3, 0],
    [-3, 1],
    [-3, 2],
    [-3, 3],
    [-3, 4],
    [-3, 5],
    [-3, 6],
    [-3, 7],
    [-3, 8],
    [-3, 9],
    [-3, 10],
    [-3, 11],
    [-3, 12],
    [-3, 13],
    [-3, 14],
    [-3, 15],
    [-3, 16],
    [-3, 17],
    [-3, 18],
    [-3, 19],
    [-3, 20],
    [-2.5, -12],
    [-2.5, -11],
    [-2.5, -10],
    [-2.5, -9],
    [-2.5, -8],
    [-2.5, -7],
    [-2.5, -6],
    [-2.5, -5],
    [-2.5, -4],
    [-2.5, -3],
    [-2.5, -2],
    [-2.5, -1],
    [-2.5, 0],
    [-2.5, 1],
    [-2.5, 2],
    [-2.5, 3],
    [-2.5, 4],
    [-2.5, 5],
    [-2.5, 6],
    [-2.5, 7],
    [-2.5, 8],
    [-2.5, 9],
    [-2.5, 10],
    [-2.5, 11],
    [-2.5, 12],
    [-2.5, 13],
    [-2.5, 14],
    [-2.5, 15],
    [-2.5, 16],
    [-2.5, 17],
    [-2.5, 18],
    [-2.5, 19],
    [-2.5, 20],
    [-2, -14],
    [-2, -13],
    [-2, -12],
    [-2, -11],
    [-2, -10],
    [-2, -9],
    [-2, -8],
    [-2, -7],
    [-2, -6],
    [-2, -5],
    [-2, -4],
    [-2, -3],
    [-2, -2],
    [-2, -1],
    [-2, 0],
    [-2, 1],
    [-2, 2],
    [-2, 3],
    [-2, 4],
    [-2, 5],
    [-2, 6],
    [-2, 7],
    [-2, 8],
    [-2, 9],
    [-2, 10],
    [-2, 11],
    [-2, 12],
    [-2, 13],
    [-2, 14],
    [-2, 15],
    [-2, 16],
    [-2, 17],
    [-2, 18],
    [-2, 19],
    [-2, 20],
    [-1.5, -17],
    [-1.5, -16],
    [-1.5, -15],
    [-1.5, -14],
    [-1.5, -13],
    [-1.5, -12],
    [-1.5, -11],
    [-1.5, -10],
    [-1.5, -9],
    [-1.5, -8],
    [-1.5, -7],
    [-1.5, -6],
    [-1.5, -5],
    [-1.5, -4],
    [-1.5, -3],
    [-1.5, -2],
    [-1.5, -1],
    [-1.5, 0],
    [-1.5, 1],
    [-1.5, 2],
    [-1.5, 3],
    [-1.5, 4],
    [-1.5, 5],
    [-1.5, 6],
    [-1.5, 7],
    [-1.5, 8],
    [-1.5, 9],
    [-1.5, 10],
    [-1.5, 11],
    [-1.5, 12],
    [-1.5, 13],
    [-1.5, 14],
    [-1.5, 15],
    [-1.5, 16],
    [-1.5, 17],
    [-1.5, 18],
    [-1.5, 19],
    [-1.5, 20],
    [-1, -19],
    [-1, -18],
    [-1, -17],
    [-1, -16],
    [-1, -15],
    [-1, -14],
    [-1, -13],
    [-1, -12],
    [-1, -11],
    [-1, -10],
    [-1, -9],
    [-1, -8],
    [-1, -7],
    [-1, -6],
    [-1, -5],
    [-1, -4],
    [-1, -3],
    [-1, -2],
    [-1, -1],
    [-1, 0],
    [-1, 1],
    [-1, 2],
    [-1, 3],
    [-1, 4],
    [-1, 5],
    [-1, 6],
    [-1, 7],
    [-1, 8],
    [-1, 9],
    [-1, 10],
    [-1, 11],
    [-1, 12],
    [-1, 13],
    [-1, 14],
    [-1, 15],
    [-1, 16],
    [-1, 17],
    [-1, 18],
    [-1, 19],
    [-1, 20],
    [-0.5, -22],
    [-0.5, -21],
    [-0.5, -20],
    [-0.5, -19],
    [-0.5, -18],
    [-0.5, -17],
    [-0.5, -16],
    [-0.5, -15],
    [-0.5, -14],
    [-0.5, -13],
    [-0.5, -12],
    [-0.5, -11],
    [-0.5, -10],
    [-0.5, -9],
    [-0.5, -8],
    [-0.5, -7],
    [-0.5, -6],
    [-0.5, -5],
    [-0.5, -4],
    [-0.5, -3],
    [-0.5, -2],
    [-0.5, -1],
    [-0.5, 0],
    [-0.5, 1],
    [-0.5, 2],
    [-0.5, 3],
    [-0.5, 4],
    [-0.5, 5],
    [-0.5, 6],
    [-0.5, 7],
    [-0.5, 8],
    [-0.5, 9],
    [-0.5, 10],
    [-0.5, 11],
    [-0.5, 12],
    [-0.5, 13],
    [-0.5, 14],
    [-0.5, 15],
    [-0.5, 16],
    [-0.5, 17],
    [-0.5, 18],
    [-0.5, 19],
    [-0.5, 20],
    [0, -24],
    [0, -23],
    [0, -22],
    [0, -21],
    [0, -20],
    [0, -19],
    [0, -18],
    [0, -17],
    [0, -16],
    [0, -15],
    [0, -14],
    [0, -13],
    [0, -12],
    [0, -11],
    [0, -10],
    [0, -9],
    [0, -8],
    [0, -7],
    [0, -6],
    [0, -5],
    [0, -4],
    [0, -3],
    [0, -2],
    [0, -1],
    [0, 0],
    [0, 1],
    [0, 2],
    [0, 3],
    [0, 4],
    [0, 5],
    [0, 6],
    [0, 7],
    [0, 8],
    [0, 9],
    [0, 10],
    [0, 11],
    [0, 12],
    [0, 13],
    [0, 14],
    [0, 15],
    [0, 16],
    [0, 17],
    [0, 18],
    [0, 19],
    [0, 20],
    [0.5, -27],
    [0.5, -26],
    [0.5, -25],
    [0.5, -24],
    [0.5, -23],
    [0.5, -22],
    [0.5, -21],
    [0.5, -20],
    [0.5, -19],
    [0.5, -18],
    [0.5, -17],
    [0.5, -16],
    [0.5, -15],
    [0.5, -14],
    [0.5, -13],
    [0.5, -12],
    [0.5, -11],
    [0.5, -10],
    [0.5, -9],
    [0.5, -8],
    [0.5, -7],
    [0.5, -6],
    [0.5, -5],
    [0.5, -4],
    [0.5, -3],
    [0.5, -2],
    [0.5, -1],
    [0.5, 0],
    [0.5, 1],
    [0.5, 2],
    [0.5, 3],
    [0.5, 4],
    [0.5, 5],
    [0.5, 6],
    [0.5, 7],
    [0.5, 8],
    [0.5, 9],
    [0.5, 10],
    [0.5, 11],
    [0.5, 12],
    [0.5, 13],
    [0.5, 14],
    [0.5, 15],
    [0.5, 16],
    [0.5, 17],
    [0.5, 18],
    [0.5, 19],
    [0.5, 20],
    [1, -29],
    [1, -28],
    [1, -27],
    [1, -26],
    [1, -25],
    [1, -24],
    [1, -23],
    [1, -22],
    [1, -21],
    [1, -20],
    [1, -19],
    [1, -18],
    [1, -17],
    [1, -16],
    [1, -15],
    [1, -14],
    [1, -13],
    [1, -12],
    [1, -11],
    [1, -10],
    [1, -9],
    [1, -8],
    [1, -7],
    [1, -6],
    [1, -5],
    [1, -4],
    [1, -3],
    [1, -2],
    [1, -1],
    [1, 0],
    [1, 1],
    [1, 2],
    [1, 3],
    [1, 4],
    [1, 5],
    [1, 6],
    [1, 7],
    [1, 8],
    [1, 9],
    [1, 10],
    [1, 11],
    [1, 12],
    [1, 13],
    [1, 14],
    [1, 15],
    [1, 16],
    [1, 17],
    [1, 18],
    [1, 19],
    [1, 20],
    [1.5, -32],
    [1.5, -31],
    [1.5, -30],
    [1.5, -29],
    [1.5, -28],
    [1.5, -27],
    [1.5, -26],
    [1.5, -25],
    [1.5, -24],
    [1.5, -23],
    [1.5, -22],
    [1.5, -21],
    [1.5, -20],
    [1.5, -19],
    [1.5, -18],
    [1.5, -17],
    [1.5, -16],
    [1.5, -15],
    [1.5, -14],
    [1.5, -13],
    [1.5, -12],
    [1.5, -11],
    [1.5, -10],
    [1.5, -9],
    [1.5, -8],
    [1.5, -7],
    [1.5, -6],
    [1.5, -5],
    [1.5, -4],
    [1.5, -3],
    [1.5, -2],
    [1.5, -1],
    [1.5, 0],
    [1.5, 1],
    [1.5, 2],
    [1.5, 3],
    [1.5, 4],
    [1.5, 5],
    [1.5, 6],
    [1.5, 7],
    [1.5, 8],
    [1.5, 9],
    [1.5, 10],
    [1.5, 11],
    [1.5, 12],
    [1.5, 13],
    [1.5, 14],
    [1.5, 15],
    [1.5, 16],
    [1.5, 17],
    [1.5, 18],
    [1.5, 19],
    [1.5, 20],
    [2, -31],
    [2, -30],
    [2, -29],
    [2, -28],
    [2, -27],
    [2, -26],
    [2, -25],
    [2, -24],
    [2, -23],
    [2, -22],
    [2, -21],
    [2, -20],
    [2, -19],
    [2, -18],
    [2, -17],
    [2, -16],
    [2, -15],
    [2, -14],
    [2, -13],
    [2, -12],
    [2, -11],
    [2, -10],
    [2, -9],
    [2, -8],
    [2, -7],
    [2, -6],
    [2, -5],
    [2, -4],
    [2, -3],
    [2, -2],
    [2, -1],
    [2, 0],
    [2, 1],
    [2, 2],
    [2, 3],
    [2, 4],
    [2, 5],
    [2, 6],
    [2, 7],
    [2, 8],
    [2, 9],
    [2, 10],
    [2, 11],
    [2, 12],
    [2, 13],
    [2, 14],
    [2, 15],
    [2, 16],
    [2, 17],
    [2, 18],
    [2, 19],
    [2, 20],
    [2.5, -27],
    [2.5, -26],
    [2.5, -25],
    [2.5, -24],
    [2.5, -23],
    [2.5, -22],
    [2.5, -21],
    [2.5, -20],
    [2.5, -19],
    [2.5, -18],
    [2.5, -17],
    [2.5, -16],
    [2.5, -15],
    [2.5, -14],
    [2.5, -13],
    [2.5, -12],
    [2.5, -11],
    [2.5, -10],
    [2.5, -9],
    [2.5, -8],
    [2.5, -7],
    [2.5, -6],
    [2.5, -5],
    [2.5, -4],
    [2.5, -3],
    [2.5, -2],
    [2.5, -1],
    [2.5, 0],
    [2.5, 1],
    [2.5, 2],
    [2.5, 3],
    [2.5, 4],
    [2.5, 5],
    [2.5, 6],
    [2.5, 7],
    [2.5, 8],
    [2.5, 9],
    [2.5, 10],
    [2.5, 11],
    [2.5, 12],
    [2.5, 13],
    [2.5, 14],
    [2.5, 15],
    [2.5, 16],
    [2.5, 17],
    [2.5, 18],
    [2.5, 19],
    [2.5, 20],
    [3, -24],
    [3, -23],
    [3, -22],
    [3, -21],
    [3, -20],
    [3, -19],
    [3, -18],
    [3, -17],
    [3, -16],
    [3, -15],
    [3, -14],
    [3, -13],
    [3, -12],
    [3, -11],
    [3, -10],
    [3, -9],
    [3, -8],
    [3, -7],
    [3, -6],
    [3, -5],
    [3, -4],
    [3, -3],
    [3, -2],
    [3, -1],
    [3, 0],
    [3, 1],
    [3, 2],
    [3, 3],
    [3, 4],
    [3, 5],
    [3, 6],
    [3, 7],
    [3, 8],
    [3, 9],
    [3, 10],
    [3, 11],
    [3, 12],
    [3, 13],
    [3, 14],
    [3, 15],
    [3, 16],
    [3, 17],
    [3, 18],
    [3, 19],
    [3, 20],
    [3.5, -21],
    [3.5, -20],
    [3.5, -19],
    [3.5, -18],
    [3.5, -17],
    [3.5, -16],
    [3.5, -15],
    [3.5, -14],
    [3.5, -13],
    [3.5, -12],
    [3.5, -11],
    [3.5, -10],
    [3.5, -9],
    [3.5, -8],
    [3.5, -7],
    [3.5, -6],
    [3.5, -5],
    [3.5, -4],
    [3.5, -3],
    [3.5, -2],
    [3.5, -1],
    [3.5, 0],
    [3.5, 1],
    [3.5, 2],
    [3.5, 3],
    [3.5, 4],
    [3.5, 5],
    [3.5, 6],
    [3.5, 7],
    [3.5, 8],
    [3.5, 9],
    [3.5, 10],
    [3.5, 11],
    [3.5, 12],
    [3.5, 13],
    [3.5, 14],
    [3.5, 15],
    [3.5, 16],
    [3.5, 17],
    [3.5, 18],
    [3.5, 19],
    [3.5, 20],
    [4, -17],
    [4, -16],
    [4, -15],
    [4, -14],
    [4, -13],
    [4, -12],
    [4, -11],
    [4, -10],
    [4, -9],
    [4, -8],
    [4, -7],
    [4, -6],
    [4, -5],
    [4, -4],
    [4, -3],
    [4, -2],
    [4, -1],
    [4, 0],
    [4, 1],
    [4, 2],
    [4, 3],
    [4, 4],
    [4, 5],
    [4, 6],
    [4, 7],
    [4, 8],
    [4, 9],
    [4, 10],
    [4, 11],
    [4, 12],
    [4, 13],
    [4, 14],
    [4, 15],
    [4, 16],
    [4, 17],
    [4, 18],
    [4, 19],
    [4, 20],
    [4.5, -14],
    [4.5, -13],
    [4.5, -12],
    [4.5, -11],
    [4.5, -10],
    [4.5, -9],
    [4.5, -8],
    [4.5, -7],
    [4.5, -6],
    [4.5, -5],
    [4.5, -4],
    [4.5, -3],
    [4.5, -2],
    [4.5, -1],
    [4.5, 0],
    [4.5, 1],
    [4.5, 2],
    [4.5, 3],
    [4.5, 4],
    [4.5, 5],
    [4.5, 6],
    [4.5, 7],
    [4.5, 8],
    [4.5, 9],
    [4.5, 10],
    [4.5, 11],
    [4.5, 12],
    [4.5, 13],
    [4.5, 14],
    [4.5, 15],
    [4.5, 16],
    [4.5, 17],
    [4.5, 18],
    [4.5, 19],
    [4.5, 20],
    [5, -11],
    [5, -10],
    [5, -9],
    [5, -8],
    [5, -7],
    [5, -6],
    [5, -5],
    [5, -4],
    [5, -3],
    [5, -2],
    [5, -1],
    [5, 0],
    [5, 1],
    [5, 2],
    [5, 3],
    [5, 4],
    [5, 5],
    [5, 6],
    [5, 7],
    [5, 8],
    [5, 9],
    [5, 10],
    [5, 11],
    [5, 12],
    [5, 13],
    [5, 14],
    [5, 15],
    [5, 16],
    [5, 17],
    [5, 18],
    [5, 19],
    [5, 20],
    [5.5, -7],
    [5.5, -6],
    [5.5, -5],
    [5.5, -4],
    [5.5, -3],
    [5.5, -2],
    [5.5, -1],
    [5.5, 0],
    [5.5, 1],
    [5.5, 2],
    [5.5, 3],
    [5.5, 4],
    [5.5, 5],
    [5.5, 6],
    [5.5, 7],
    [5.5, 8],
    [5.5, 9],
    [5.5, 10],
    [5.5, 11],
    [5.5, 12],
    [5.5, 13],
    [5.5, 14],
    [5.5, 15],
    [5.5, 16],
    [5.5, 17],
    [5.5, 18],
    [5.5, 19],
    [5.5, 20],
    [6, -4],
    [6, -3],
    [6, -2],
    [6, -1],
    [6, 0],
    [6, 1],
    [6, 2],
    [6, 3],
    [6, 4],
    [6, 5],
    [6, 6],
    [6, 7],
    [6, 8],
    [6, 9],
    [6, 10],
    [6, 11],
    [6, 12],
    [6, 13],
    [6, 14],
    [6, 15],
    [6, 16],
    [6, 17],
    [6, 18],
    [6, 19],
    [6, 20],
    [6.5, -1],
    [6.5, 0],
    [6.5, 1],
    [6.5, 2],
    [6.5, 3],
    [6.5, 4],
    [6.5, 5],
    [6.5, 6],
    [6.5, 7],
    [6.5, 8],
    [6.5, 9],
    [6.5, 10],
    [6.5, 11],
    [6.5, 12],
    [6.5, 13],
    [6.5, 14],
    [6.5, 15],
    [6.5, 16],
    [6.5, 17],
    [6.5, 18],
    [6.5, 19],
    [6.5, 20],
    [7, 3],
    [7, 4],
    [7, 5],
    [7, 6],
    [7, 7],
    [7, 8],
    [7, 9],
    [7, 10],
    [7, 11],
    [7, 12],
    [7, 13],
    [7, 14],
    [7, 15],
    [7, 16],
    [7, 17],
    [7, 18],
    [7, 19],
    [7, 20],
    [7.5, 6],
    [7.5, 7],
    [7.5, 8],
    [7.5, 9],
    [7.5, 10],
    [7.5, 11],
    [7.5, 12],
    [7.5, 13],
    [7.5, 14],
    [7.5, 15],
    [7.5, 16],
    [7.5, 17],
    [7.5, 18],
    [7.5, 19],
    [7.5, 20],
    [8, 9],
    [8, 10],
    [8, 11],
    [8, 12],
    [8, 13],
    [8, 14],
    [8, 15],
    [8, 16],
    [8, 17],
    [8, 18],
    [8, 19],
    [8, 20],
    [8.5, 13],
    [8.5, 14],
    [8.5, 15],
    [8.5, 16],
    [8.5, 17],
    [8.5, 18],
    [8.5, 19],
    [8.5, 20],
    [9, 16],
    [9, 17],
    [9, 18],
    [9, 19],
    [9, 20],
    [9.5, 19],
    [9.5, 20]
    ], dtype='float64')
    #
    low_gc = crdna.constants.MIN_GC 
    high_gc = crdna.constants.MAX_GC
    # TODO: may want to revisit this 700 minimum ceiling
    count_ceiling = max(700, np.percentile(genome_profile, 95))
    gc_bin_selector = which(
        (genome_merged_gc >= low_gc) & 
        (genome_merged_gc <= high_gc) &
        (genome_profile <= count_ceiling))
    counts_in_gc_bin = genome_profile[gc_bin_selector]
    gc0 = crdna.constants.GC_ORIGIN
    g = genome_merged_gc[gc_bin_selector] - gc0
    g2 = g * g
    #
    #print('frequency:')
    #print(frequency)
    #print('bin_counts:')
    #print(bin_counts)
    #
    best_entropy = np.finfo('f').max
    best_linear = 0.0
    best_quadratic = 0.0
    bins = np.linspace(0, count_ceiling, 141)
    n_grid_points = grid_points.shape[0]
    #print(grid_points)
    for i in xrange(n_grid_points):
        linear, quadratic = grid_points[i, :]
        #
        try:
            parabola = 1 + linear * g + quadratic * g2
            transformed_profile = counts_in_gc_bin / parabola
            #frequency, bin_counts, patches = plt.hist(
            #        transformed_profile, bins=bins)
            frequency, bin_counts = np.histogram(
                transformed_profile, bins=bins)
            g_step = np.nanmean(np.diff(bins))
            frequency = frequency.astype('f')
            # Must shift the frequency because of a bug in numpy:
            # When bins are not of equal length, frequencies are sometimes negative
            # (see https://github.com/numpy/numpy/issues/9222)
            baseline = np.nanmin(frequency)
            frequency -= baseline
            #print('min(frequency)=%f, max(frequency)=%f' % (
            #        np.nanmin(frequency), np.nanmax(frequency)))
            norm = np.nansum(frequency) * g_step
            frequency /= norm
            e = entropy(frequency, g_step)
            #
            #fig = plt.figure(figsize=[12, 5])
            #plt.plot(bin_counts[:-1], frequency)
            #plt.title('Linear: %f, Quadratic: %f, Entropy: %f' % (linear, quadratic, e))
            #print('frequency:')
            #print(frequency)
            #print('Norm: %f, Linear: %f, Quadratic: %f, Entropy: %f' % (norm, linear, quadratic, e))
            #
        except Exception as error:
            # TODO: fix failures (if any)
            print('minimize_entropy() caught an exception: %s' % repr(error))
            print('Exception thrown. Setting target to infinity')
            e = np.finfo('f').max
        # try/except
        if e < best_entropy:
            best_entropy = e
            best_linear = linear
            best_quadratic = quadratic
            #print('i=%d, best_linear=%d, best_quadratic=%f, best_entropy=%f' % 
            #      (i, best_linear, best_quadratic, best_entropy))
        # if target
    # for for i, (linear, quadratic)
    #
    #fig = plt.figure(figsize=[12, 5])
    #plt.scatter(g, counts_in_gc_bin)
    #x = np.linspace(0.3, 0.65, 100) - crdna.constants.GC_ORIGIN
    #plt.plot(x, 200 * (1 + best_linear * x + best_quadratic * x * x), color='r', linewidth=4)
    #
    ### Just out of curiosity:
    ##linear_fit, goodness_of_linear_fit = curve_fit(
    ##    lambda g, a, b: a + b * g, g, counts_in_gc_bin, p0=[200, 0])
    ##print('linear_fit:')
    ##print(linear_fit)
    #
    #plt.title('Linear: %f, Quadratic: %f, Entropy: %f' % (best_linear, best_quadratic, best_entropy))
    #
    result = {
        'Entropy': best_entropy,
        'Linear': best_linear,
        'Quadratic': best_quadratic
        }
    return(result)
# minimize_entropy

#.................................................................................
#def normalize_profile(merged_cell_counts, merged_gc, results):
#    normalized_cell_counts = []
#    n_chrom = get_n_chrom(merged_cell_counts)
#    scale = (results['NormalizationResults']['InitialNormalizationParameters']['Scale'] *
#             results['NormalizationResults']['RenormalizationParameters']['Scale'])
#    lin =   (results['NormalizationResults']['InitialNormalizationParameters']['Linear'] +
#             results['NormalizationResults']['RenormalizationParameters']['Linear'])
#    quad =  (results['NormalizationResults']['InitialNormalizationParameters']['Quadratic'] +
#             results['NormalizationResults']['RenormalizationParameters']['Quadratic'])
#    for chrom_index in xrange(n_chrom):
#        chr_profile = (merged_cell_counts[chrom_index] /
#            gc_bias(merged_gc[chrom_index], scale, lin, quad))
#        chr_profile[chr_profile < 0] = np.nan
#        normalized_cell_counts.append(chr_profile)
#    # for chrom_index
#    return(normalized_cell_counts)
## normalize_profile

#.................................................................................
def estimate_n_merge(raw_cell_counts, target_bin_count=200.0, confident_genome_fraction=1.0, 
                     fixed_n_merge=50, fix_n_merge=True):
    if fix_n_merge:
        return(fixed_n_merge)
    # if fix_n_merge
    mean_read_count = np.nanmean(np.concatenate(raw_cell_counts))
    n_merge = round(confident_genome_fraction * target_bin_count / mean_read_count)
    # Special case: high coverage depth, already exceeds desired target_bin_count
    n_merge = np.nanmax([n_merge, 1])
    return(n_merge)
# estimate_n_merge

#.................................................................................
def apply_normalization_parameters_single_cell(
        raw_profiles, cell, 
        mask, gc, 
        scale, linear, quadratic,
        merge_bins=False,
        n_merge=None,
        target_bin_count=200.0, 
        confident_genome_fraction=1.0):
    #
    raw_cell_counts = get_single_cell_counts(raw_profiles, cell)
    n_chrom = get_n_chrom(raw_profiles) # TODO: enable sex chromosomes, any reference
    if merge_bins:
        if n_merge is None:
            # Determine bin resolution
            n_merge = estimate_n_merge(
                raw_cell_counts, target_bin_count, confident_genome_fraction)
        # if n_merge
        #
        # Merge bins
        merged_cell_counts = merge_bins(
            [np.array([raw_profiles[chrom][cell, :]])
            for chrom in range(n_chrom)], n_merge)
        #print(len(merged_cell_counts))
        merged_gc = merge_bins([np.array([gc[chrom]]) 
            for chrom in range(n_chrom)], n_merge)
        merged_gc = [merged_gc[chrom_index][0] / n_merge 
            for chrom_index in range(n_chrom)]
        normalized_cell_counts = [merged_cell_counts[chrom_index].ravel()
            for chrom_index in range(n_chrom)]
    else:
        n_merge = 1
        normalized_cell_counts = raw_cell_counts
        merged_gc = gc
    # if merge_bins else
    #
    # Normalize and scale
    for chrom_index in xrange(n_chrom):
        #print('chrom_index=%d' % chrom_index)
        #print('normalized bins: %d' % len(normalized_cell_counts[chrom_index]))
        #print('gc bins: %d' % len(merged_gc[chrom_index]))
        assert(len(normalized_cell_counts[chrom_index]) == len(merged_gc[chrom_index]))
        #
        denominator = extrapolated_gc_bias(
            merged_gc[chrom_index], n_merge, linear, quadratic)
        normalized_cell_counts[chrom_index] /= denominator
        #
        # Set negative values to zero:
        normalized_cell_counts[chrom_index][
            normalized_cell_counts[chrom_index] < 0] = 0.0
    # for chrom_index
    #
    return(normalized_cell_counts)
# apply_normalization_parameters_single_cell

#.................................................................................
def add_profile(normalized_profiles, cell, normalized_cell_counts):
    n_chrom = get_n_chrom(normalized_profiles)
    for chrom_index in xrange(n_chrom):
        normalized_profiles[chrom_index][cell, :] = normalized_cell_counts[chrom_index]
    # for chrom
# add_profile

#.................................................................................
def generate_normalized_profiles(
        raw_profiles, mask, bin_parameters, 
        scale, linear, quadratic, chrom_names):
    n_cells = get_n_cells(raw_profiles)
    assert(raw_profiles[0].shape[0] == n_cells), "num profiles != num cells"
    assert(len(scale) == n_cells), "need scale for every cell"
    assert(len(linear) == n_cells), "need GC params for every cell"
    assert(len(quadratic) == n_cells), "need GC params for every cell"
    normalized_profiles = raw_profiles
    gc = get_gc(bin_parameters, chrom_names)
    for cell in xrange(n_cells):
        normalized_cell_counts = apply_normalization_parameters_single_cell(
            raw_profiles, cell,
            mask, gc,
            scale[cell], linear[cell], quadratic[cell],
            merge_bins=False)
        add_profile(normalized_profiles, cell, normalized_cell_counts)
    # for cell
    return(normalized_profiles)
# generate_normalized_profiles

#.................................................................................
def estimate_normalization_parameters_single_cell(
        genome_profile,
        genome_merged_gc,
        cell,
        mask,
        method='entropy'):
    scale = 1.0
    linear = 0.0
    quadratic = 0.0
    delta_gc_cv_filtered = 0.0
    if method == 'entropy':
        normalization_parameters = minimize_entropy(
            genome_profile,
            genome_merged_gc,
            cell)
        linear = normalization_parameters['Linear']
        quadratic = normalization_parameters['Quadratic']
        go = crdna.constants.GC_ORIGIN
        g = genome_merged_gc - go
        parabola = np.clip(1 + g * linear + g * g * quadratic, 0.1, 10.0)
        transformed_profile = genome_profile / parabola
        #
        delta_gc_cv_filtered = get_delta_gc_cv_filtered(
            raw_profile=genome_profile,
            normalized_profile=transformed_profile,
            genome_merged_gc=genome_merged_gc,
            mappability_mask=mask)
        #
        if delta_gc_cv_filtered < 0:
            linear = 0.0
            quadratic = 0.0
            delta_gc_cv_filtered = 0.0
        # if delta_gc_cv_filtered
        #
        # Use a robust procedure to find scale:
        try:
            scale = optimize_scale(transformed_profile)
        except Exception as error:
            # TODO: Understand why exception is being thrown,
            # handle it properly.
            print('estimate_normalization_parameters_single_cell() caught an error: %s' %
                repr(error))
            scale = 1.0 # 0.5 * np.nanmean(transformed_profile)
        # try/except
    else:
        print('Unknown method: %s' % method)
    # if use else
    #
    results = {
        'Scale': scale,
        'Linear': linear,
        'Quadratic': quadratic,
        'DeltaGCCVFiltered': delta_gc_cv_filtered}
    return(results)
# estimate_normalization_parameters_single_cell

#.................................................................................
def get_cv(profile):
    std = np.nanstd(profile)
    mean = np.nanmean(profile)
    if mean > 0:
        return(std / mean)
    else:
        return(np.nan)
    # if mean else
# get_cv

#.................................................................................
def get_delta_gc(raw_cv, norm_cv):
    tmp = raw_cv * raw_cv - norm_cv * norm_cv
    delta_gc = np.sign(tmp) * np.sqrt(np.abs(tmp))
    return(delta_gc)
# get_delta_gc

#.................................................................................
def get_delta_gc_cv_filtered(
    raw_profile,
    normalized_profile,
    genome_merged_gc,
    mappability_mask):
    gc_mask = np.logical_and(
        genome_merged_gc >= crdna.constants.MIN_GC,
        genome_merged_gc <= crdna.constants.MAX_GC)
    final_mask = np.logical_and(gc_mask, mappability_mask)
    #
    raw_cv = get_cv(raw_profile[final_mask])
    norm_cv = get_cv(normalized_profile[final_mask])
    delta_gc_cv = get_delta_gc(raw_cv, norm_cv)
    return(delta_gc_cv)
    #
# get_delta_gc_cv_filtered

#.................................................................................
def process_single_cell(
        raw_profiles,
        gc,
        cell,
        mask,
        target_bin_count=200,
        n_chrom=22,
        confident_genome_fraction=1.0
        ):
    # Determine bin resolution
    raw_cell_counts = get_single_cell_counts(raw_profiles, cell)
    #
    # Merge bins
    n_merge = estimate_n_merge(raw_cell_counts, target_bin_count, confident_genome_fraction)
    merged_cell_counts = merge_bins(
        [np.array([raw_profiles[chrom][cell, :]]) 
        for chrom in range(n_chrom)], n_merge)
    #print(len(merged_cell_counts))
    for chrom_index in xrange(n_chrom):
        # print(merged_cell_counts[chrom_index].shape)
        merged_cell_counts[chrom_index] = merged_cell_counts[chrom_index].ravel()
    # for chrom_index
    #
    # Normalize and scale
    genome_profile = np.concatenate(merged_cell_counts)
    #
    merged_gc = merge_bins([np.array([gc[chrom]]) for chrom in range(n_chrom)], n_merge)
    merged_gc = [merged_gc[chrom_index][0] / n_merge for chrom_index in range(n_chrom)]
    genome_merged_gc = np.concatenate(merged_gc)
    #
    merged_mask = merge_bins([np.array([mask[chrom].astype(int)]) for chrom in range(n_chrom)], n_merge)
    merged_mask = [merged_mask[chrom_index][0] / n_merge for chrom_index in range(n_chrom)]
    merged_mask = np.concatenate(merged_mask).astype(bool)
    #
    results = estimate_normalization_parameters_single_cell(
        genome_profile,
        genome_merged_gc, 
        cell,
        merged_mask)
    # Report results
    results['Scale'] /= float(n_merge)
    #
    return(results)
# process_single_cell

#.................................................................................
def get_mean_bin_counts(raw_profiles, cell):
    genome_profile = np.hstack(get_single_cell_counts(raw_profiles, cell))
    mean_counts = np.nanmean(genome_profile)
    return(mean_counts)
# get_mean_bin_counts

#.................................................................................
def get_n_cells(profiles):
    n_cells = profiles[0].shape[0]
    return(n_cells)
# get_n_cells

#.................................................................................
def process_profiles(raw_profiles, mask, chrom_names, bin_parameters):
    gc = get_genome_gc(bin_parameters, chrom_names)
    #    
    target_bin_count = 200
    n_cells = get_n_cells(raw_profiles)
    scale = [1.0] * n_cells
    linear = [0.0] * n_cells
    quadratic = [0.0] * n_cells
    for cell in xrange(n_cells):
        results = process_single_cell(
            raw_profiles,
            gc,
            cell,
            mask,
            target_bin_count=target_bin_count,
            n_chrom=len(chrom_names))
        scale[cell] = results['Scale']
        linear[cell] = results['Linear']
        quadratic[cell] = results['Quadratic']
        print('Completed cell %d/%d: scale=%.1f, linear=%.1f, quadratic=%.1f' % 
              (cell, n_cells, scale[cell], linear[cell], quadratic[cell]))
    # for cell
    return((scale, linear, quadratic))
# process_profiles

#.................................................................................
def optimize_scale(profile, max_ploidy=6, multiplier=3, n_scale_steps=10):
    count_frequency = se.get_count_frequency(profile)   
    count_bins = count_frequency['Count']
    frequency = count_frequency['Frequency']
    frequency[:10] = 0 # Must remove the peak at zero 
    max_counts = count_bins[np.argmax(frequency)]
    delta = multiplier * np.sqrt(max_counts)
    low_counts = max_counts - delta
    high_counts = max_counts + delta
    ploidies = np.arange(2, max_ploidy + 1)
    #
    best_target = np.finfo('f').max
    best_scale = 0.5 * max_counts
    #best_ploidy = 2
    for ploidy in ploidies:
        min_scale = low_counts / ploidy
        max_scale = high_counts / ploidy
        for scale in np.linspace(min_scale, max_scale, n_scale_steps):
            #target = evaluate_scale_target(scale, profile)
            target = gaussian_scale_target(scale, profile, count_bins, frequency)
            if target < best_target:
                best_target = target
                best_scale = scale
                #best_ploidy = ploidy
                #print('best_scale=%f, best_target=%f, best_ploidy=%d' % 
                #      (best_scale, best_target, best_ploidy))
            # if target
        # for scale
    # for ploidy
    #
    #fig = plt.figure(figsize=[12, 5])
    #plt.plot(count_bins, frequency, linewidth=2)
    #plt.gca().axvline(x=low_counts, color='r')
    #plt.gca().axvline(x=high_counts, color='r')
    #for k in range(6):
    #    plt.gca().axvline(x=(k * best_scale), color='k')
    ## for k
    #
    best_scale = double_check_scale(best_scale, profile)
    #
    return(best_scale)
# optimize_scale

#.................................................................................
def gaussian_scale_target(scale, profile, count_bins, frequency, min_bin_count=5):
    target = np.finfo('f').max
    # Classify bins:
    bin_class = np.round(profile / scale).astype('int')
    # 
    n_bins = count_bins.shape[0]
    gaussian_distribution = np.zeros(n_bins)
    n_classes = np.max(bin_class)
    for class_i in np.arange(1, n_classes + 1):
        selector_i = np.where(bin_class == class_i)[0]
        n_i = float(len(selector_i))
        if n_i > min_bin_count:
            mu_i = scale * class_i
            #sigma_i = np.sqrt(mu_i)
            sigma_i = np.nanstd(profile[selector_i])
            gaussian_distribution = gaussian_distribution + gauss(
                count_bins, np.sqrt(n_i), mu_i, sigma_i)
            #print('n_i=%f, mu_i=%f, sigma_i=%f' % (n_i, mu_i, sigma_i))
        # if n_i
    # for class
    #
    # Evaluate target
    tmp = frequency - gaussian_distribution
    tmp *= tmp
    target = np.nansum(tmp)
    #
    #fig = plt.figure(figsize=[12, 5])
    #plt.plot(count_bins, frequency, color='b')
    #plt.plot(count_bins, gaussian_distribution, color='r')
    #plt.title('scale=%f, target=%f' % (scale, target))
    #
    return(target)
# gaussian_scale_target

#.................................................................................
def double_check_scale(scale, profile, many_more=100):
    multiplier = 1.0
    scaled_profile = profile / scale
    wholes = scaled_profile.astype(float) / (1.0 + np.round(scaled_profile))
    n = len(wholes)
    #
    # Double check:
    by_twos = (wholes % 2.0).astype(int)
    evens = len(which(by_twos == 0))
    odds = n - evens
    if evens > (odds * many_more):
        if (np.round(np.nanmean(scaled_profile) / 2.0) > 1):
            multiplier *= 2.0
        # if round
    # if evens
    scale *= multiplier
    #
    return(scale)
# double_check_scale

#.................................................................................
def triple_check_scale(scale, profile, many_more=100):
    multiplier = 1.0
    scaled_profile = profile / scale
    wholes = scaled_profile.astype(float) / (1.0 + np.round(scaled_profile))
    n = len(wholes)
    #
    # Triple check:
    by_threes = (wholes % 3.0).astype(int)
    thirds = len(which(by_threes == 0))
    feconds = n - thirds
    if thirds > (feconds * many_more):
        multiplier *= 3.0
    # if thirds
    scale *= multiplier
    #
    return(scale)
# triple_check_scale


#.................................................................................
