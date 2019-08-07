#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Compute pre-/post-normalization metrics on variation in depth of coverage 
#
import pandas as pd
import numpy as np
import shutil
import scipy.stats

import crdna.profiles_data as crdna_profiles
from longranger.cnv.copy_number_tools import get_dpcvs, get_local_dpcvs, aggregate_matrix
import martian
from longranger.cnv import coverage_matrix

__MRO__ = """
stage REPORT_NORMALIZATION(
    in  string   reference_path,
    in  h5       raw_profiles,
    in  h5       normalized_profiles,
    in  int      normalization_window_size,
    in  bool     is_singlecell,
    out h5       normalized_profiles,
    src py       "stages/reporter/report_normalization",
) split using () 
"""

#
# when calculating local_dpcv, rebin dpcv coverage matrix by this factor
#
LOCAL_DPCV_REBINNING_FACTOR = 10
#
# when calculating mapd, use this window_size
# This is different from our standard window size because there
# is a precedent for comparing at 500kb
#
MAPD_WINDOW_SIZE = 500000

#
# Currently this code is very fast and doesn't require parallelization
# All of the profiles (raw/normalized) are read into memory but that's not a huge hit

def split(args):
    try:
        mat_size_gb = coverage_matrix.get_genome_matrix_size_gb(args.raw_profiles)
    except KeyError:
        mat_size_gb = 1.0
    mem_gb = int(np.ceil(3*mat_size_gb + 1))
    return {'chunks': [], 'join': {'__mem_gb': mem_gb}}

def join(args, outs, chunk_defs, chunk_outs):
    args.coerce_strings()
    outs.coerce_strings()

    raw_dpcv, raw_local_dpcv, raw_cv, raw_mu, raw_vpd, raw_local_cv, raw_profiles = \
        calculate_dpcv_for_profiles( args.raw_profiles,
            args.reference_path, args.normalization_window_size )
    if args.is_singlecell:
        ncells = raw_profiles.get_num_cells()
    else:
        ncells = (raw_profiles.get_num_cells()+1)/2
    raw_mapd, raw_dimapd = calculate_mapd_for_profiles( args.raw_profiles,
            args.reference_path, MAPD_WINDOW_SIZE, ncells )

    norm_dpcv, norm_local_dpcv, norm_cv, norm_mu, norm_vpd, norm_local_cv, _ = \
        calculate_dpcv_for_profiles( args.normalized_profiles,
            args.reference_path, args.normalization_window_size )
    norm_mapd, norm_dimapd = calculate_mapd_for_profiles( args.normalized_profiles,
            args.reference_path, MAPD_WINDOW_SIZE, ncells )

    # set min dimapd = 0
    norm_dimapd = np.clip(norm_dimapd, 0, None)
    raw_dimapd = np.clip(raw_dimapd, 0, None)

    high_dimapd_cells = define_high_dimapd_cells(norm_dimapd)

    def delta_gc(raw,norm):
        return np.sign(raw-norm) * np.sqrt(np.abs(raw*raw - norm*norm))
    # delta_gc

    delta_gc_dpcv = delta_gc(raw_dpcv,norm_dpcv)
    delta_gc_local_dpcv = delta_gc(raw_local_dpcv,norm_local_dpcv)
    delta_gc_cv = delta_gc(raw_cv, norm_cv)
    #
    # percent CV change is implemented with the opposite sign of delta_gc_dpcv,
    # delta_gc_local_dpcv
    #`
    pct_cv_change = 100.0 * np.where(raw_cv>0,(norm_cv-raw_cv)/raw_cv,0.0)
    pct_vpd_change = 100.0 * np.where(raw_vpd>0,(norm_vpd-raw_vpd)/raw_vpd,0.0)

    d = pd.DataFrame( {\
        "raw_dpcv": raw_dpcv,
        "raw_local_dpcv": raw_local_dpcv,
        "raw_cv": raw_cv,
        "raw_mu": raw_mu,
        "raw_local_cv": raw_local_cv,
        "raw_vpd": raw_vpd,
        "raw_mapd": raw_mapd,
        "raw_dimapd": raw_dimapd,
        "norm_dpcv": norm_dpcv,
        "norm_local_dpcv": norm_local_dpcv,
        "norm_cv": norm_cv,
        "norm_mu": norm_mu,
        "norm_local_cv": norm_local_cv,
        "norm_vpd": norm_vpd,
        "norm_mapd": norm_mapd,
        "norm_dimapd": norm_dimapd,
        "is_high_dimapd": high_dimapd_cells,
        "delta_gc_dpcv": delta_gc_dpcv,
        "delta_gc_local_dpcv": delta_gc_local_dpcv,
        "delta_gc_cv": delta_gc_cv,
        "pct_cv_change": pct_cv_change,
        "pct_vpd_change": pct_vpd_change,
        }, 
        columns = ['raw_dpcv','norm_dpcv','delta_gc_dpcv',
        'raw_local_dpcv','norm_local_dpcv','delta_gc_local_dpcv',
        'raw_cv','norm_cv','raw_mu','norm_mu','pct_cv_change','delta_gc_cv',
        'raw_local_cv','norm_local_cv',
        'raw_vpd','norm_vpd','pct_vpd_change',
        'raw_mapd','norm_mapd','raw_dimapd','norm_dimapd', 'is_high_dimapd']
    )

    # d.to_csv( outs.normalization_metrics, index=False, na_rep="nan" )

    shutil.copyfile(args.normalized_profiles, outs.normalized_profiles)
    norm_store = pd.HDFStore(outs.normalized_profiles, mode='a')
    norm_store['/normalization_metrics'] = d
    norm_store.close()

    # logging
    print "# normalization metrics"
    print "value, mean, min, p25, p50, p75, max"
    for col in d.keys():
        if len(d[col]) > 0 and not np.isnan(d[col]).all():
            mean = np.nanmean(d[col])
            vals = np.nanpercentile(d[col], [0, 25, 50, 75, 100])
        else:
            mean = float('nan')
            vals = [float('nan')]*5
        print ("{}, " + ", ".join(["{:.3g}"]*6)).format(col, mean, *vals)

def define_high_dimapd_cells(dimapd):
    """ Label cells with high dimapd. """
    def noisy_filter(dimapd, mu, sig):
        dev_pvalue = scipy.stats.norm.sf(dimapd, mu, sig)
        mapd_noisy = dev_pvalue < 0.01
        return mapd_noisy

    ## use MAD as initial estimate of std
    mmad = np.median(np.abs(dimapd - np.median(dimapd)))
    sig1 = 1.4826*mmad
    mu1 = np.median(dimapd)
    mapd_noisy = noisy_filter(dimapd, mu1, sig1)

    ## filter outliers and re-estimate mu, sigma and outliers
    mu2 = dimapd[~mapd_noisy].mean()
    sig2 = dimapd[~mapd_noisy].std()
    mapd_noisy2 = noisy_filter(dimapd, mu2, sig2)
    return mapd_noisy2.astype(int)

def _read_stitched_coverage( profiles_h5, reference_path, window_size, mask_data=None ):
    #
    # load profiles into memory
    #
    profiles = crdna_profiles.ProfilesData2( 
        profiles_h5, reference_path, load_conf_filter=False, reuse_mask_from=mask_data )
    #
    # apply mappability mask
    #
    nbins, n_unmasked = profiles.apply_mask(use_default_mask=True,use_conf_filter=False)
    martian.log_info("%s: (%d/%d) unmasked bins"%(profiles_h5,n_unmasked,nbins))
    #
    # rebin data to requested window size
    #
    norm_factor = window_size / profiles.get_window_size()
    profiles.aggregate(norm_factor)
    #
    # calculate dpcv and local_dpcv
    # only include autosomes
    #
    _, _, coverage = profiles.get_stitched_coverage(allow_sex_chromosomes=False)

    return coverage, profiles

def calculate_dpcv_for_profiles( profiles_h5, reference_path, normalization_window_size, mask_data=None ):
    """Returns per cell dpcv,local_dpcv,cv,mean,vpd,ProfilesData2 for the profiles in profiles.h5
   
    Statistical quantities are calculated using a binning of normalization_window_size

    If mask_data is supplied then the mask is lifted over from the supplied object and not
    read from profiles_h5
    """
    coverage, profiles = _read_stitched_coverage(profiles_h5,reference_path,
        normalization_window_size,mask_data)
    
    dpcv_per_cell = get_dpcvs(coverage)
    #
    # calculate CV and mean
    #
    cv, mu = get_cv_and_mean(coverage)
    martian.log_info("%s: (cv,mu)=(%f,%f)"%(profiles_h5,cv.mean(),mu.mean()))
    #
    # calculate VPD
    #
    vpd = get_vpd(coverage)
    #
    # calculate local_dpcv
    #
    G_coarse = aggregate_matrix( coverage, LOCAL_DPCV_REBINNING_FACTOR )
    _, local_dpcv_per_cell = get_local_dpcvs( G_coarse, dpcv_per_cell, LOCAL_DPCV_REBINNING_FACTOR )
    
    local_cv = get_local_cv( coverage, mu )

    return dpcv_per_cell, local_dpcv_per_cell, cv, mu, vpd, local_cv, profiles
# calculate_dpcv_for_profiles

def calculate_mapd_for_profiles(profiles_h5, reference_path, normalization_window_size, ncells, mask_data=None):
    """Returns per cell MAPD for the profiles in profiles.h5

    Statistical quantities are calculated using a binning of normalization_window_size

    If mask_data is supplied then the mask is lifted over from the supplied object and not
    read from profiles_h5
    """
    coverage, profiles = _read_stitched_coverage(profiles_h5,reference_path,
        normalization_window_size,mask_data)
    cv, mu = get_cv_and_mean(coverage)
    martian.log_info("%s: (window_size,cv,mu)=(%d,%f,%f)"%\
        (profiles_h5,normalization_window_size,cv.mean(),mu.mean()))
    mapd = get_mapd(coverage,mu)
    dimapd = get_dimapd(mapd,mu,ncells)

    return mapd, dimapd

#
#..............................................................................
def get_local_cv( counts, mu ):
    scaled_profiles = np.divide(counts, np.expand_dims(mu, 1))
    scaled_diff = np.diff(scaled_profiles, axis=1)
    local_cv = scaled_diff.std(axis=1)
    return( local_cv )

def get_cv_and_mean( counts ):
    sigma = counts.std(axis=1)
    mu = counts.mean(axis=1)
    cv = np.where(mu==0.0,0.0,sigma/mu)
    return cv, mu

def _interpolated_median( counts ):
    """This version of an interpolated median is specifically designed for a integer-spaced
    counts matrix (ncells x nbins)
    
    The definition of interpolated median can be found here: https://en.wikipedia.org/wiki/Median
    :returns: np.array of interpolated median for each cell"""
    result = np.zeros(counts.shape[0])
    # n.b. w is the width of the discretization interval
    # if we estimate this from the data it would look like 0.5, but in truth the vast majority
    # is discretized at 1.0 and we should use that
    w = 1.0
    meds = np.median(counts,axis=1)
    d = (np.sign(counts-meds[:,np.newaxis])+1).astype('int')
    for i in xrange(counts.shape[0]):
        e = np.bincount(d[i])
        if len(e)!=3:
            result[i] = meds[i]
        else:
            result[i] = meds[i] + w/2.0 * float(e[2]-e[0])/float(e[1]) if e[1] > 0 else meds[i]
    return result

def get_mapd( counts, mu ):
    """MAPD is median of absolute deviation of pairwise differences.
    Uses interpolated median for more robust implementation of
    median when dealing with highly discretized data (like
    low count data)"""
    diff_counts = np.diff(counts,axis=1)
    median_diff_counts = np.median(diff_counts,axis=1)[:,np.newaxis]
    return np.where(mu==0,0.0,_interpolated_median(np.abs(diff_counts-median_diff_counts))/mu)

def _remove_linear_trend( x, y, frac_discard=0.05 ):
    """Removes a linear trend from y, using a robust regression technique:
    1) fit a linear trend to x,y
    2) discard all points with absolute residuals >= (1-frac_discard) percentile
    3) fit an improved linear trend to x',y'
    4) remove the linear trend from y, leaving the mean intact
    """
    if len(x)<5:
        return y
    c = np.polyfit(x,y,1)
    y2 = c[1] + c[0]*x
    resid = y2-y
    resid_max = np.percentile(np.abs(resid),100.0*(1.0-frac_discard))
    i = np.where( np.abs(resid) < resid_max )[0]
    if len(i)==0:
        return y
    c2 = np.polyfit(x[i],y[i],1)
    return y - c2[0]*x

def get_dimapd( mapd, mu, ncells ):
    """Returns a depth-independent version of MAPD using the property that 
    MAPD*sqrt(mu) is relatively constant across depth (with a small linear trend
    attributed to uncaught PCR duplications or other background exponential noise
    
    Until we can correct this trend, nodes behave differently than single cells and
    DIMAPD is only calculated for single cells.
    """
    result = np.zeros(len(mapd))
    result[:ncells] = _remove_linear_trend(mu[:ncells],(mapd*np.sqrt(mu))[:ncells],frac_discard=0.2)
    result[ncells:] = np.nan
    return result

def get_vpd_trimmed( counts ):
    """
    VPD is "variance of pairwise differences", a customer-facing
    evenness metric intended to have invariance to CNVs.
    
    It is defined as (var(diff(Y))-2*mean(Y))/(2*mean(Y)^2)

    To reduce impact of CNVs, this function trims away spikes 
    (both positive and negative) from the differential profile 
    and evaluates vpd on sequential differences within segments.
    The trimming is done on the highest and lowest 5% data points,
    corresponding to CNV start and end breakpoints, respectively.

    The trimmed version of VPD is currently not included in the 
    normalization metric table.
    """
    diff_counts = (np.diff(counts,axis=1))
    # Trim extreme diff_counts values:
    n_trim = int(np.round(0.05 * diff_counts.shape[1]))
    diff_counts = np.sort(diff_counts, axis=1)[:, n_trim:-n_trim]
    var_pd = diff_counts.var(axis=1)
    mu = counts.mean(axis=1)
    vpd = (var_pd - 2.0 * mu) / (2.0 * mu * mu)
    return vpd
# get_vpd_trimmed
#
#..............................................................................
def get_vpd( counts ):
    """VPD is "variance of pairwise differences", a customer-facing
    evenness metric intended to have invariance to CNVs.
    
    It is defined as (var(diff(Y))-2*mean(Y))/(2*mean(Y)^2)"""
    var_pd = (np.diff(counts,axis=1)).var(axis=1)
    mu = counts.mean(axis=1)
    vpd = (var_pd - 2.0 * mu) / (2.0 * mu * mu)
    return vpd
# get_vpd
#
#..............................................................................
