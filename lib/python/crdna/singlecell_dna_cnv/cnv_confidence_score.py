#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import math
import numpy as np
import crdna.constants
import longranger.cnv.contig_manager as contig_manager

def get_segment_scores(raw_profile, poisson_expectations, mask,
        start, end, ploidy):
    """Calculate the per-bin confidence scores (log(posterior)) for the provided
    profile under the hypothesis that true ploidy==ploidy

    Bins where mask==False are provided with a score of 0.0
    """
    observed_count = raw_profile[start:end]
    observed_count = observed_count[ mask[start:end] ]
    if len(observed_count)==0:
        return np.zeros(end-start,dtype='float32')
    count = np.round(observed_count).astype(int)
    mu = poisson_expectations[start:end]
    mu = mu[ mask[start:end] ]
    ## ploidy zero emission rate
    EPSILON = 1e-2

    if ploidy == 0:
        expectation_tested = EPSILON * mu
        expectation_above = mu
        lm = np.log(mu)
        let = math.log(EPSILON) + lm
        lea = lm
        num = -expectation_tested + count*let
        denom = np.logaddexp(num, -expectation_above+count*lea)
    else:
        expectation_below, expectation_tested, expectation_above = [(ploidy + a) * mu for a in (-1, 0, 1)]
        lm = np.log(mu)
        if ploidy == 1:
            expectation_below = EPSILON * mu
            leb = math.log(EPSILON) + lm
            let, lea = [ math.log(float(ploidy+a)) + lm for a in (0,1) ]
        else:
            leb, let, lea = [ math.log(float(ploidy+a)) + lm for a in (-1,0,1) ]
        num = -expectation_tested + count*let
        denom = np.logaddexp(-expectation_below+count*leb,
            np.logaddexp(num, -expectation_above+count*lea))
    q_values = - np.log1p(-np.exp(num-denom)) / math.log(10.0)
    #
    # reinflate to all bins
    #
    q_values_all = np.zeros(end-start, dtype='float32')
    q_values_all[ mask[start:end] ] = q_values
    return q_values_all

def calculate_logposterior_matrix( raw_profiles, poisson_expectations, mask, 
        cnv_calls, reference_path, bin_size ):
    """Create a ncell x nbins matrix of log(posterior) values
    """
    ref = contig_manager.contig_manager(reference_path)
    chrom_names = ref.primary_contigs(allow_sex_chromosomes=True)

    logp = []
    ncells = 0
    for chrom_index in xrange(len(raw_profiles)):
        if ncells==0:
            ncells = raw_profiles[chrom_index].shape[0]
        nbins = raw_profiles[chrom_index].shape[1]
        logp.append( np.zeros( ( ncells, nbins ), dtype='float32' ) )

    for cnv_call in cnv_calls.itertuples():
        #
        # the create cnv tracks module already sets confidence to zero for masked bins
        # just use that confidence if it's already set
        #
        if cnv_call.Confidence==0.0:
            continue
        chrom_name = cnv_call.Chr
        chrom_index = chrom_names.index(chrom_name)
        ploidy = int(round(cnv_call.CopyNumber))
        assert(ploidy >= 0), 'Negative ploidy: %s' % repr(cnv_call)
        start = int(round(cnv_call.Start / bin_size))
        # end in BED file is exclusive
        end = int(round(cnv_call.End / bin_size))
        cell = cnv_call.NodeID
        #
        start = max(0, start)
        n_bins = raw_profiles[chrom_index].shape[1]
        end = min(end, n_bins)
        scores = get_segment_scores(raw_profiles[chrom_index][cell,:],
            poisson_expectations[chrom_index][cell,:],
            mask[chrom_index], start, end, ploidy)
        logp[chrom_index][cell,start:end] = scores
    # for cnv_call
    return logp

def estimate_cnv_confidence_score_v2( raw_profiles, cnv_calls, reference_path, logp, bin_size ):
    """
    Calculates a CNV confidence score (log(posterior)) for each CNV call using the pre-computed
    logp matrix of per-bin confidence scores.
    """
    ref = contig_manager.contig_manager(reference_path)
    chrom_names = ref.primary_contigs(allow_sex_chromosomes=True)
    PER_BIN_MAX_SCORE = 100.0
    scores = np.zeros( len(cnv_calls), dtype='int32' )

    for i, cnv_call in enumerate(cnv_calls.itertuples()):
        #
        # the create cnv tracks module already sets confidence to zero for masked bins
        # just use that confidence if it's already set
        #
        if cnv_call.Confidence==0.0:
            continue
        chrom_name = cnv_call.Chr
        chrom_index = chrom_names.index(chrom_name)
        start = int(round(cnv_call.Start / bin_size))
        # end in BED file is exclusive
        end = int(round(cnv_call.End / bin_size))
        cell = cnv_call.NodeID
        #
        start = max([0, start])
        n_bins = raw_profiles[chrom_index].shape[1]
        end = min([end, n_bins])
        # start can == end in the case where a CNV call happens only on the terminal bin
        # this was found by randomly breaking up a reference such that a mappable bin is
        # cut in two. This is extremely unlikely in a 'real' reference since bins at the
        # end of contigs will unmappable and/or have the same ploidy as the neighboring
        # bins. The pipeline steps that will lead to this are in
        # CREATE_CNV_TRACKS_AND_BED which converts cluster_data.h5 into cnv_calls.bed
        if start == end:  # special case: CNV call on single bin
            score = np.nansum( logp[chrom_index][cell,start] )
            score = np.clip(score, 0, PER_BIN_MAX_SCORE)
            scores[i] = min(np.round(score * 100), np.iinfo("uint8").max)
        else:
            score = np.nansum( logp[chrom_index][cell,start:end] )
            score = np.clip(score, 0, PER_BIN_MAX_SCORE*(end-start))
            scores[i] = min(np.round(score/(end-start)*100), np.iinfo("uint8").max)
    cnv_calls['Confidence'] = scores
    return cnv_calls

#............................................................................
def gc_bias(scale, linear, quadratic, gc_gc0, gc_gc02):
    return scale * np.clip(1.0 + linear * gc_gc0 + quadratic * gc_gc02,
                           crdna.constants.MIN_CURVE, crdna.constants.MAX_CURVE)
# gc_bias

#............................................................................
def get_haploid_poisson_expectations(
        scale, linear, quadratic, 
        raw_profiles, bin_parameters):
    """For each chromosome in each cell, compute the ploidy=1 expected number
    of reads per bin accounting for GC bias and scaling."""
    epsilon = 1e-6
    n_chrom = len(raw_profiles)
    n_cells = raw_profiles[0].shape[0]
    haploid_poisson_expectations = [None] * n_chrom
    for chrom_index in xrange(n_chrom):
        gc = bin_parameters[chrom_index][3, :]
        left_selector = np.where(gc < crdna.constants.MIN_GC)[0]
        right_selector = np.where(gc > crdna.constants.MAX_GC)[0]
        gc_gc0 = gc - crdna.constants.GC_ORIGIN
        gc_gc02 = gc_gc0 * gc_gc0
        left_gc_gc0 = crdna.constants.MIN_GC - crdna.constants.GC_ORIGIN
        left_gc_gc02 = left_gc_gc0 * left_gc_gc0
        right_gc_gc0 = crdna.constants.MAX_GC - crdna.constants.GC_ORIGIN
        right_gc_gc02 = right_gc_gc0 * right_gc_gc0
        n_bins = raw_profiles[chrom_index].shape[1]
        tmp = np.zeros((n_cells,n_bins))
        for cell in xrange(n_cells):
            trend = gc_bias(scale[cell], linear[cell], quadratic[cell], gc_gc0, gc_gc02)
            left_trend = gc_bias(scale[cell], linear[cell], quadratic[cell], left_gc_gc0, left_gc_gc02)
            right_trend = gc_bias(scale[cell], linear[cell], quadratic[cell], right_gc_gc0, right_gc_gc02)
            trend[left_selector] = left_trend
            trend[right_selector] = right_trend
            trend[trend <= 0.0] = epsilon
            tmp[cell, :] = trend
        # for cell
        haploid_poisson_expectations[chrom_index] = tmp
    # for chrom_index
    return haploid_poisson_expectations
# get_haploid_poisson_expectations

#............................................................................
def process_cnv_calls(
        raw_profiles, mask, bin_parameters, 
        reference_path, sex,
        scale, linear, quadratic, cnv_calls, bin_size,
        logp=None):
    """For the provided raw profiles and cnv_calls calculate CNV confidence
    scores for each CNV call.

    An intermediate result of this calculation is the logp matrix of per-bin
    confidence scores.

    The returned result is a tuple of logp, cnv_calls 
    :returns: logp is an array (size ncontigs) of nparrays (ncells x nbins)
    :returns: cnv_calls the original cnv calls with confidence supplied
    """
    poisson_expectations = get_haploid_poisson_expectations(
        np.array(scale), linear, quadratic, 
        raw_profiles, bin_parameters)

    if logp is None:
        logp = calculate_logposterior_matrix(raw_profiles, poisson_expectations, mask, 
            cnv_calls, reference_path, bin_size)
    #
    # calculate CNV scores
    #
    cnv_calls = estimate_cnv_confidence_score_v2(raw_profiles, cnv_calls,
        reference_path, logp, bin_size )

    return logp, cnv_calls
# process_cnv_calls

#............................................................................
