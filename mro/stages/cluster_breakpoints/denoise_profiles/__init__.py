#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import pandas as pd
import numpy as np
from crdna.utils import load_genome_data
import crdna.constants
from longranger.cnv import contig_manager
from collections import Counter
from crdna.breakpoints import get_breakpoint_positions, \
    break_segments_at_points, validate_segment_intervals, find_best_scale_v14,\
    get_ploidy_vector, parabola
import sys
import json

__MRO__ = '''
stage DENOISE_PROFILES(
    in  string   reference_path,
    in  h5       tracks,
    in  h5       profiles,
    in  h5       ll_ratios,
    in  map      params,
    in  json     gc_norm_params,
    out h5       denoised_profiles,
    src py       "stages/cluster_breakpoints/denoise_profiles",
) split using
(
    in  int      dummy,
)
'''

## Trivial chunking implemented to specify __mem_gb for chunk
def split( args ):
    return {'chunks': [], 'join': {'__mem_gb': 64}}

## Trivial join
def join( args, outs, chunk_defs, chunk_outs ):
    args.coerce_strings()
    outs.coerce_strings()

    ref = contig_manager.contig_manager(args.reference_path)
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)

    ## load genome data and cell profiles
    X, gmask, bdy = load_genome_data( args.profiles, args.tracks, chroms,
        mappability_threshold=crdna.constants.MAPPABILITY_THRESHOLD )

    # compute chromosome boundaries after masking by gmask
    cbdy = np.zeros_like(bdy)
    for i in xrange(1, len(bdy)):
        cbdy[i] = (gmask[0:bdy[i]].sum())

    ## load GC info and create GC emission track
    gctrack = []
    store = pd.HDFStore(args.tracks, "r")
    for chrom in chroms:
        gctrack.extend(store["/GC/"+chrom].values)
    gctrack = np.array(gctrack)[gmask]
    store.close( )
       
    store = pd.HDFStore( args.ll_ratios, "r" )
    llrs = store["/llrs"].values
    store.close( )

    nbins = gmask.sum( )
    ncells = X.shape[0]

    ## Heuristics to define breakpoints
    ll_threshold    = 5
    delta_threshold = 0.10

    Y_quant       = np.zeros((ncells, nbins), dtype="int8")
    scale_factor  = np.zeros(ncells)
    windows_per_cell = []

    gc_norm_params = json.load(open(args.gc_norm_params, "r"))
    print "Starting loop over cells"
    sys.stdout.flush( )
    for i in xrange(ncells):
        print "-"*80
        print "Cell", i
        sys.stdout.flush( )

        ## genome profile
        y  = X[i][gmask]
        ## log likelihood ratio profile
        ll = llrs[i]
        
        ## GC coefficients
        gc_linear = gc_norm_params["linear"][i]
        gc_quadratic = gc_norm_params["quadratic"][i]

        ## GC correction track for cell
        xi = parabola(gctrack, crdna.constants.GC_ORIGIN, 
            gc_linear, gc_quadratic)
        xi_low = parabola(crdna.constants.MIN_GC, 
            crdna.constants.GC_ORIGIN, gc_linear, gc_quadratic)
        xi_high = parabola(crdna.constants.MAX_GC, 
            crdna.constants.GC_ORIGIN, gc_linear, gc_quadratic)
        xi[gctrack < crdna.constants.MIN_GC] = xi_low
        xi[gctrack > crdna.constants.MAX_GC] = xi_high

        ## Define breakpoints
        ##
        
        bp_cands2 = get_breakpoint_positions(y, ll, xi, 
            ll_threshold=ll_threshold, delta_threshold=delta_threshold)

        assert bp_cands2[0] == 0, "genome start must be breakpoint"
        assert bp_cands2[-1] == y.shape[0], "genome end must be breakpoint"

        ## define segments using breakpoints
        segment_bdy = []
        for j in xrange(len(bp_cands2)-1):
            segment_bdy.append((bp_cands2[j], bp_cands2[j+1]))

        ## add chromosome boundaries as mandatory breakpoints
        segment_bdy = break_segments_at_points(segment_bdy, cbdy, verbose=False)
        validate_segment_intervals(segment_bdy, cbdy)

        ## aggregate bins within a segment to resolution given by window
        ## and compute segment mean read counts and lengths
        window = int(np.round(crdna.constants.BREAKPOINT_READ_THRESHOLD/
            np.median(y[y > 0])))
        window = np.clip(window, 1, None)
        windows_per_cell.append(window)

        segment_means = []
        segment_lengths = []
        for s, e in segment_bdy:
            segment = y[s:e]
            xi_piece = xi[s:e]
            length = e - s
            agg = []
            xi_agg = []
            j = 0
            while j < length:
                piece = segment[j:j+window]
                assert len(piece) > 0, "%d, %d-%d"%(j, s, e)
                corr  = float(window)/len(piece)
                agg.append(corr*piece.sum())
                xi_agg.append( xi_piece[j:j+window].mean( ) )
                j += window
            agg = np.array(agg)
            xi_agg = np.array(xi_agg)
            ## remove outliers
            med = np.median(agg)
            mad = np.abs( agg - med)
            mmad = np.median(mad)
            outlier_mask = mad <= 5*mmad
            segment_means.append( 
                np.sum(agg[outlier_mask])/np.sum(xi_agg[outlier_mask]) )
            segment_lengths.append( e-s )
        segment_means = np.array(segment_means)
        segment_lengths = np.array(segment_lengths)
        
        ## Find the scaling factor to produce integer ploidies

        ## Heuristics
        
        # max ploidy to assign to initially chosen long segment
        max_ploidy_long=10
        
        # max value of segment mean to consider "zero ploidy"
        zero_ploidy_count=crdna.constants.BREAKPOINT_READ_THRESHOLD/4.0
        
        # longest segment with segment mean > zero_ploidy_count
        # that we will push to zero ploidy
        max_segment_push_to_zero=200

        # prior params
        prior_params = {"prior_mean"   : args.params.get("prior_mean", 2.0),
                        "prior_std"    : args.params.get("prior_std", 1.0)}
        min_ploidy = args.params.get("min_ploidy", None)
        max_ploidy = args.params.get("max_ploidy", None)
        lam_best = find_best_scale_v14( y, segment_bdy, segment_means, 
            segment_lengths, window, max_ploidy_long=max_ploidy_long,
            zero_ploidy_count=zero_ploidy_count, prior_params=prior_params,
            max_segment_push_to_zero=max_segment_push_to_zero,
            min_ploidy=min_ploidy, max_ploidy=max_ploidy, verbose=True )
            
        print "Scaling factor:"
        print lam_best

        assert lam_best > 0.001, "run away to zero"
        scale_factor[i] = lam_best/window

        ## Compute the ploidy vector

        ploidy = get_ploidy_vector(y, segment_means, segment_bdy, lam_best)
        
        ## Set max ploidy at 127
        ploidy = np.clip(ploidy, None, np.iinfo("int8").max)

        print "Ploidies encountered"
        print Counter(ploidy).most_common( )
            
        Y_quant[i, :] = ploidy.astype("int8")
    
    ## store in output h5
    out_store = pd.HDFStore( outs.denoised_profiles, "w")
    out_store["/constants"] = pd.Series({"segment_windows" : 
        int(np.median(windows_per_cell))})
    out_store["/scale_factor"] = pd.Series(scale_factor)
    out_store["/quantized"] = pd.DataFrame(Y_quant)
    out_store.close( )

