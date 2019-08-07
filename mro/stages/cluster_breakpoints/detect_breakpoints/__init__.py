#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Computes the log likelihood ratio that tests for whether a given position
# in a cell profile is a breakpoint or not. Computed for all cells and bins.
#

import pandas as pd
import numpy as np
import json
import sys

import martian
from longranger.cnv import contig_manager, coverage_matrix
from crdna.utils import load_genome_data, load_h5
from crdna.breakpoints import parabola, call_cnvs, compute_corr_dist_to_all, \
    get_segment_bdy_from_index, compute_segment_data, get_ploidy_vector, \
    get_segment_window_size
import crdna.constants
from plot_utils import aggregate_matrix

__MRO__ = """
stage DETECT_BREAKPOINTS(
    in  string   reference_path,
    in  h5       tracks,
    in  h5       profiles,
    in  map      params,
    in  json     gc_norm_params,
    in  bool     is_singlecell,
    out h5       denoised_profiles,
    src py       "stages/cluster_breakpoints/detect_breakpoints",
) split using (
    in  map      chunk,
)
"""

def split(args):
    MAX_CHUNKS = 30
    MIN_CELLS_PER_CHUNK = 100

    ## TODO : store ncells in the profiles.h5 as a constant so we don't have
    ## to do this to get the number of cells
    ref = contig_manager.contig_manager(args.reference_path)
    chrom = ref.primary_contigs(allow_sex_chromosomes=True)[0]
    
    store = pd.HDFStore(args.profiles, "r")
    ncells, _ = store["/contigs/"+chrom].shape
    store.close( )

    ## no cells, do nothing!
    if ncells < 1:
        return {'chunks': [], 'join': {}}
    
    nchunks = np.clip(ncells/MIN_CELLS_PER_CHUNK +
        int(ncells % MIN_CELLS_PER_CHUNK != 0), 1, MAX_CHUNKS)
    cells_per_chunk = ncells/nchunks + int(ncells % nchunks != 0)
    
    mat_size_gb = coverage_matrix.get_genome_matrix_size_gb(args.profiles)
    chunk_mem_gb = int(np.ceil(4*mat_size_gb/ncells*cells_per_chunk + 1))
    join_mem_gb =  int(np.ceil(4*mat_size_gb + 1))

    ## if this is a multi species sample do nothing
    if len(ref.list_species( )) > 1:
        return {'chunks': [{'chunk': {'start': 0,
                                      'end': ncells,
                                      'ncells': ncells}}],
                'join': {'__mem_gb': join_mem_gb}}

    chunk_defs = [{'chunk': {'start': i, 'end': min(i+cells_per_chunk, ncells),
        'ncells': ncells}, '__mem_gb': chunk_mem_gb} 
        for i in xrange(0, ncells, cells_per_chunk)]
    
    return {'chunks': chunk_defs, 'join': {'__mem_gb': join_mem_gb}}

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()
    
    start_cell = args.chunk["start"]
    end_cell   = args.chunk["end"]
    
    ref = contig_manager.contig_manager( args.reference_path )
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)

    ## load genome data and cell profiles
    X, gmask, bdy = load_genome_data( args.profiles, args.tracks, chroms,
        start_cell=start_cell, end_cell=end_cell, integer=False, rounding=False, 
        mappability_threshold = crdna.constants.MAPPABILITY_THRESHOLD )
    
    # compute chromosome boundaries after masking by gmask
    cbdy = np.zeros_like(bdy)
    for i in xrange(1, len(bdy)):
        cbdy[i] = (gmask[0:bdy[i]].sum())

    ## load GC info and create GC track
    gctrack = []
    store = pd.HDFStore(args.tracks, "r")
    for chrom in chroms:
        gctrack.extend(store["/GC/"+chrom].values)
    gctrack = np.array(gctrack)[gmask]
    store.close( )
    
    ## if calling on nodes use scale factors from cells
    if args.is_singlecell:
        scale_guess_chunk = [None for _ in xrange(start_cell, end_cell)]
        cell_offset = 0
    else:
        scale_guess = load_h5(args.profiles, "scale_guess")
        scale_guess_chunk = [[s] for s in scale_guess[start_cell:end_cell]]
        ## num cells = num internal nodes + 1
        cell_offset = args.chunk["ncells"]+1
        


    ncells = X.shape[0]
    nbins  = gmask.sum( )
    P = 2*np.ones((ncells, nbins), dtype="int8")
    S = np.zeros((ncells, nbins), dtype=bool)
    sdfs = []
    scale_factors = np.zeros(ncells)
    pconf = np.zeros(ncells)
    windows = np.zeros(ncells, dtype=int)

    gc_norm_params = json.load(open(args.gc_norm_params, "r"))

    ## initialize parameters
    read_threshold = crdna.constants.BREAKPOINT_READ_THRESHOLD
    heuristics = crdna.constants.BREAKPOINT_CALLER_HEURISTICS
    ## override/augment heuristics by supplied params
    if args.params is not None:
        for k, v in args.params.iteritems():
            heuristics[k] = v

    ## log heuristics used
    martian.log_info("Heuristics used:")
    for k, v in heuristics.iteritems():
        martian.log_info("%s: %s"%(str(k), str(v)))
    
    debug_out = open(outs.debug, "w")

    if len(ref.list_species( )) == 1:
        for i in xrange(ncells):
            debug_out.write("-"*80+"\n")
            debug_out.write("Cell %d\n"%(cell_offset+start_cell+i))

            ## GC coefficients
            gc_linear = gc_norm_params["linear"][start_cell+i]
            gc_quadratic = gc_norm_params["quadratic"][start_cell+i]

            ## GC correction track for cell
            xi = parabola(gctrack, crdna.constants.GC_ORIGIN, 
                gc_linear, gc_quadratic)
            xi_low = parabola(crdna.constants.MIN_GC, 
                crdna.constants.GC_ORIGIN, gc_linear, gc_quadratic)
            xi_high = parabola(crdna.constants.MAX_GC, 
                crdna.constants.GC_ORIGIN, gc_linear, gc_quadratic)
            xi[gctrack < crdna.constants.MIN_GC] = xi_low
            xi[gctrack > crdna.constants.MAX_GC] = xi_high
           
            y = X[i][gmask]
            
            ## do the CNV calling
            ploidy, S[i], gap, sdf, sf = call_cnvs(y, xi, ref, cbdy,
                scale_guess=scale_guess_chunk[i], log_func=debug_out.write,
                **heuristics)

            scale_factors[i] = sf
            sdfs.append(sdf)
            P[i] = np.clip(ploidy, 0, np.iinfo("int8").max-1)
            pconf[i] = gap
            windows[i] = get_segment_window_size(y, read_threshold)
            debug_out.flush()
    debug_out.close()

    out = pd.HDFStore( outs.denoised_profiles, "w" )
    out["/quantized"] = pd.DataFrame(P)
    out["/segment_index"] = pd.DataFrame(S)
    out["/scaling_data"] = pd.Series(sdfs)
    out["/scale_factor"] = pd.Series(scale_factors)
    out["/ploidy_conf"] = pd.Series(pconf)
    out["/windows"] = pd.Series(np.clip(windows, 1, None))
    out.close( )

## join chunks
def join(args, outs, chunk_defs, chunk_outs):
    ref = contig_manager.contig_manager( args.reference_path )
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)

    ## load genome data and cell profiles
    X, gmask, bdy = load_genome_data( args.profiles, args.tracks, chroms,
        integer=False, rounding=False, 
        mappability_threshold = crdna.constants.MAPPABILITY_THRESHOLD )
    nbins = gmask.sum()
    
    # compute chromosome boundaries after masking by gmask
    cbdy = np.zeros_like(bdy)
    for i in xrange(1, len(bdy)):
        cbdy[i] = (gmask[0:bdy[i]].sum())

    ## load GC info and create GC track
    gctrack = []
    store = pd.HDFStore(args.tracks, "r")
    for chrom in chroms:
        gctrack.extend(store["/GC/"+chrom].values)
    gctrack = np.array(gctrack)[gmask]
    store.close( )
    gc_norm_params = json.load(open(args.gc_norm_params, "r"))

    ## Aggregate data structures from individual chunks
    P = np.zeros((0, nbins), dtype="int8") # ploidy per cell
    S = np.zeros((0, nbins), dtype=bool)   # segment index per cell
    sdfs = []                              # scaling dataframes per cell
    scale_factors = []                     # scale factor per cell
    pconf = np.zeros((0,), dtype=float)    # scaling confidence per cell
    windows = np.zeros((0,), dtype=int)    # segment window size per cell

    logger = sys.stdout.write 

    ## add logging info from chunks
    for chunk_out in chunk_outs:
        with open(chunk_out.debug, "r") as debug_in:
            for line in debug_in:
                logger(line)
    logger("\n"+"*"*80+"\n")

    for chunk_out, chunk_def in zip(chunk_outs, chunk_defs):
        start_cell = chunk_def.chunk["start"]
        end_cell   = chunk_def.chunk["end"]
        chunk_store = pd.HDFStore( chunk_out.denoised_profiles, "r" )
        p_chunk = chunk_store["/quantized"].values
        s_chunk = chunk_store["/segment_index"].values
        sf_chunk = chunk_store["/scale_factor"].values
        pc_chunk = chunk_store["/ploidy_conf"].values
        sd_chunk = list(chunk_store["/scaling_data"])
        w_chunk = chunk_store["/windows"].values
        chunk_store.close( )
       
        if P.shape[0] == 0:
            ncells = chunk_def.chunk["ncells"]
            nbins  = p_chunk.shape[1]
            P = np.zeros((ncells, nbins), dtype="int8")
            S = np.zeros((ncells, nbins), dtype=bool)
            scale_factors = np.zeros(ncells, dtype=float)
            pconf = np.zeros(ncells, dtype=float)
            windows = np.zeros(ncells, dtype=int)
        P[start_cell:end_cell, :] = p_chunk
        S[start_cell:end_cell, :] = s_chunk
        scale_factors[start_cell:end_cell] = sf_chunk
        sdfs.extend(sd_chunk)
        pconf[start_cell:end_cell] = pc_chunk
        windows[start_cell:end_cell] = w_chunk
    
    ## Find cells with low scaling confidence

    fix_scaling = np.zeros(0, dtype=int)
    if args.is_singlecell:
        cell_offset = 0
        fix_scaling = np.where((pconf >= 0) & (pconf <= 0.02))[0]
    else:
        cell_offset = X.shape[0]+1

    good_cells = np.where(np.logical_or(pconf == -2, pconf > 0.10))[0]
    
    agg_window = int(np.median(windows if len(windows) > 0 else [sum(gmask)]))
    X_agg = aggregate_matrix(X[:, gmask], agg_window)
    
    ## initialize parameters
    heuristics = crdna.constants.BREAKPOINT_CALLER_HEURISTICS
    ## override/augment heuristics by supplied params
    if args.params is not None:
        for k, v in args.params.iteritems():
            heuristics[k] = v
    logger("%d cells with low ploidy confidence\n"%len(fix_scaling))

    for cell in fix_scaling:
        if len(good_cells) == 0:
            continue
        logger("-"*80+"\n")
        logger("Fixing cell %d\n"%(cell+cell_offset))

        ## GC coefficients
        gc_linear = gc_norm_params["linear"][cell]
        gc_quadratic = gc_norm_params["quadratic"][cell]

        ## GC correction track for cell
        xi = parabola(gctrack, crdna.constants.GC_ORIGIN, 
            gc_linear, gc_quadratic)
        xi_low = parabola(crdna.constants.MIN_GC, 
            crdna.constants.GC_ORIGIN, gc_linear, gc_quadratic)
        xi_high = parabola(crdna.constants.MAX_GC, 
            crdna.constants.GC_ORIGIN, gc_linear, gc_quadratic)
        xi[gctrack < crdna.constants.MIN_GC] = xi_low
        xi[gctrack > crdna.constants.MAX_GC] = xi_high
 
        y = X[cell][gmask]

        ## find the correlation distance to all cells that were scaled
        ## confidently. Then take all matches with > 90% correlation and 
        ## compute the median ploidy over these cells. Find the closest
        ## scaling solution to the median and declare that the answer.
        all_corrs = compute_corr_dist_to_all(X_agg[cell][np.newaxis, :], X_agg)
        good_corrs = all_corrs[good_cells]
        best_matches = good_cells[good_corrs > 0.90]
        if len(best_matches) == 0:
            continue
        best_guess_ploidy = np.median(P[best_matches, :].mean(axis=1))

        best_scaling_soln = np.argmin(np.abs(sdfs[cell]["aploidy"].values 
            - best_guess_ploidy))
        lam_best = sdfs[cell].loc[best_scaling_soln]["lam"]    
   
        sindex = S[cell]
        segment_bdy2 = get_segment_bdy_from_index(sindex)
        window = get_segment_window_size(y, heuristics["ll_read_threshold"]) 

        segment_means2, _ = compute_segment_data(y, xi, segment_bdy2, window) 
        ploidy = get_ploidy_vector(y, segment_means2, segment_bdy2, lam_best)
        delta_ploidy = np.abs(P[cell].mean() - ploidy.mean())
        logger("Ploidy: %.2f -> %.2f\n"%(P[cell].mean(),
            ploidy.mean()))

        P[cell] = np.clip(ploidy, 0, np.iinfo("int8").max-1).astype("int8")
        if delta_ploidy > 0.1:
            pconf[cell] = -4

    ## Compute read depth
    depth = np.zeros_like(pconf)
    for cell in xrange(len(pconf)):
        depth[cell] = X[cell][gmask].mean()

    ## Write data to h5
    store = pd.HDFStore( outs.denoised_profiles, "w" )
    store["/quantized"] = pd.DataFrame(P)
    store["/scale_factor"] = pd.Series(scale_factors)
    store["/reads_per_bin"] = pd.Series(depth) 
    store["/segment_index"] = pd.DataFrame(S)
    store["/ploidy_conf"] = pd.Series(pconf)
    store["/scaling_data"] = pd.Series(sdfs)
    store["/windows"] = pd.Series(windows)
    segment_windows = int(np.median(windows)) if len(windows) else 1
    constants = load_h5(args.profiles, "constants").to_dict()
    constants["segment_windows"] = segment_windows
    store["constants"] = pd.Series(constants)
    store.close( )


