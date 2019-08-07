#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# aggregate single cell and internal node data
import martian
import pandas as pd
import json
import os
import sys
import subprocess
import numpy as np

from crdna.clustering import compute_heterogeneity
from crdna.constants import MAPPABILITY_THRESHOLD
from crdna.utils import load_h5
from longranger.cnv import contig_manager

__MRO__ = """
stage AGGREGATE_NODES(
    in  string reference_path,
    in  h5     tracks,
    in  h5     internal_cnv_tracks,
    in  h5     sc_cnv_tracks,
    in  bed    internal_cnv_calls,
    in  bed    sc_cnv_calls,
    in  bed    internal_unmerged_cnv_calls,
    in  bed    sc_unmerged_cnv_calls,
    in  json   sc_gc_params,
    in  json   internal_gc_params,
    in  h5     sc_norm_profiles,
    in  h5     internal_norm_profiles,
    in  h5     tree_data,
    out json   node_gc_params,
    out h5     norm_node_profiles,
    out h5     node_cnv_tracks,
    out bed    node_cnv_calls,
    out bed    node_unmerged_cnv_calls,
    out h5     tree_data,
    src py     "stages/cluster_breakpoints/aggregate_nodes",
)
"""

def split(args):
    ref = contig_manager.contig_manager( args.reference_path )
    constants = load_h5(args.sc_norm_profiles, "constants").to_dict()
    ncells = constants["ncells"]
    window_size = constants["window_size"]
    # maximum memory usage is the maximum of these four values:
    # sumbins = sum(len(c) for c in primary_contigs)/window_size
    # maxbins = max(len(c) for c in all_contigs)/window_size
    # X + Q + H = ((2*sizeof(i8) + 2*sizeof(f32)) * ncells * sumbins)
    # occupancy = sizeof(f32) * levels(=6) * (ncells - 1) * sumbins / nchunks(=100)
    # het = X + Q + H + occupancy
    # X + Y + Z = ((2*sizeof(float)) * ncells * maxbins)
    # merged_bed = sc_cnv_calls_bed + internal_cnv_calls_bed
    # unmerged_bed = sc_unmerged_cnv_calls_bed + internal_unmerged_cnv_calls_bed
    # * NOTE: ask for the double the matrix sizes to acct for intermediate values
    f32sz = 4
    sumbins = sum(ref.contig_lengths[c]/window_size+1 for c in ref.primary_contigs())
    maxbins = max(ref.contig_lengths[c]/window_size+1 for c in ref.list_all_contigs())
    XQH_mem_gb = float((2 + 2*f32sz) * ncells * sumbins)/1e9
    occ_mem_gb = float(f32sz * 6 * (ncells - 1) * sumbins/100)/1e9
    het_mem_gb = XQH_mem_gb + occ_mem_gb
    XYZ_mem_gb = 2 * float(f32sz * ncells * maxbins) / 1e9
    merged_bed_gb = os.path.getsize(args.sc_cnv_calls)/1e9 + \
                    os.path.getsize(args.internal_cnv_calls)/1e9 + 1
    unmerged_bed_gb = os.path.getsize(args.sc_unmerged_cnv_calls)/1e9 + \
                      os.path.getsize(args.internal_unmerged_cnv_calls)/1e9 + 1
    mem_gb = int(np.ceil(max(het_mem_gb, XYZ_mem_gb, merged_bed_gb, unmerged_bed_gb))) + 3
    return {'chunks': [], 'join': {'__mem_gb': mem_gb}}

def join(args, outs, chunk_defs, chunk_outs):
    ## merge gc params jsons
    node_gc_params = {}
    sc_gc_params = json.load(open(args.sc_gc_params, "r"))
    internal_gc_params = json.load(open(args.internal_gc_params, "r"))

    ncells = len(sc_gc_params['linear'])
    nnodes = 2*ncells - 1

    for key in ["scale", "linear", "quadratic"]:
        node_gc_params[key] = sc_gc_params[key] + internal_gc_params[key]
    with open(outs.node_gc_params, "w") as out:
        json.dump(node_gc_params, out, indent=4)

    ref = contig_manager.contig_manager(args.reference_path)
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)
    index_chrom = dict([(str(i), c) for i, c in enumerate(chroms)])
    chrom_index = dict([(c, str(i)) for i, c in enumerate(chroms)])
    tmp = martian.make_path('tmp.bed')
    tmp_dir = os.path.dirname(tmp)
    tmp_sorted = martian.make_path('tmp_sorted.bed')
    calls = [[args.sc_cnv_calls, args.internal_cnv_calls],
             [args.sc_unmerged_cnv_calls, args.internal_unmerged_cnv_calls]]
    out_calls = [outs.node_cnv_calls, outs.node_unmerged_cnv_calls]
    for calls, out in zip(calls, out_calls):
        with open(tmp, 'w') as outf:
            for f in calls:
                for l in open(f):
                    fields = l.split()
                    # offset internal node indices by ncells
                    if f == calls[1]:
                        fields[3] = str(int(fields[3]) + ncells)
                    # fix type of confidence field to integer
                    fields[-1] = str(int(float(fields[-1])))
                    # replace index number at start for sorting
                    fields[0] = chrom_index[fields[0]]
                    outf.write('\t'.join(fields) + '\n')

        no_unicode = dict(LC_ALL='C')
        tmp_mem_gib = max(1, int(np.ceil(float(os.path.getsize(tmp)) / (1024**3))))
        try:
            subprocess.check_call(['sort', '-k1,1n', '-k2,2n', '-k3,3n',
                                   '--parallel=1',  # force sort to use 1 thread
                                   '-S', '{}G'.format(tmp_mem_gib),
                                   '-T', tmp_dir,
                                   '-o', tmp_sorted, tmp],
                                  env=no_unicode, stderr=sys.stderr)
        # on some systems, --parallel is unavailable
        except subprocess.CalledProcessError:
            subprocess.check_call(['sort', '-k1,1n', '-k2,2n', '-k3,3n',
                                   # will by default only use 1 thread
                                   '-S', '{}G'.format(tmp_mem_gib),
                                   '-T', tmp_dir,
                                   '-o', tmp_sorted, tmp],
                                  env=no_unicode, stderr=sys.stderr)

        # strip index column into outfile
        with open(out, 'w') as outf:
            version = martian.get_pipelines_version()
            outf.write("#cellranger-dna {}\n".format(version))
            outf.write("#reference genome: {}\n".format(args.reference_path))
            outf.write("#chrom\tstart\tend\tid\tcopy_number\tevent_confidence\n")
            for l in open(tmp_sorted):
                l = l.split('\t')
                l[0] = index_chrom[l[0]]
                outf.write('\t'.join(l))

    os.remove(tmp)
    os.remove(tmp_sorted)

    ## cnv tracks file
    sc_windows = load_h5(args.sc_cnv_tracks, "windows")
    internal_windows = load_h5(args.internal_cnv_tracks, "windows")
    windows = sc_windows.append(internal_windows).values
    constants = load_h5(args.sc_cnv_tracks, "constants")
    
    sc_ploidy_conf = scale_confidence_score(load_h5(args.sc_cnv_tracks, 
        "ploidy_conf").values)
    internal_ploidy_conf = scale_confidence_score(load_h5(
        args.internal_cnv_tracks, "ploidy_conf").values)
    
    sc_scale_factor= load_h5(args.sc_cnv_tracks, "scale_factor")
    internal_scale_factor = load_h5(args.internal_cnv_tracks, "scale_factor")

    sc_rpb= load_h5(args.sc_cnv_tracks, "reads_per_bin")
    internal_rpb= load_h5(args.internal_cnv_tracks, "reads_per_bin")
    
    X = load_h5(args.sc_cnv_tracks, "cnv_tracks").values
    nbins = X.shape[1]
    Q = np.zeros((nnodes, nbins), dtype=X.dtype)
    Q[0:ncells, :] = X
    del X
    Q[ncells:, :] = load_h5(args.internal_cnv_tracks, "cnv_tracks").values

    store = pd.HDFStore(outs.node_cnv_tracks, "w")
    store["constants"] = constants
    store["windows"] = sc_windows.append(internal_windows)
    store["ploidy_conf"] = sc_ploidy_conf.append(internal_ploidy_conf)
    store["scale_factor"] = sc_scale_factor.append(internal_scale_factor)
    store["reads_per_bin"] = sc_rpb.append(internal_rpb)
    store["cnv_tracks"] = pd.DataFrame(Q)
    store.close()
    
    ## Compute heterogeneity and store in tree_data
    ref = contig_manager.contig_manager(args.reference_path)
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)
    if args.tracks is None:
        gmask = np.ones(nbins, dtype=bool)
    else:
        gmask = []
        maptrack = pd.HDFStore(args.tracks, "r")
        for chrom in chroms:
            gmask.extend(maptrack["/map/"+chrom].values > MAPPABILITY_THRESHOLD)
        maptrack.close( )
        gmask = np.array(gmask)

    ## update tree data
    # load tree
    store = pd.HDFStore( args.tree_data, "r" )
    Z = store["/Z"].values
    distances = store["/distances"].values
    constants = store["/constants"]
    store.close( )

    # Compute the heterogeneity at every *internal* node of the tree
    # obviously the heterogeneity is zero at every leaf, so don't
    # store a bunch of zeros
    levels = 6
    het = compute_heterogeneity(Q, Z, gmask, windows, levels=levels)

    del Q

    # dump to disk
    store = pd.HDFStore( outs.tree_data, "w" )
    store["Z"] = pd.DataFrame(Z)
    store["het"] = pd.DataFrame(het)
    store["distances"] = pd.Series(distances)
    store["windows"] = pd.Series(windows)
    store["constants"] = constants
    store.close( )

    del het

    ## normalized profiles
    sc_store = pd.HDFStore(args.sc_norm_profiles, "r")
    internal_store = pd.HDFStore(args.internal_norm_profiles, "r")
    out_store = pd.HDFStore(outs.norm_node_profiles, "w")
    out_store["/constants"] = sc_store["/constants"]
    for chrom in chroms:
        ## first do the /contigs
        X = sc_store["/contigs/"+chrom].values
        Y = internal_store["/contigs/"+chrom].values
        assert X.shape[1] == Y.shape[1]
        nbins = X.shape[1]
        Z = np.zeros((2*ncells-1, nbins), dtype=X.dtype)
        Z[:ncells, :] = X
        Z[ncells:, :] = Y
        out_store["/contigs/"+chrom] = pd.DataFrame(Z)
        del X, Y, Z

        ## next do the /masks
        out_store["/masks/"+chrom] = sc_store["/masks/"+chrom]
    ## gc params
    for key in ["scale", "linear", "quadratic"]:
        out_store["/gc_params/"+key] = pd.concat([sc_store["/gc_params/"+key],
            internal_store["/gc_params/"+key]], ignore_index=True)

    ## do the normalization metrics
    out_store["/normalization_metrics"] =sc_store["normalization_metrics"].append(internal_store["/normalization_metrics"], ignore_index=True)

    out_store.close()
    sc_store.close()
    internal_store.close()

def scale_confidence_score( pconf ):
    """ Scale the confidence factor so it can be stored as int8. """
    pos_score = pconf >= 0
    pconf_score = np.zeros_like(pconf, dtype="int8")
    pconf_score[pos_score] = np.clip(np.round(100*pconf[pos_score]), None,
        np.iinfo("int8").max).astype("int8")
    pconf_score[~pos_score] = pconf[~pos_score].astype("int8")
    return pd.Series(pconf_score, dtype="int8")
