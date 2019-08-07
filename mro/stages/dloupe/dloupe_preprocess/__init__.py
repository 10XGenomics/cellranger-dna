#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

from longranger.cnv import contig_manager
import martian
import subprocess
import math
import json
import os
import tenkit.log_subprocess as tk_subproc
import tenkit.reference as tk_ref

__MRO__ = """
stage DLOUPE_PREPROCESS (
    in  string sample_id,
    in  string sample_desc,
    in  string reference_path,
    in  int    profile_bin_size,
    in  h5     normalized_node_profiles,
    in  bed    node_cnv_calls,
    in  h5     tree_data,
    in  h5     tracks,
    in  csv    per_cell_summary_metrics,
    out dloupe output_for_dloupe,
    out json   contig_info_json,
    src py     "stages/dloupe/dloupe_preprocess",
) split using (
)
"""

def get_contig_info(args):
    manager = contig_manager.contig_manager(args.reference_path)
    contig_info = {
        "contig_order": {},
        "contig_lengths": {}
    }
    contig_lengths = manager.get_contig_lengths()
    for idx, contig in enumerate(manager.contigs["primary_contigs"]):
        contig_info["contig_order"][contig] = idx
        contig_info["contig_lengths"][contig] = contig_lengths[contig]
    contig_info["species"] = manager.list_species()
    return contig_info


def split(args):
    contig_info = get_contig_info(args)
    length_bases = [v for v in contig_info["contig_lengths"].values()]
    total_bins = 0
    for bases in length_bases:
        bins, over = divmod(bases, args.profile_bin_size)
        total_bins += bins
        if over > 0:
            total_bins += 1

    with open(args.per_cell_summary_metrics) as f:
        summary_line_count = len(f.readlines())
    total_nodes = 2*(summary_line_count-1)-1
    martian.log_info("Bins: %d, Nodes: %d" % (total_bins, total_nodes))

    # observed worst-case dlconverter memory scenario is when the pipeline is
    # copying two elements: the complete float64 float read depth
    # 2*(8 x (nodes*2(cells)-1) and the int8 ploidy (in theory,
    # 2*nodes*2(cells)-1, for a total of 18*(total_nodes*total_bins),  but
    # observed is closer to 21x.  Going higher (24x) to be safe.
    #
    # Estimated hit for a 1000-cell dataset is 7GB.
    # Update 5/29/18: Test dataset JD-100_79 broke these assumptions, and requested
    # slightly more memory than expected. The easiest solution is to just add a bit
    # of a buffer. Changing it to 26x and adding a 4GB minimum. This is an ad-hoc
    # solution, but it *should* be safe generally. Better to err on the side of
    # not crashing things.
    mem_bytes = total_nodes*total_bins*26
    mem_gb = int(math.ceil(mem_bytes / (1024.0*1024.0*1024.0))) + 6

    martian.log_info("Bins: %d, Nodes: %d, Bytes: %d" % (total_bins, total_nodes, mem_bytes))

    return {'chunks': [], 'join': {'__mem_gb': mem_gb, '__threads': 2}}

def join(args, outs, chunk_defs, chunk_outs):
    contig_info = get_contig_info(args)
    with open(outs.contig_info_json, 'w') as outfile:
        json.dump(contig_info, outfile)

    call = ["dlconverter",
            args.sample_id,
            "--output", outs.output_for_dloupe,
            "--description", args.sample_desc,
            "--node-profile-h5", args.normalized_node_profiles,
            "--contig-info-json", outs.contig_info_json,
            "--merged-bed", args.node_cnv_calls,
            "--tree-data", args.tree_data,
            "--tracks", args.tracks,
            "--per-cell-summary", args.per_cell_summary_metrics]

    gene_annotation_path = tk_ref.get_loupe_genes(args.reference_path)
    if os.path.exists(gene_annotation_path):
        call.extend(["--gene-annotations", gene_annotation_path])

    # the sample desc may be unicode, so send the whole
    # set of args str utf-8 to check_output
    unicode_call = [arg.encode('utf-8') for arg in call]

    martian.log_info("Running dlconverter: %s" % " ".join(call))
    try:
        results = tk_subproc.check_output(unicode_call)
        martian.log_info("dlconverter output: %s" % results)
    except subprocess.CalledProcessError, e:
        outs.output_for_dloupe = None
        martian.throw("Could not generate .dloupe file: \n%s" % e.output)
