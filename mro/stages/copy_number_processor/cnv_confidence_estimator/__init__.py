#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import json
import os
import pandas as pd
import shutil
import martian
import numpy as np

from longranger.cnv import contig_manager, coverage_matrix
from crdna.singlecell_dna_cnv import cnv_confidence_score as ccs
from crdna.singlecell_dna_cnv import vesna

__MRO__ = """
stage ESTIMATE_CNV_CONFIDENCE(
    in h5     raw_profiles,
    in h5     tracks,
    in string reference_path,
    in string sex,
    in bed    cnv_calls,
    in bed    unmerged_cnv_calls,
    in json   gc_norm_params,
    #
    out bed   cnvs,
    out bed   unmerged_cnvs,
    #
    src py    "stages/copy_number_processor/cnv_confidence_estimator",
) split using (
    in  map   chunk,
)       
"""

################################################################################
# Globals:
COLUMN_NAMES = ['Chr', 'Start', 'End', 'NodeID', 'CopyNumber', 'Confidence']
COLUMN_TYPES = [str, int, int, int, int, float]

################################################################################
## 
def split(args):
    with open(args.cnv_calls,'r') as infile:
        nodes = { l.rstrip().split('\t')[3] for l in infile }
    num_nodes = len(nodes)

    store = pd.HDFStore(args.raw_profiles)
    ref = contig_manager.contig_manager(args.reference_path)
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)
    max_chrom_nbins = max(store["/contigs/"+chrom].shape[1] for chrom in chroms)
    store.close()

    MAX_CHUNKS = 30
    MIN_NODES_PER_CHUNK = 5

    nchunks = np.clip(np.ceil(1.0*num_nodes/MIN_NODES_PER_CHUNK), 1, MAX_CHUNKS)
    nodes_per_chunk = max(1, int(np.ceil(1.0*num_nodes/nchunks)))

    chromsz_gb = 1.0 * max_chrom_nbins * max(1, num_nodes) / 1e9
    matsize_gb = (1.0 * coverage_matrix.get_genome_matrix_size_gb(args.raw_profiles) *
            nodes_per_chunk / max(1, num_nodes))
    unmerged_gb = int(np.ceil(os.path.getsize(args.unmerged_cnv_calls)/1e9))
    chunk_mem_gb = int(np.ceil(6*max(matsize_gb, chromsz_gb) + 2))
    join_mem_gb = int(np.ceil(6*unmerged_gb + 2))

    chunk_defs = [
        {
            'chunk':    {'start': i, 'end': min(i+nodes_per_chunk, num_nodes)},
            '__mem_gb': chunk_mem_gb
        }
        for i in xrange(0, num_nodes, nodes_per_chunk)
    ]

    return {'chunks': chunk_defs, 'join': {'__mem_gb': join_mem_gb}}

# split

################################################################################
def main(args, outs):
    """Compute a CNV confidence score from the profile for a specific choice of cluster
    and contig."""
    martian.log_info('Entering __init__.main()')
    node_start = args.chunk['start']
    # exclusive end
    node_end = args.chunk['end']
    raw_profiles, mask = coverage_matrix.load_matrix(args.raw_profiles,
        args.reference_path, start_cell=node_start, end_cell=node_end)

    bin_size = coverage_matrix.get_bin_size(args.raw_profiles)

    ## read in CNV data for nodes of interest
    node_column = COLUMN_NAMES.index("NodeID")
    cnv_calls = read_cnv_data(args.cnv_calls, node_start, node_end, node_column)
    #
    scale = get_scaling_factors(raw_profiles,cnv_calls)
    with open(args.gc_norm_params, "r") as handle:
        gc_norm_params = json.load(handle)
    linear = gc_norm_params["linear"]
    quadratic = gc_norm_params["quadratic"]
    #
    ref = contig_manager.contig_manager(args.reference_path)
    #
    # Get mappability, GC content:
    bin_parameters = []
    vesna.load_track_parameters(args.tracks, bin_parameters, ref)
    #
    logp, cnv_calls2 = ccs.process_cnv_calls(
        raw_profiles, mask, bin_parameters, 
        args.reference_path, args.sex,
        scale, linear, quadratic, cnv_calls,
        bin_size)

    export_segments(outs.cnvs, cnv_calls2, node_start)

    # free some memory
    del cnv_calls
    del cnv_calls2

    #
    # Compute confidence values for unmerged, broken-up CNV calls
    #
    unmerged_cnv_calls = read_cnv_data(args.unmerged_cnv_calls, node_start,
        node_end, node_column)
    
    _, unmerged_cnv_calls2 = ccs.process_cnv_calls(
        raw_profiles, mask, bin_parameters, 
        args.reference_path, args.sex,
        scale, linear, quadratic, unmerged_cnv_calls,
        bin_size, logp=logp)

    export_segments(outs.unmerged_cnvs, unmerged_cnv_calls2, node_start)
    # 
    # Debugging:
    #
    martian.log_info('Leaving __init__.main()')
    martian.log_info('.' * 80)
# main

################################################################################
def join(args, outs, chunk_defs, chunk_outs):
    """Join step.  Each chunk returns a bed file. This just concatenates
    the bed file from each chunk and sorts them."""
    #
    def write_sorted_bed(chunk_getter,outfilename):
        with open(outfilename,'w') as out_file:
            #
            for chunk in chunk_outs:
                if not os.path.exists(chunk_getter(chunk)):
                    continue
                # if !exists
                with open(chunk_getter(chunk),'r') as in_file:
                    shutil.copyfileobj(in_file, out_file, 1024 * 1024)
            # for chunk
        ref = contig_manager.contig_manager(args.reference_path)
        chroms = ref.primary_contigs(allow_sex_chromosomes=True)
        chrom_index = dict([(c, i) for i, c in enumerate(chroms)])
        
        cnv_df = pd.read_csv(outfilename, sep="\t", names=COLUMN_NAMES)
        cnv_df["chrom_index"] = cnv_df["Chr"].apply(chrom_index.get)
        cnv_df.sort_values(by=["chrom_index", "Start", "End"], inplace=True)
        cnv_df.to_csv(outfilename, sep="\t", columns=COLUMN_NAMES,
            header=False, index=False)

    ## write sorted cnvs.bed
    write_sorted_bed(lambda c:c.cnvs, outs.cnvs)
    write_sorted_bed(lambda c:c.unmerged_cnvs, outs.unmerged_cnvs)

# join

################################################################################
##
def read_cnv_data(bed_file, node_start, node_end, node_column):
    cnv_data = []
    with open(bed_file, "r") as cnv_in:
        for line in cnv_in:
            fields = line.strip().split("\t")
            row = map(lambda x, y: x(y), COLUMN_TYPES, fields)
            if row[node_column] < node_start or row[node_column] >= node_end:
                continue
            cnv_data.append(row)
    cnv_calls = pd.DataFrame(cnv_data, columns=COLUMN_NAMES)
    cnv_calls['NodeID'] = cnv_calls['NodeID'] - node_start
    return cnv_calls

def export_segments(file_name, segments, node_offset):
    output_file = open(file_name, 'a')
    for segment in segments.itertuples():
        ploidy = segment.CopyNumber
        chrom_name = segment.Chr
        node_id = segment.NodeID + node_offset
        start = segment.Start
        end = segment.End
        confidence = segment.Confidence
        try:
            out_text = '%s\t%i\t%i\t%s\t%i\t%d\n' % (
                chrom_name, start, end, node_id, ploidy, confidence)
            output_file.write(out_text)
        except Exception as error:
            print('export_segments() caught an error: %s' % repr(error))
            print('chrom_name, start, end, ploidy, confidence, node_id:')
            print(chrom_name, start, end, ploidy, confidence, node_id)
        # if ploidy
    # for segment
    output_file.close()
# export_segments

def get_scaling_factors(raw_profiles,cnv_calls):
    """Returns list of scaling factors in the same order
    as the cells (nodes) in raw_profiles."""
    assert len(raw_profiles)>0, "empty raw_profiles supplied to get_scaling_factors"
    avg_ploidy = get_node_avg_ploidy(cnv_calls)
    scaling_factors = []
    for node_id in xrange(raw_profiles[0].shape[0]):
        num, denom = 0.0, 0.0
        for p in raw_profiles:
            num += p[node_id,:].sum()
            denom += p.shape[1]
        mean_count = num / denom
        scaling_factors.append( mean_count / avg_ploidy[node_id] if avg_ploidy[node_id]>0.0 else 1.0 )
    return scaling_factors

def get_node_avg_ploidy(calls):
    """Returns dict of node_id -> avg. ploidy"""
    avg_ploidy = {}
    for node_id, rows in calls.groupby(calls['NodeID']):
        total_len = 0.0
        total_ploidy = 0.0
        for i,row in rows.iterrows():
            l = row['End']-row['Start']
            ploidy = row['CopyNumber']
            total_len += l
            total_ploidy += l * ploidy
        avg_ploidy[node_id] = total_ploidy/total_len
    return avg_ploidy
