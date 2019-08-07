#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
# Create per-cell and per-analysis summaries for high-level metrics
#
import os
import json
import pandas as pd
import numpy as np

import martian
import tenkit.safe_json as tk_json
import tenkit.stats as tk_stats

import longranger.cnv.contig_manager as contig_manager
from crdna.utils import load_h5, load_mask_bdy

__MRO__ = """
stage REPORT_BASIC_SUMMARY(
    in  string   reference_path,
    in  csv      barnyard,
    in json      singlecell_summary,
    in  bed      node_cnv_calls,
    in  h5       norm_node_profiles,
    in  h5       node_cnv_tracks,
    in  json     report_basic,
    out csv      per_cell_summary_metrics,
    out json     summary,
    out csv      summary_cs,
    src py       "stages/reporter/report_basic_summary",
) 
"""

#
# These parameters are used to determine the num_altploidy_events_qv300 metric
#
DEFAULT_CONFIDENCE = 300.0
DEFAULT_PLOIDY = 2
#
# Only count flat clusters with this size or greater (expressed as a fraction of ncells)
#
CLUSTER_MIN_SIZE = 0.01

def split(args):
    martian.exit("No split")

def join(args, outs, chunk_defs, chunk_outs):
    martian.exit("No join")

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    stats = pd.read_csv(args.barnyard)
    stats = stats[stats['cell_id']!='None'].copy()
    ncells = len(stats)
    martian.log_info( 'Subsetting per-barcode statistics to %d cells' % ncells )

    ref = contig_manager.contig_manager(args.reference_path)
    contig_lengths = [ ref.contig_lengths[k] for k in ref.primary_contigs(allow_sex_chromosomes=True) ]
    tot_ref_len = float(sum(contig_lengths))
    martian.log_info( 'Reference sequence at %s has %d bp'%(args.reference_path,tot_ref_len) )

    #
    # Accumulate per-cell summary stats
    #

    PER_CELL_HEADER = ['barcode', 'cell_id', 'total_num_reads',
        'num_unmapped_reads', 'num_lowmapq_reads', 'num_duplicate_reads',
        'num_mapped_dedup_reads', 'frac_mapped_duplicates',
        'effective_depth_of_coverage', 'effective_reads_per_1Mbp',
        'raw_mapd', 'normalized_mapd', 'raw_dimapd',
        'normalized_dimapd', 'mean_ploidy',
        'ploidy_confidence', 'is_high_dimapd', 'is_noisy']

    # unusable reads are those that are non-cell barcodes that are also any of mapped, low mapq, nor dups
    # no barcode are reads whose barcode is not on the whitelist
    num_dups = stats['dups']
    num_lowmapq = stats['low_mapq_lt_30']
    num_unmapped = stats['no_barcode'] + stats['unusable_read'] + stats['unmapped']
    num_mapped = stats['mapped']

    assert all(num_unmapped + num_dups + num_lowmapq + num_mapped == stats['denominator'])

    per_cell = pd.DataFrame(columns=PER_CELL_HEADER)
    per_cell['barcode'] = stats['BC']
    per_cell['cell_id'] = np.arange(0,ncells)
    per_cell['num_mapped_dedup_reads'] = num_mapped
    per_cell['frac_mapped_duplicates'] = stats['dups_frac']
    per_cell['num_unmapped_reads'] = num_unmapped
    per_cell['num_lowmapq_reads'] = num_lowmapq
    per_cell['num_duplicate_reads'] = num_dups
    per_cell['total_num_reads'] = stats['denominator']
    per_cell['effective_depth_of_coverage'] = stats['num_mapped_bases'].astype(float)/tot_ref_len
    per_cell['effective_reads_per_1Mbp'] = np.round(stats['mapped']/(tot_ref_len/1e6)).astype(int)

    flat_metrics = load_h5(args.norm_node_profiles, 'normalization_metrics')
    per_cell['raw_mapd'] = flat_metrics['raw_mapd'].iloc[0:ncells].values
    per_cell['normalized_mapd'] = flat_metrics['norm_mapd'].iloc[0:ncells].values
    per_cell['raw_dimapd'] = flat_metrics['raw_dimapd'].iloc[0:ncells].values
    per_cell['normalized_dimapd'] = flat_metrics['norm_dimapd'].iloc[0:ncells].values

    mean_ploidy, num_altevents = process_cnv_metrics(ncells,args.node_cnv_calls,DEFAULT_CONFIDENCE)
    per_cell['mean_ploidy'] = mean_ploidy
    
    # per cell confidence score
    pconf = load_h5(args.node_cnv_tracks, "ploidy_conf").values[0:ncells]
    per_cell['ploidy_confidence'] = pconf

    # is noisy cell flag
    high_dimapd = flat_metrics['is_high_dimapd'].iloc[0:ncells].values
    per_cell['is_high_dimapd'] = high_dimapd
    # cells with low confidence, or cells whose ploidy estimate was 
    # overruled using high confidence cells
    low_ploidy_conf = np.logical_or(pconf == -4,
        (pconf >= 0) & (pconf <= 2))
    per_cell['is_noisy'] = np.logical_or(high_dimapd == 1, 
        low_ploidy_conf).astype(int)
    with open(outs.per_cell_summary_metrics,'w') as outfile:
        per_cell.to_csv(outfile,columns=PER_CELL_HEADER,index=False)

    #
    # Accumulate per-analysis summary stats
    #

    # combine mask vectors to calculate genome-wide mappability
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)
    masks, _ = load_mask_bdy(args.norm_node_profiles, chroms)

    with open(args.report_basic, 'r') as infile:
        report_basic = json.load(infile)

    with open(args.singlecell_summary, 'r') as infile:
        singlecell_summary = json.load(infile)

    per_analysis = {
        'total_num_bases_R1': report_basic['r1_tot_bases'],
        'total_num_bases_R1_Q30': report_basic['r1_q30_bases'],
        'total_num_bases_R2': report_basic['r2_tot_bases'],
        'total_num_bases_R2_Q30': report_basic['r2_q30_bases'],
        'frac_bases_R1_Q30': tk_stats.robust_divide(report_basic['r1_q30_bases'], report_basic['r1_tot_bases']),
        'frac_bases_R2_Q30': tk_stats.robust_divide(report_basic['r2_q30_bases'], report_basic['r2_tot_bases']),
        'total_num_reads': report_basic['num_reads'],
        'total_num_reads_in_cells': per_cell['total_num_reads'].sum(),
        'total_num_mapped_dedup_reads_in_cells': per_cell['num_mapped_dedup_reads'].sum(),
        'mean_mapped_dedup_reads_per_cell': per_cell['num_mapped_dedup_reads'].mean(),
        'median_frac_mapped_duplicates_per_cell': np.median(per_cell['frac_mapped_duplicates']),
        'num_cells': ncells,
        'median_effective_reads_per_1Mbp': np.median(per_cell['effective_reads_per_1Mbp']),
        'frac_mappable_bins': tk_stats.robust_divide(sum(masks), len(masks)),
        'frac_noisy_cells': tk_stats.robust_divide(sum(per_cell['is_noisy']), len(per_cell['is_noisy'])),
        'shortest_primary_contig': min(contig_lengths),
        'frac_non_cell_barcode': singlecell_summary['frac_waste_non_cell_barcode'],
        'correct_bc_rate': report_basic['correct_bc_rate'],
        'median_unmapped_frac': singlecell_summary['median_unmapped_frac']}

    for prefix in ["normalized", "raw"]:
        for metric in ["mapd", "dimapd"]:
            per_cell_key = "%s_%s"%(prefix, metric)
            for perc in [25, 50, 75]:
                summary_key = "%s_%s_p%d"%(prefix, metric, perc)
                per_analysis[summary_key] = tk_stats.robust_percentile(
                    per_cell[per_cell_key], perc)
    
    for cutoff in (25, 50, 75):
        k = 'mean_ploidy_p{:.2g}'.format(cutoff)
        per_analysis[k] = tk_stats.robust_percentile(per_cell['mean_ploidy'], cutoff)
    per_analysis['median_ploidy'] = per_analysis['mean_ploidy_p50']

    with open(outs.summary,'w') as outfile:
        outfile.write(tk_json.safe_jsonify(per_analysis, pretty=True)+os.linesep)
    SUMMARY_METRICS = [
        'total_num_reads',
        'frac_bases_R1_Q30',
        'frac_bases_R2_Q30',
        'correct_bc_rate',
        'frac_non_cell_barcode',
        'shortest_primary_contig',
        'frac_mappable_bins',
        'num_cells',
        'total_num_reads_in_cells',
        'total_num_mapped_dedup_reads_in_cells',
        'median_frac_mapped_duplicates_per_cell',
        'mean_mapped_dedup_reads_per_cell',
        'median_effective_reads_per_1Mbp',
        'median_unmapped_frac',
        'mean_ploidy_p25',
        'mean_ploidy_p50',
        'mean_ploidy_p75',
        'raw_mapd_p25',
        'raw_mapd_p50',
        'raw_mapd_p75',
        'normalized_mapd_p25',
        'normalized_mapd_p50',
        'normalized_mapd_p75',
        'normalized_dimapd_p25',
        'normalized_dimapd_p50',
        'normalized_dimapd_p75',
        'raw_dimapd_p25',
        'raw_dimapd_p50',
        'raw_dimapd_p75',
        'frac_noisy_cells'
    ]
    with open(outs.summary_cs,'w') as outfile:
        values = [per_analysis[key] for key in SUMMARY_METRICS]
        outfile.write(",".join(SUMMARY_METRICS) + "\n")
        outfile.write(",".join(map(str, values)) + "\n")

def process_cnv_metrics(ncells,cnv_calls_bed,min_confidence):
    """Calculates mean ploidy and number of alternative CNV events for
    the CNV calls in the provided BED file."""
    ploidies = np.zeros(ncells,dtype='float')
    tot_length = np.zeros(ncells,dtype='float')
    num_altevents = np.zeros(ncells,dtype='int')
    with open(cnv_calls_bed,'r') as infile:
        for l in infile:
            if l[0] == "#":
                continue
            values = l.rstrip().split('\t')
            cell_id = int(values[3])
            if cell_id>=ncells:
                continue
            start = int(values[1])
            end = int(values[2])
            length = end - start
            ploidy = int(values[4])
            conf = float(values[5])
            ploidies[cell_id] += ploidy * length
            tot_length[cell_id] += length
            if conf >= min_confidence and ploidy!=DEFAULT_PLOIDY:
                num_altevents[cell_id] += 1
    mean_ploidy = np.where( tot_length>0, ploidies/tot_length, 0.0 )
    return mean_ploidy, num_altevents
