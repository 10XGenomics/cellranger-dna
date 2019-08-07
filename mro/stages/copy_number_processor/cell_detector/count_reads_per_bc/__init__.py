#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import tenkit.bam as tk_bam
import crdna.bio_io as crdna_io
import barcodes.utils as bc_utils

__MRO__ = """
stage COUNT_READS_PER_BC(
    in  bam    bcsorted_bam,
    in  string barcode_whitelist,
    out json   reads_per_barcode,
    src py     "stages/copy_number_processor/cell_detector/count_reads_per_bc",
) split using (
    in  string chunk_start,
    in  string chunk_end,
)
"""

#...............................................................................
def chunk_split_func(r):
    return crdna_io.get_read_barcode(r)
# chunk_split_func

#...............................................................................
def split(args):
    if args.bcsorted_bam is None or args.barcode_whitelist is None:
        chunk_defs = [{'chunk_start':"0", 'chunk_end':"0"}]
        return {'chunks': chunk_defs}
    # if args.input

    # Some R&D bc sets have very small diversity -- don't run on them
    barcode_whitelist = bc_utils.load_barcode_whitelist(args.barcode_whitelist)
    if len(barcode_whitelist) < 100:
        chunk_defs = [{'chunk_start':"0", 'chunk_end':"0"}]
        return {'chunks': chunk_defs}
    # if barcode_whitelist

    min_chunks = 4
    if len(barcode_whitelist) > 1e6:
        min_chunks = 8
    # if barcode_whitelist

    bam_in = tk_bam.create_bam_infile(args.bcsorted_bam)
    chunks = tk_bam.chunk_bam_records(bam_in, chunk_split_func, chunk_size_gb=8.0, min_chunks=min_chunks)
    for c in chunks:
        c['__mem_gb'] = 12
    # for c

    return {'chunks': chunks, 'join': {'__mem_gb': 32}}
# split

#...............................................................................
def main(args, outs):
    print("detect_cells args: barcode_whitelist=%s\n" % args.barcode_whitelist)
    pass
#    min_insert_size = 0
#    max_insert_size = 1e4
#    # TODO: Make this a variable or infer automatically?
#    sc_ratio = 100.0
#    # TODO: Add human genome size to tenkit/constants.py
#    human_genome_size = 3.1e9
#    
#    unique_min_mapq = 60.0
#
#    args.coerce_strings()
#    outs.coerce_strings()
#
#    window_sizes = args.window_sizes
#
#    # Bail out if there no valid barcodes
#    if args.barcode_whitelist is None or args.input is None:
#        outs.summary = None
#        return
#
#    bam_in = tk_bam.create_bam_infile(args.input)
#    bam_chunk = tk_bam.read_bam_chunk(bam_in, (args.chunk_start, args.chunk_end))
#    # Skip reads without a barcode
#    bam_chunk_filt = itertools.ifilter(read_has_barcode, bam_chunk)
#    bc_read_iter = itertools.groupby(bam_chunk_filt, lambda x: tk_io.get_read_barcode(x))
#
#    has_species_info = False
#    refs = bam_in.references
#    if '_' in refs[0]:
#        has_species_info = True
#
#    # For each barocode, count # per each contig, number per each window (for each window size)
#    # number per species (if available in contig), number per species
#    # TODO: Add detailed matrix by contigs, windows output
#    num_sc_bcs = 0
#    num_qual_reads = 0
#    num_sc_reads = 0
#    sc_counts = []
#
#    barnyard_file = open(outs.barnyard, 'w')
#    all_species = set()
#    if has_species_info:
#        for ref in refs:
#            all_species.add(ref.split('_')[0])
#            
#    barnyard_file.write(','.join(['BC'] + 
#        [s + '.reads' for s in list(all_species)] + 
#        [s + '.contigs' for s in list(all_species)] + 
#        ['mapped'] +
#        ['low.mapq'] +
#        ['unmapped'] +
#        ['duplicates'] +
#        ['mapped.fraction'] +
#        ['amp.rate'] +
#        ['dup.ratio'] +
#        ['amplicon.length'] +
#        ['num.pairs']
#        ) + '\n')
#
#    paired_multiplier = 2
#    ploidy = 2
#    bc_hist = {}
#    for bc, reads in bc_read_iter:
#        bc_count = 0
#        num_per_species = {}
#        contigs_per_species = {}
#        num_per_contig = {}
#        num_per_windows = {window_size:{} for window_size in window_sizes}
#        num_unmapped = 0
#        num_low_mapq = 0
#        num_dups = 0
#        amplicon_length = 0.0
#        num_pairs = 0;
#        for r in reads:
#            bc_count += 1
#
#            if r.is_unmapped: 
#                num_unmapped += 1
#                continue
#            elif not(r.mapq >= unique_min_mapq): 
#                num_low_mapq += 1
#                continue
#
#            if r.is_duplicate:
#                num_dups += 1
#
#            if r.is_paired:
#                delta = abs(r.reference_start - r.next_reference_start)
#                if ((r.reference_id == r.next_reference_id) and 
#                    (delta > min_insert_size) and
#                    (delta < max_insert_size)):
#                    insert_size = tk_bam.get_insert_size(r)
#                    if insert_size != float('Inf'): 
#                        amplicon_length += insert_size
#                        num_pairs += 1
#        
#            contig = refs[r.tid]
#            if has_species_info:
#                species = contig.split('_')[0]
#                num_per_species[species] = num_per_species.get(species, 0) + 1
#                species_contigs = contigs_per_species.setdefault(species, set())
#                species_contigs.add(contig)
#                
#            num_per_contig[contig] = num_per_contig.get(contig, 0) + 1
#
#            for window_size in window_sizes:
#                window = r.pos/window_size
#                num_per_windows[window_size][(contig, window)] = num_per_windows[window_size].get((contig, window), 0) + 1
#
#        num_pairs = num_pairs / 2 # because every pair is counted twice, once for each read in the pair
#        amplicon_length = amplicon_length / 2 # because every pair contributes the same insert size twice, once for each read in the pair
#        if num_pairs > 0:
#            amplicon_length = float(amplicon_length) / float(num_pairs)
#        else:
#            amplicon_length = float('NaN')
#
#        #print("ZDZ DEBUG: num_pairs=%d, amplicon_length=%f" % (num_pairs, amplicon_length))
#
#        num_mapped = sum(num_per_contig.values())
#        denominator = float(num_mapped - num_dups)
#        if denominator > 0:
#            dup_ratio = float(num_dups) / denominator
#        else:
#            dup_ratio = float('NaN')
#        # TODO: this is a rough amp rate estimate:
#        denominator = float(2 * dup_ratio * ploidy * human_genome_size)
#        if denominator > 0:
#            num_amplicons = num_mapped / paired_multiplier
#            amp_rate = float(num_amplicons * amplicon_length) / denominator
#            # amp_rate = float(num_amplicons * 250 ) / denominator
#            #print("ZDZ DEBUG: num_amplicons=%d, num_pairs=%d, num_mapped=%d" % (num_amplicons, num_pairs, num_mapped))
#        else:
#            amp_rate = float('NaN')
#        bc_hist[bc] = bc_count
#        num_qual_reads += num_mapped
#        denominator = float(num_unmapped + num_low_mapq + num_mapped)
#        if denominator > 0:
#            map_rate = float(num_mapped) / denominator
#        else:
#            map_rate = float('NaN')
#
#        barnyard_file.write(','.join([bc] + 
#            [str(num_per_species.get(s, 0)) for s in all_species] + 
#            [str(len(contigs_per_species.get(s, set()))) for s in all_species] + 
#            [str(num_mapped)] +
#            [str(num_low_mapq)] +
#            [str(num_unmapped)] +
#            [str(num_dups)] +
#            [str(map_rate)] +
#            [str(amp_rate)] +
#            [str(dup_ratio)] +
#            [str(amplicon_length)] +
#            [str(num_pairs)]
#            ) + '\n')
#        if has_species_info:
#            counts_by_species = sorted(num_per_species.values(), reverse=True)
#            if len(counts_by_species) < 2:
#                species_ratio = float('NaN')
#            else:
#                species_ratio = float(counts_by_species[0])/float(counts_by_species[1])
#
#            if species_ratio >= sc_ratio:
#                num_sc_bcs += 1
#                num_sc_reads += sum(num_per_contig.values())
#
#                num_contigs = len(num_per_contig.keys())
#                num_windows = {}
#                for window_size in window_sizes:
#                    num_windows[window_size] = len(num_per_windows.get(window_size, {}).keys())
#
#                sc_counts.append((num_contigs, num_windows))
#                
#    # Summarize (among single cell barcodes)
#    # For now, only infer single cells from barnyard.
#    # TODO: Use non-barnyard single cell inference
#    # Summary metrics:
#    # - Total number qualifying reads
#    # - Number barcodes appear barnyard single cell (if species info available)
#    # - Number of reads in banyard single cell
#    # - Total number of contigs covered by at least 1 read across all sc bcs
#    # - Total number of windows covered by at least 1 read for each window size across all sc bcs
#    tot_sc_contigs = 0.0
#    tot_sc_windows = {}
#    for num_contigs, num_windows in sc_counts:
#        tot_sc_contigs += num_contigs
#
#        for window_size in window_sizes:
#            tot_sc_windows[window_size] = tot_sc_windows.get(window_size, 0.0) + num_windows.get(window_size, 0.0)
#
#    summary_info = {}
#    summary_info['num_sc_bcs'] = num_sc_bcs
#    summary_info['num_sc_qual_reads'] = num_qual_reads
#    summary_info['num_sc_reads'] = num_sc_reads
#    summary_info['tot_sc_contigs'] = tot_sc_contigs
#    summary_info['tot_sc_windows'] = tot_sc_windows
#
#    barnyard_file.close()
#    
#    with open(outs.summary, 'w') as summary_file:
#        summary_file.write(json.dumps(summary_info))
#
#    with open(outs.barcode_histogram, 'w') as bc_hist_file:
#        bc_hist_file.write(json.dumps(bc_hist))

def join(args, outs, chunk_defs, chunk_outs):
    pass
#    num_sc_bcs = 0
#    num_qual_reads = 0
#    num_sc_reads = 0
#    tot_sc_contigs = 0.0
#    tot_sc_windows = {}
#
#    barnyard_file = open(outs.barnyard, 'w')
#
#    bc_counts = {}
#    for j,chunk_out in enumerate(chunk_outs):
#        if chunk_out.summary is None: continue
#
#        chunk_barnyard_file = open(chunk_out.barnyard)
#        if not(j == 0):
#            chunk_barnyard_file.readline()
#
#        for line in chunk_barnyard_file:
#            barnyard_file.write(line)
#
#        chunk_summary = json.loads(open(chunk_out.summary).read())
#        num_sc_bcs += chunk_summary['num_sc_bcs']
#        num_qual_reads += chunk_summary['num_sc_qual_reads']
#        num_sc_reads += chunk_summary['num_sc_reads']
#        tot_sc_contigs += chunk_summary['tot_sc_contigs']
#        
#        chunk_sc_windows = chunk_summary['tot_sc_windows']
#        for window_size in args.window_sizes:
#            tot_sc_windows[window_size] = tot_sc_windows.get(window_size,0) + chunk_sc_windows.get(str(window_size),0)
#
#        chunk_bc_counts_file = open(chunk_out.barcode_histogram)
#        chunk_bc_counts = json.loads(chunk_bc_counts_file.read())
#        bc_counts.update(chunk_bc_counts)
#
#    n_reads = np.array(bc_counts.values())
#    max_val = np.percentile(n_reads, 99.99) * 1.3
#    min_val = n_reads.min()
#    num_bins = 400
#    step = math.ceil((max_val - min_val)/num_bins)
#    bins = np.arange(min_val, max_val, step)
#    (hist, edges) = np.histogram(n_reads, bins=bins)
#    bc_hist = {int(edges[i]):hist[i] for i in range(len(bins)-1)}
#                                
#    mean_sc_contigs = tk_stats.robust_divide(tot_sc_contigs, num_sc_bcs)
#    mean_sc_windows = {}
#    for window_size in args.window_sizes:
#        mean_sc_windows[window_size] = tk_stats.robust_divide(tot_sc_windows.get(window_size,0), num_sc_bcs)
#
#    summary_info = {}
#    summary_info['num_sc_bcs'] = num_sc_bcs
#    summary_info['num_sc_qual_reads'] = num_qual_reads
#    summary_info['num_sc_reads'] = num_sc_reads
#    summary_info['fract_sc_reads'] = tk_stats.robust_divide(num_sc_reads, num_qual_reads)
#    summary_info['mean_sc_contigs'] = mean_sc_contigs
#    summary_info['mean_sc_windows'] = mean_sc_windows
#
#    barnyard_file.close()
#    with open(outs.summary, 'w') as summary_file:
#        summary_file.write(json.dumps(summary_info))
#
#    with open(outs.barcode_histogram, 'w') as bc_hist_file:
#        bc_hist_file.write(json.dumps(bc_hist))
# main

##...............................................................................
def read_has_barcode(r):
    bc = crdna_io.get_read_barcode(r)
    if bc is None:
        return False
    else:
        return True
# read_has_barcode
