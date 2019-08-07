#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import json
import tenkit.bam as tk_bam
import tenkit.alarms as tk_alarms
import crdna.bio_io as crdna_io
import numpy as np
import longranger.cnv.contig_manager as contig_manager
import martian
from crdna.constants import CELL_DETECT_MAPQ_THRESHOLD, MAX_CELLS, CRDNA_ALARMS

__MRO__ = """
stage REPORT_SINGLE_CELL(
    in  string reference_path,
    in  int[]  window_sizes,
    in  bam    input,
    in  string barcode_whitelist,
    in  int    force_cells,
    out json   summary,
    out csv    barnyard,
    out json   barcode_histogram,
    src py     "stages/reporter/report_singlecell",
) split using (
    in  string chunk_start,
    in  string chunk_end,
)
"""

def chunk_split_func(r):
    return crdna_io.get_read_barcode(r)

def split(args):
    if args.input is None or args.barcode_whitelist is None:
        chunk_defs = [{'chunk_start':"0", 'chunk_end':"0", '__mem_gb': 1}]
        return {'chunks': chunk_defs, 'join': {'__mem_gb': 1}}

    ref = contig_manager.contig_manager(args.reference_path)
    species_list = ref.list_species()
    if (args.force_cells is not None and args.force_cells > 0 and
        len(species_list) > 1):
        martian.exit("force_cells can only be used for single species reference.")
    min_chunks = 10
    bam_in = tk_bam.create_bam_infile(args.input)
    chunks = tk_bam.chunk_bam_records(bam_in, chunk_split_func, 
        chunk_size_gb=8.0, min_chunks=min_chunks)

    # 0.03 =~ 26meg = 1M bcs * (sizeof(int64) + 18)
    join_mem_gb = int(np.ceil(0.03*(len(chunks) + 1) + 1))
    return {'chunks': chunks, 'join': {'__mem_gb': join_mem_gb}}

def main(args, outs):
    ref = contig_manager.contig_manager(args.reference_path)
    args.coerce_strings()
    outs.coerce_strings()
    
    # Bail out if there no valid barcodes
    if args.barcode_whitelist is None or args.input is None:
        outs.summary = None
        return

    bam_in = tk_bam.create_bam_infile(args.input)
    bam_chunk = tk_bam.read_bam_chunk(bam_in,(args.chunk_start, args.chunk_end))
    
    # Skip reads without a barcode
    bam_chunk_filt = itertools.ifilter(read_has_barcode, bam_chunk)
    bc_read_iter = itertools.groupby(bam_chunk_filt,
        lambda x: crdna_io.get_read_barcode(x))

    counts = {}

    for bc, reads in bc_read_iter:
        for r in reads:
            contig = bam_in.references[r.tid]
            species = ref.species_from_contig(contig)
            if not species in counts:
                counts[species] = {}
            if not bc in counts[species]:
                counts[species][bc] = 0
            if r.is_secondary or r.is_supplementary:
                ## we are ignoring alternate alignments
                continue
            if (r.is_unmapped or
                r.mapping_quality < CELL_DETECT_MAPQ_THRESHOLD
                or r.is_duplicate):
                ## if read is unmapped, poor mapping quality or dup
                continue
            counts[species][bc] += 1
    outs.counts = counts

def join(args, outs, chunk_defs, chunk_outs):
    cell_barcodes = {}
    full_counts = {}
    for chunk_out in chunk_outs:
        if not chunk_out.counts:
            continue
        for (species, bc_counts) in chunk_out.counts.iteritems():
            if not species in full_counts:
                full_counts[species] = {}
            for (bc, count) in bc_counts.iteritems():
                full_counts[species][bc] = count

    for species in full_counts.iterkeys():
        sorted_data = sorted(full_counts[species].items(),
            key=lambda item: item[1], reverse=True)
        bc_counts = [item[1] for item in sorted_data]

        ## Define cells by first taking one log width. Then refine by choosing
        ## one log width from the 99th percentile amongst the cells
        max_count = bc_counts[0]
        min_count = float(max_count)/np.power(10, args.log_width)

        count_99th = np.percentile([x for x in bc_counts if x >= min_count], 99)
        min_count = float(count_99th)/np.power(10, args.log_width)

        # implement force_cells if supplied to overrule min_count
        if args.force_cells is not None and args.force_cells > 0:
            index = min(args.force_cells, len(bc_counts))-1
            min_count = max(bc_counts[index], 1)
            martian.log_info("Using force_cells")

        if not species in cell_barcodes:
            cell_barcodes[species] = {}
        for i, (bc, count) in enumerate(sorted_data, start=1):
            if count < min_count:
                break
            cell_barcodes[species][bc] = count
            if i >= MAX_CELLS:
                martian.log_info("%s: hit maximum number of cells "\
                    "(%d)"%(species, MAX_CELLS))
                min_count = count
                break

        martian.log_info("%s: max count %d, min count %d"%(species, max_count,
            min_count))

        # some logging
        ncell = len(cell_barcodes[species])
        nobs = len(bc_counts)
        if len(cell_barcodes[species]) > 0:
            mean = np.mean(cell_barcodes[species].values())
            median = np.median(cell_barcodes[species].values())
        else:
            mean = 0
            median = 0
        print ("{}: {} cells of {} obs, cell barcode reads: "
                "mean = {:.2f}, median = {:.1f}").format(
                species, ncell, nobs, mean, median)

    # alarm user
    with open(CRDNA_ALARMS) as af:
        alarms = json.load(af)
    # filter alarms
    alarms = [alarm for alarm in alarms if alarm['id'] in ['not_enough_cells',
        'too_many_cells']]
    num_cells = sum([len(bc_cts) for bc_cts in cell_barcodes.itervalues()])
    alarm_results = tk_alarms.evaluate_alarms(alarms,
        {'num_cells': num_cells})
    for a in alarm_results:
        martian.alarm("%s is %s. %s\n" % (a["title"], a["formatted_value"],
            a["message"]))

    outs.cell_barcodes = cell_barcodes

def read_has_barcode(r):
    bc = crdna_io.get_read_barcode(r)
    if bc is None:
        return False
    else:
        return True

