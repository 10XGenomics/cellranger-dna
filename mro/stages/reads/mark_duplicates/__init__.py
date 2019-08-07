#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Mark PCR duplicates in a BAM file
#
from crdna.duplicates import DupSummary, broadcast

import json
import itertools
import tenkit.dict_utils
import crdna.bio_io as crdna_io
import tenkit.bam as tk_bam
import tenkit.lane as tk_lane
import tenkit.coverage

import martian
import os
from crdna.bam import BamTemplateShim
from crdna.constants import SELF_FIVE_PRIME_POS_TAG

__MRO__ = """
stage MARK_DUPLICATES(
    in  bam     input,
    in  int     perfect_read_count,
    in  bed     targets_file,
    in  json    diffusion_dup_summary,
    in  bool    write_bam,
    out bam     output,
    out bam.bai index,
    out json    duplicate_summary,
    src py      "stages/reads/mark_duplicates",
) split using (
    in  map     lane_map,
    in  string  chunk_start,
    in  string  chunk_end,
    in  int     chunk_index,
    in  float   diffusion_threshold,
)
"""

def chunk_bound_func(read):
    # Since the reads are sorted by SELF_FIVE_PRIME_POS_TAG tag, use that
    # for the chunk boundary
    if not read.is_unmapped:
        return (read.reference_id, read.get_tag(SELF_FIVE_PRIME_POS_TAG))
    else:
        return None


def split(args):
    # Chunk bam to get 1GB per chunk
    bam_in = tk_bam.create_bam_infile(args.input)
    chunk_defs = tk_bam.chunk_bam_records(bam_in, chunk_bound_func, chunk_size_gb=0.75)

    for i, chunk in enumerate(chunk_defs):
        chunk['chunk_index'] = i
        chunk['__mem_gb'] = 3

    lane_coord_sys = tk_lane.LaneCoordinateSystem()

    # Reopen BAM for estimating tile extents
    bam_in = tk_bam.create_bam_infile(args.input)
    lane_coord_sys.estimate_tile_extents(bam_in)

    with open(args.diffusion_dup_summary) as f:
        data = json.load(f)
        threshold = data['diffusion']['threshold']

    for chunk in chunk_defs:
        chunk['lane_map'] = lane_coord_sys.to_dict()
        chunk['diffusion_threshold'] = threshold

    return {'chunks': chunk_defs, 'join': {'__mem_gb': 1, '__threads': 6}}

def main(args, outs):
    """
    Mark exact duplicate reads in the BAM file. Duplicates have the same read1 start site and read2 start site
    """

    lane_coord_sys = tk_lane.LaneCoordinateSystem.from_dict(args.lane_map)

    args.coerce_strings()
    outs.coerce_strings()

    bam_in = tk_bam.create_bam_infile(args.input)
    template = BamTemplateShim(bam_in, keep_comments=(args.chunk_index==0))
    
    if args.write_bam:
        bam_prefix, ext = os.path.splitext(outs.output)
        out_bam_name = bam_prefix + '_five_prime_pos_sorted' + ext
        bam_out, _ = tk_bam.create_bam_outfile(out_bam_name, None, None, template=template,
                                               pgs=[tk_bam.make_pg_header(martian.get_pipelines_version(),
                                                                          "mark_duplicates")])
        outs.index = None # chunk bams don't get indexed
    else:
        bam_out = None
        outs.output = None
        outs.index = None

    # Determine whether the BAM has 10x barcodes
    bam_in.reset()
    has_barcodes = [crdna_io.read_has_barcode(x) for x in itertools.islice(bam_in, 1000)]
    have_barcodes = (float(sum(has_barcodes)) / len(has_barcodes)) > 0.1

    # All read duplicate marking - these dup decisions are written to bam_out
    # the output bam has BC aware dup marking if available.
    # Ensure the summary key indicates what kind of dup marking was actually performed.
    if have_barcodes:
        no_filter_dups_bcs =    DupSummary(False, 1.0, True,  "no_filter_full_use_bcs", lane_coord_sys, output_bam=bam_out, threshold=args.diffusion_threshold)
        no_filter_dups_no_bcs = DupSummary(False, 1.0, False, "no_filter_full_ignore_bcs", lane_coord_sys, threshold=args.diffusion_threshold)
    else:
        no_filter_dups_bcs =    DupSummary(False, 1.0, True,  "no_filter_full_use_bcs", lane_coord_sys, threshold=args.diffusion_threshold)
        no_filter_dups_no_bcs = DupSummary(False, 1.0, False, "no_filter_full_ignore_bcs", lane_coord_sys, output_bam=bam_out, threshold=args.diffusion_threshold)


    # Dup marking on all perfect reads
    full_dups_bcs = DupSummary(True, 1.0, True, "full_use_bcs", lane_coord_sys, threshold=args.diffusion_threshold, tag_counts=True)
    full_dups_no_bcs = DupSummary(True, 1.0, False, "full_ignore_bcs", lane_coord_sys, threshold=args.diffusion_threshold)

    dup_sums = [full_dups_bcs, full_dups_no_bcs, no_filter_dups_bcs, no_filter_dups_no_bcs]

    # Now broadcast the selected reads to the summarizers
    # We can't do the points the require a sample_rate > 1.0 so, skip those.
    # If we don't have barcodes, don't run the set that are split by barcode.
    consumers = [x.read_consumer() for x in dup_sums if x.sample_rate <= 1.0 and ((not x.split_bcs) or have_barcodes)]

    source = tk_bam.read_bam_chunk(bam_in, (args.chunk_start, args.chunk_end))
    broadcast(source, consumers)

    # We close the BAM
    if bam_out:
        bam_out.close()
        # Note - the indexing happens in join
        bam_prefix, _ = os.path.splitext(outs.output)
        tk_bam.sort(out_bam_name, bam_prefix)

    # Package up the summaries:
    dup_results = {}
    for x in dup_sums:
        (dups, optical_dups, diff_dups, custom_diff_dups) = x.result
        desc = x.description
        dup_results[desc] = dups
        optical_desc = "optical_" + desc
        dup_results[optical_desc] = optical_dups
        diff_desc = "diffusion_old_" + desc
        dup_results[diff_desc] = diff_dups
        custom_diff_desc = "diffusion_" + desc
        dup_results[custom_diff_desc] = custom_diff_dups

    if outs.duplicate_summary:
        with open(outs.duplicate_summary, 'w') as f:
            json.dump(dup_results, f, indent=4)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    # combine the duplicate summary counts
    dup_summaries = [json.load(open(out.duplicate_summary)) for out in chunk_outs]
    combined_dups = reduce(lambda x,y: tenkit.dict_utils.add_dicts(x,y,2), dup_summaries)

    diffusion_summary = json.load(open(args.diffusion_dup_summary))

    combined_dups['read_counts'] = {}
    combined_dups['read_counts']['perfect_read_count'] = args.perfect_read_count

    for k, v in diffusion_summary.items():
        combined_dups[k] = v

    # TODO: Remove null_* observed_* ?

    with open(outs.duplicate_summary, 'w') as f:
        json.dump(combined_dups, f, indent=4)

    # combine & index the chunks of the BAM
    if args.write_bam:
        tk_bam.merge(outs.output, [c.output for c in chunk_outs], args.__threads)
        tk_bam.index(outs.output)
        outs.index = outs.output + '.bai'
    else:
        outs.output = None
        outs.index = None
