#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Compute confident windows from BED file
#

import tenkit.stats as tk_stats
import tenkit.bio_io as tk_io
import tenkit.regions as tk_regions
from tenkit.chunk_utils import generate_chrom_loci, pack_loci 
import tenkit.reference
import longranger.cnv.contig_manager as contig_manager

__MRO__ = """
stage REPORT_CONFIDENT_WINDOWS(
    in  string reference_path,
    in  bed    confident_regions,
    in  int    window_size,
    out bed    confident_windows,
    src py     "stages/reporter/report_confident_windows",
) split using (
    in  string[] loci,
)
"""

def split(args):
    ref = contig_manager.contig_manager(args.reference_path)
    contig_lengths = ref.get_contig_lengths( )

    target_regions = None
    all_loci = []
    for (chrom_name, chrom_size) in contig_lengths.iteritems():
        all_loci.extend(generate_chrom_loci(target_regions, chrom_name, chrom_size,
            tenkit.constants.PARALLEL_LOCUS_SIZE))

    locus_sets = pack_loci(all_loci)

    chunk_defs = [{'loci': loci, '__mem_gb': 12} for loci in locus_sets]
    return {'chunks': chunk_defs, 'join': {'__mem_gb': 12}}


def join(args, outs, chunk_defs, chunk_outs):
    outfile = open(outs.confident_windows, 'w')
    for chunk_out in chunk_outs:
        with open(chunk_out.confident_windows) as infile:
            for line in infile:
                outfile.write(line)
    outfile.close( )

def main(args, outs):
    args.coerce_strings()
    outs.coerce_strings()

    if args.confident_regions is None:
        confident_regions = None
    else:
        confident_regions = tk_io.get_target_regions(open(args.confident_regions))
    
    outfile = open(outs.confident_windows, "w")
    for (chrom, start, end) in (tk_io.get_locus_info(l) for l in args.loci):
        conf_regions = get_conf_regions(chrom, confident_regions)
        location = start
        while location < end:
            region = tk_regions.Regions(regions=[(location, location+args.window_size)])
            isect = region.intersect(conf_regions)
            size = isect.get_total_size()
            percent = tk_stats.robust_divide(float(size), float(args.window_size))
            row = [chrom, location, location+args.window_size, percent]
            outfile.write("\t".join(map(str, row)) + "\n")
            location += args.window_size
    outfile.close( )

def get_conf_regions(chrom, conf_regions_genome):
    if conf_regions_genome == None:
        return tk_regions.Regions([(0,10**10)])
    return tk_regions.Regions(conf_regions_genome.get(chrom, []))
