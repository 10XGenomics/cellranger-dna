#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Sort BAM file by sorting buckets, then concatenating bucket files
#

from __future__ import absolute_import, division, print_function

import numpy as np
import tenkit.bam as tk_bam

from crdna.constants import PROCESSED_BARCODE_TAG as BC_TAG
from tenkit.bam import merge_by_tag

__MRO__ = """
stage SORT_BY_BC(
    in  map     bc_buckets,
    in  bam[]   non_bc_bams,
    in  bam     possorted_bam,
    out int     total_reads,
    out bam     bc_sorted_bam,
    src py      "stages/reads/sort_reads_by_bc",
) split using (
    in  string  prefix,
    in  bam[]   bucket,
)
"""

def mem_for_bucket(bucket):
    gb_per_samfile = 0.012  # 12mb per Samfile
    return max(1, int(np.ceil(len(bucket) * gb_per_samfile + 0.5)))

def split(args):
    chunk_defs = [{'prefix': None, 'bucket': args.non_bc_bams,
                  '__mem_gb': mem_for_bucket(args.non_bc_bams), '__threads': 4}]
    for prefix, bucket in args.bc_buckets.iteritems():
        chunk_defs.append({'prefix': prefix, 'bucket': bucket,
                           '__mem_gb': mem_for_bucket(bucket), '__threads': 4})
    return {'chunks': chunk_defs}

def main(args, outs):
    outs.coerce_strings()
    merge_by_tag(outs.bcsorted_bam, args.bucket, tag=BC_TAG, name=True, threads=4)

def join(args, outs, chunk_defs, chunk_outs):
    buckets = []
    for _, out in sorted(zip(chunk_defs, chunk_outs), key=lambda x: x[0].prefix):
        buckets.append(out.bcsorted_bam)
    tk_bam.concatenate(outs.bcsorted_bam, buckets)
