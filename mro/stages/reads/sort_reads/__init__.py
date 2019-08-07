#!/usr/bin/env python
#
# Copyright (c) 2015-18 10X Genomics, Inc. All rights reserved.
#
# Sort BAM file
#
import os.path as op
import tenkit.bam as tk_bam
import crdna.read_filter
from tenkit.log_subprocess import check_call
from crdna.constants import SELF_FIVE_PRIME_POS_TAG
from tenkit.bam import merge_by_tag

__MRO__ = """
stage SORT_BY_POS(
    in  bam[] inputs,
    out bam,
    src py  "stages/reads/sort_reads",
) split using (
    in  bam chunk_input,
)
"""

def split(args):
    chunk_defs = [{'chunk_input': x} for x in args.inputs if op.isfile(x)]
    return {'chunks': chunk_defs, 'join': {'__threads': 6}}

def main(args, outs):
    args.coerce_strings()
    bam_prefix, ext = op.splitext(outs.default)

    # Sort based on the five prime position tag
    sort_args = ["samtools", "sort",
                 "-t", SELF_FIVE_PRIME_POS_TAG,
                 "-o", outs.default,
                 args.chunk_input]
    check_call(sort_args)

    perfect_read_count = 0
    bam = tk_bam.create_bam_infile(str(args.chunk_input))
    while True:
        try:
            read = bam.next()
            if crdna.read_filter.stringent_read_filter(read, True):
                perfect_read_count += 1
        except StopIteration:
            break
    outs.perfect_read_count = perfect_read_count

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.default) for chunk in chunk_outs]
    merge_by_tag(outs.default, input_bams, tag=SELF_FIVE_PRIME_POS_TAG, threads=int(args.__threads))
    outs.perfect_read_count = sum([chunk.perfect_read_count for chunk in chunk_outs])
