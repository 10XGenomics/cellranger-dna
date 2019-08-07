#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Make a barcoded-sorted bam file from a bunch of sorted buckets.
#

import tenkit.bam as tk_bam
import os

__MRO__ = """
stage MERGE_BC_BAM(
    in  bam      bc_sorted_bams,
    out bam      bc_sorted_bam,
    src py       "stages/reads/merge_bc_bam",
)
"""

def main(args, outs):
    filtered_inputs = [a for a in args.bc_sorted_bams if os.path.isfile(a)]
    tk_bam.concatenate(outs.bc_sorted_bam, filtered_inputs)
