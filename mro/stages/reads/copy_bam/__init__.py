#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Copy BAM file to output
#
import os
import shutil
import pysam

__MRO__ = """
stage COPY_BAM(
    in  bam     input,
    in  bam.bai input_index,
    out bam     possorted_bam,
    out bam.bai possorted_bam_index,
    src py      "stages/reads/copy_bam",
)
"""

def link_or_copy_file(src, dst):
    try:
        os.link(src, dst)
    except:
        shutil.copyfile(src, dst)

def main(args, outs):
    link_or_copy_file(args.duplicate_summary, outs.duplicate_summary)
    link_or_copy_file(args.input, outs.output)

    outs.output_index = outs.output+ ".bai"
    link_or_copy_file(args.input_index, outs.output_index)

    bam = pysam.AlignmentFile(outs.output, "rb", check_sq=False)
    total = bam.mapped + bam.unmapped
    print "mapped reads:   {} ({:.3g}%)".format(bam.mapped, 100.0*bam.mapped/total)
    print "unmapped reads: {} ({:.3g}%)".format(bam.unmapped, 100.0*bam.unmapped/total)
    bam.close()
