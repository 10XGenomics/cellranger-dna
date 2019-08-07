#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Put BAM file into buckets by barcode prefix
#
import os.path
import shutil
import martian
import tenkit.bam as tk_bam
import crdna.bio_io as crdna_io
import tenkit.seq as tk_seq
from subprocess import check_call
from crdna.bam import BamTemplateShim
from crdna.constants import PROCESSED_BARCODE_TAG as BC_TAG

__MRO__ = """
stage BUCKET_BY_BC(
    in  int      nbases,
    in  bam      input,
    out map      buckets,
    out bam[]    non_bc_bams,
    src py       "stages/reads/bucket_reads_by_bc",
) split using (
    in  string   chunk_start,
    in  string   chunk_end,
    in  int      chunk_index,
)
"""

def split(args):
    bam_in = tk_bam.create_bam_infile(args.input)
    chunk_defs = tk_bam.chunk_bam_records(bam_in, chunk_bound_key=None, max_chunks=120)
    for i, chunk_def in enumerate(chunk_defs):
        chunk_def['chunk_index'] = i
        chunk_def['__mem_gb'] = 6
        chunk_def['__threads'] = 4
    return {'chunks': chunk_defs, 'join': {'__mem_gb': 1}}

def main(args, outs):
    chunk_start = args.chunk_start
    chunk_end = args.chunk_end
    chunk_index = args.chunk_index

    prefixes = get_seqs(args.nbases)
    bam_in = tk_bam.create_bam_infile(args.input)
    template = BamTemplateShim(bam_in, keep_comments=(chunk_index==0))

    bams_out = {}
    for prefix in prefixes:
        filename = martian.make_path("bc_{}.bam".format(prefix))
        bams_out[prefix], _ = tk_bam.create_bam_outfile(filename, None, None, template=template)

    non_bc_bam = martian.make_path("bc_{}.bam".format(None))
    non_bc_bam_out, _ = tk_bam.create_bam_outfile(non_bc_bam, None, None, template=template)
    for read in tk_bam.read_bam_chunk(bam_in, (chunk_start, chunk_end)):
        barcode = crdna_io.get_read_barcode(read)
        if barcode is None:
            non_bc_bam_out.write(read)
        else:
            prefix = barcode[:args.nbases]
            bams_out[prefix].write(read)
    bam_in.close()

    non_bc_bam_out.close()
    sort_bam(non_bc_bam)
    outs.non_bc_bams = [non_bc_bam]

    outs.buckets = {}
    for prefix in prefixes:
        filename = bams_out[prefix].filename
        bams_out[prefix].close()
        sort_bam(filename)
        outs.buckets[prefix] = filename

def get_seqs(l):
    if l == 1:
        return tk_seq.NUCS
    old_seqs = get_seqs(l-1)
    new_seqs = []
    for old_seq in old_seqs:
        for base in tk_seq.NUCS:
            new_seqs.append(old_seq + base)
    return new_seqs

def sort_bam(filename):
    tmp_bam, _ = os.path.splitext(filename)
    tmp_bam += '.sorted'
    check_call(['samtools', 'sort',
                '-@', '4',           # 4 threads
                '-n', '-t', BC_TAG,  # sort by tag then name
                '-o', tmp_bam, filename])
    shutil.move(tmp_bam, filename)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    outs.non_bc_bams = []
    outs.buckets = {}
    for out in chunk_outs:
        outs.non_bc_bams.extend(out.non_bc_bams)
        for prefix, filename in out.buckets.iteritems():
            if prefix not in outs.buckets:
                outs.buckets[prefix] = []
            outs.buckets[prefix].append(filename)
