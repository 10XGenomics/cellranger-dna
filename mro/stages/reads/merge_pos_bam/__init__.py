#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Make a position-sorted bam file from a bunch of sorted buckets.
#

import subprocess
import json
import heapq
import tenkit.bam as tk_bam
from sets import Set
import tenkit.reference as tk_reference
__MRO__ = """
stage MERGE_POS_BAM(
    in  json     position_chunks,
    out bam      pos_sorted_bam,
    out bam.bai  pos_sorted_bam_index,
    src py       "stages/reads/merge_pos_bam",
) split using (
    in  string[] position_chunk,
)
"""
def split(args):
    chunk_defs = []
    with open(args.position_chunks) as position_chunks_file:
        position_chunks = json.load(position_chunks_file)
    pos_chunks = {int(key):val for (key,val) in position_chunks.iteritems()}
    for key, position_chunk in sorted(pos_chunks.iteritems()):
        chunk_defs.append({'position_chunk':position_chunk, "__mem_gb": 6})
    return {'chunks':chunk_defs}

def main(args, outs):
    samtools_merge_args = ['samtools','merge',outs.pos_sorted_bam+"tmp.bam"]
    samtools_merge_args.extend(args.position_chunk)
    subprocess.check_call(samtools_merge_args)
    update_mapqs(outs.pos_sorted_bam+"tmp.bam", outs.pos_sorted_bam+"tmp2.bam", args.reference_path)
    subprocess.check_call(['rm',outs.pos_sorted_bam+"tmp.bam"])
    subprocess.check_call(['samtools','sort', outs.pos_sorted_bam+"tmp2.bam", outs.pos_sorted_bam[:-4]])
    subprocess.check_call(['rm',outs.pos_sorted_bam+"tmp2.bam"])

def join(args, outs, chunk_defs, chunk_outs):

    samtools_cat_args = ['samtools','cat']
    position_bams = [chunk.pos_sorted_bam for chunk in chunk_outs]
    samtools_cat_args.extend(position_bams)
    with open(outs.pos_sorted_bam,'w') as outfile:
        subprocess.check_call(samtools_cat_args,stdout=outfile)

    subprocess.check_call(['samtools','index',outs.pos_sorted_bam])
    outs.pos_sorted_bam_index = outs.pos_sorted_bam + '.bai'

def update_mapqs(bamfilename, outfile, reference_path):
    bam = tk_bam.create_bam_infile(bamfilename)
    bam_out, _ = tk_bam.create_bam_outfile(outfile, None, None, template=bam)
    variant_heap = []
    variant_map = {}
    read_heap = []
    primary_contigs = tk_reference.load_primary_contigs(reference_path)
    for read in bam:
        tags = [(key, value) for (key, value) in dict(read.tags).iteritems()]
        tags.append(('OM',int(read.mapq)))
        read.tags = tags
        if bam.references[read.tid] not in primary_contigs:
            read.tags = [(key, value) for (key, value) in read.tags if key != 'AC' and key != 'XC']
            bam_out.write(read)
            continue
        add_variant_counts(read, variant_heap, variant_map)
        heapq.heappush(read_heap, (read.pos, read))
        update_updatable(read_heap, read.pos, variant_heap, variant_map, bam_out)
    update_updatable(read_heap, 500000000, variant_heap, variant_map, bam_out, empty_me = True)
    bam_out.close()


def add_variant_counts(read, variant_heap, variant_map):
    for variant in get_variants_best(read):
        variant_id = get_variant_id(variant, read)
        if variant_id in variant_map:
            variant_map[variant_id] += 1
        else:
            variant_map[variant_id] = 1
            heapq.heappush(variant_heap, variant_id)

def update_updatable(read_heap, pos, variant_heap, variant_map, bam_out, empty_me = False):
    while len(read_heap) > 0 and ((pos - read_heap[0][0] > 300 or empty_me) or len(read_heap) > 3000):
        _, read = heapq.heappop(read_heap)
        v_best = Set()
        v_best_full = {}
        for variant_in_best in get_variants_best(read):
            bases = bases_used(variant_in_best)
            v_best.add(bases)
            v_best_full[bases] = variant_in_best
        v_second_best = Set()
        for variant_in_second_best in get_variants_second_best(read):
            v_second_best.add(bases_used(variant_in_second_best))
        diff1 = v_best - v_second_best
        best_candidates = [v_best_full[b] for b in diff1]
        for variant_in_best in best_candidates:
            count = variant_map[get_variant_id(variant_in_best, read)]
            if count > 1:
                penalty = get_penalty(variant_in_best)
                modified_penalty = float(penalty)/float(count)
                AS = float(dict(read.tags).get('AS'))
                calculated_xs = AS - (read.mapq/10.0) # use this instead of the actual XS as this has already taken into account multiple alternate alignments
                new_mapq = (AS - penalty + modified_penalty - calculated_xs) * 10.0
                new_mapq = int(min(new_mapq, 60.0))
                read.mapq = new_mapq
        read.tags = [(key, value) for (key, value) in read.tags if key != 'AC' and key != 'XC']
        bam_out.write(read)
    while len(variant_heap) > 0 and (pos - variant_heap[0][0] > 400 or empty_me):
        variant_map.pop(variant_heap[0], None)
        heapq.heappop(variant_heap)

def get_penalty(variant):
    varlength = variant[2]
    if varlength == 1:
        return -2
    else:
        return -3

def get_variants_best(read):
    return get_variants(read, "AC")

def get_variants_second_best(read):
    return get_variants(read, "XC")

def get_variants(read, tag):
    var_string = dict(read.tags).get(tag)
    if var_string is None:
        return []
    else:
        variants = var_string.split(";")
        ret = []
        for v in variants:
            if len(v) > 1:
                info = v.split(',')
                if len(info) < 3:
                    print v
                    print tag
                ret.append((int(info[0]), int(info[1]), int(info[2]))) # format is pos, start base in read, length
    return ret

def bases_used(variant):
    return (variant[1],variant[2]) # returning start base in read, length

def get_variant_id(variant,read):
    return (variant[0], variant[2], read.seq[variant[1]: variant[1] + max(0, variant[2])]) # returning pos in reference, length, bases (empty string for deletion)
