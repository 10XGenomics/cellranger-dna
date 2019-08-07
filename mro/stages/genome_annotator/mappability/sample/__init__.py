#   
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

# Sample from a fasta reference to create a fasta file suitable for alignment. 
# We sample using emperical read lengths, insert sizes, and QV scores. 


import martian
import numpy as np
import json

import tenkit.bio_io as tk_io
import tenkit.fasta as fasta
import tenkit.reference as reference
import tenkit.seq as tk_seq
import longranger.cnv.contig_manager as contig_manager
from tenkit.chunk_utils import generate_chrom_loci, pack_loci
from tenkit.constants import ILLUMINA_QUAL_OFFSET

__MRO__="""
stage SAMPLE_REFERENCE(
    in  string   reference_path,
    in  json     basic_stats,
    in  json     insert_sizes,
    in  int      window_size,
    in  float    target_coverage,
    out map[]    chunks,
    out fasta    tmp,
    out int      samples_per_bin,
    src py       "stages/genome_annotator/mappability/sample",
) split using(
    in  string[] loci,
)"""

# Split on slices of the reference
def split(args):
    ref = contig_manager.contig_manager(args.reference_path)
    contig_lengths = ref.get_contig_lengths( )

    target_regions = None
    all_loci = []
    for (chrom_name, chrom_size) in contig_lengths.iteritems():
        all_loci.extend(generate_chrom_loci(target_regions, chrom_name, chrom_size,
            100000000))
        
    locus_sets = pack_loci(all_loci)

    chunk_defs = [{'loci': loci} for loci in locus_sets]
    return {'chunks': chunk_defs}

def main(args, outs):
    """For each slice produce a fasta file sampling reads from that slice. 
    We split our section of the genome into a bunch of 20kb chunks. For each
    chunk we sample an identical number of paired end reads. The name of each
    read encodes the true position that it was sampled from."""

    # Grab basic stats for the read lengths and quality scores
    stats_fp = open(args.basic_stats)
    stats = json.load(stats_fp)
    
    # Fix the random seed
    np.random.seed(0)

    # Info is a map we use everywhere to track the sampling parameters. 
    # r1_len: the length of read1
    # r2_len: the length of read2
    # insert_size_map: a map of insert-size (as a string) to frequency
    # q_score_map a map of quality score (as a string) to frequency


    info={'r1_len':stats['r1_len'],
          'r2_len':stats['r2_len']
         }
    

    info['q_score_map'] = {
            '30': stats['bc_q30_bases'],
            '20': stats['bc_q20_bases'] - stats['bc_q30_bases'],
            '0':  stats['bc_tot_bases'] - stats['bc_q20_bases']}
     
    stats_is_fp = open(args.insert_sizes)
    info['insert_size_map'] = json.load(stats_is_fp)['60']

    # How many samples will we make from each window?
    samples = int(round(2.0*args.target_coverage*(float(args.window_size) / (stats['r1_len']+stats['r2_len']))))

    martian.log_info("Using %i samples per %i bin" %(samples, args.window_size))
    output_path= martian.make_path("chnk.fasta")
    output=open(output_path, "w")

    ref = reference.open_reference(args.reference_path)
    #Loop over every window in every loci.
    for (chrom, start, end) in (tk_io.get_locus_info(l) for l in args.loci):
        cur = start
        while (cur < end):
            # Sample |samples| reads from chrom:cur-chrom:cur+window_size and put
            # the results in the output file
            perbin(chrom, cur, ref, output, info, args.window_size, samples)
            cur+=args.window_size
    outs.tmp=output_path
    outs.samples_per_bin = samples
    output.close()

def join(args, outs, chunk_defs, chunk_outs):
    """Join doesn't actual do any work. Its output is a array of chunks suitable for
    injestion by a subsequent ALIGN stage."""
    a = [{'barcode': None,
          'barcode_reverse_complement': False,
          'gem_group': 1,
          'read1': x.tmp,
          'read2': None,
          'read_group': "test:test:test:test:test",
          'reads_interleaved': True,
          'sample_index': None} for x in chunk_outs]
    outs.chunks=a
    outs.samples_per_bin=chunk_outs[0].samples_per_bin

# Compute a normalized, weighted histogram array from a map of frequencies
def weighted_histogram_from_map(histogram_map):
    keys = [(x) for x in histogram_map.keys()]
    a = np.zeros(len(keys))
    for i in range(len(keys)):
        a[i] = histogram_map[keys[i]]

    a/=a.sum()

    return (keys,a)

# Given a QV score, return the expected error roate. 
# TODO: This function should be adjusted to reflect the fact that QV scores
# systematically underestimate the error rate.
def error_rate_from_qv(qv, info):
    return 10.0**(-qv/10.0)

# Perform work for one bin. This makes |samples| synthetic reads from the |ref|
# from |chrom|:|bin_start| to |chrom|:|bin_start|+window size. It adds
# errors to each sampled read and puts the result in the output file
def perbin(chrom, bin_start, ref, output, info, window_size, samples):
    
    ref[chrom].as_string=False

    q_scores, q_probs=weighted_histogram_from_map(info['q_score_map'])

    i_length, i_probs=weighted_histogram_from_map(info['insert_size_map'])

    r1_len = info['r1_len']
    r2_len = info['r2_len']
    for _ in range(0,samples):
        # Decide where R1 and R2 should start
        start1 = bin_start + np.random.randint(window_size)
        start2 = start1 + int(np.random.choice(i_length, p=i_probs))
        
        if (start1+r1_len >= len(ref[chrom]) or start2+r2_len >= len(ref[chrom])):
            continue
        # Grab the sequence for R1 and R2
        seq1 = ref[chrom][start1:start1+r1_len].copy()
        seq2 = [c for c in tk_seq.get_rev_comp("".join(ref[chrom][start2:start2+r2_len]))]
        
        # Make up quality scores
        qscore1 = np.random.choice(q_scores, p=q_probs, size=r1_len)
        qscore2 = np.random.choice(q_scores, p=q_probs, size=r2_len)

        # Insert errors in each sequence
        #insert_errors(seq1, qscore1, info)
        #insert_errors(seq2, qscore2, info)

        # Name the sequence. The format of the name is important, ANALYZE_OUTPUT
        # parses this to figure out where a sequence SHOULD have mapped.
        name = "X:%s:%i:%i" % (chrom, start1, start2)

        qscore1_ascii = "".join([qscore_map(x) for x in qscore1])
        qscore2_ascii = "".join([qscore_map(x) for x in qscore2])

        fasta.write_read_pair_fastq(output, name, "".join(seq1), qscore1_ascii, name, "".join(seq2), qscore2_ascii)


# Return the ascii PHRED score for a QV.
def qscore_map(x):
    return(chr(int(x)+ILLUMINA_QUAL_OFFSET))


# Insert errors into a sequence, based on an array of quality scores
def insert_errors(sequence, qscore, info):
    for i in range(0,len(sequence)):
        p = error_rate_from_qv(int(qscore[i]), info)
        if (sequence[i] != 'N'):
            if np.random.rand() < p:
                # XXX This isn't actually right. There is a 25% change that we
                # change sequence[i] back to its self. However I don't think this
                # matters very much.
                sequence[i]=np.random.choice(['A','C','G','T'])


#def rc(seq):
#    """Reverse complement a sequence"""
#    return seq
#    rcmap={
#            'a': 't',
#            'A': 'T',
#            'c': 'g',
#            'C': 'G',
#            'g': 'c',
#            'G': 'g',
#            't': 'a',
#            'T': 'A',
#            'n': 'n',
#            'N': 'N'
#    }
#
#    ns = ([rcmap[x] for x in seq])
#    ns.reverse()
#    return ns

