#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Setup read chunks.
#
import os.path
import tenkit.fasta as tk_fasta
import itertools
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import tenkit.preflight as tk_preflight
import tenkit.safe_json
import tenkit.bam as tk_bam
import barcodes.utils as bc_utils
import martian
import gzip
import re

__MRO__ = """
stage SETUP_CHUNKS(
    in  string    sample_id          "id of the sample",
    in  map[]     sample_def         "list of dictionary specifying input data",
    in  string    input_mode         "configuration of the input fastqs",
    in  string    barcode_whitelist,
    in  map       downsample         "map specifies either subsample_rate (float) or gigabases (int)",
    out map[]     chunks             "map has barcode, barcode_reverse_complement, sample_index, read1, read2, gem_group, and read_group fields",
    out string[]  read_groups        "list of strings representing read groups",
    out json      downsample_info    "info about downsampling"
    out txt       barcode_whitelist_path,
    out int       requested_read_pairs,
    src py        "stages/reads/setup_chunks",
)
"""

'''
== Subsampling Documentation ==

SETUP_CHUNKS and the following stages BUCKET_FASTQS and SORT_FASTQS
work together to  implement various kinds of subsampling of reads and barcode.

Subsampling is controlled by the `downsample` map parameter to SETUP_CHUNKS.

The overall data quantity used by the pipeline can be specified by including on of the following key-value pairs:
 'gigabases': float  -- target final amount of GB of sequence data (sum of lengths of R1 and R2) to use
 'target_reads': int -- target final number of single-end reads to use
 'subsample_rate': float -- target fraction of full data sets to use

If the request in 'gigabases' or 'target_reads' is not achievable given the amount of data available and the bc subsampling rate (see below),
the pipeline will use all the available data in the sampled set of barcode.

You can instruct the pipeline to only use data from fraction of the barcodes in the experiment by adding the following key-value pair:
 'bc_subsample_rate': float -- select a pseduo-random set of barcodes. Each barcode has the given probability of be included, regardless of how many reads it contains. The pseudo-random number generator is seeded with a constant value, so repeated runs of the pipeline will always use the same subset of barcodes. Subsampling is based on the corrected barcode sequence of each read. Reads without a corrected barcode are also subsampled at the bc_subsample_rate. 


= Implementation Notes =

SETUP_CHUNKS makes a coarse estimate of the amount of input data available by scanning the beginning of each fastq.gz file and making extrapolations. 
It computes an 'first-pass' subsampling rate based on the requested parameters, and includes 20% extra as a margin of error, due to possible inaccuracies in the FASTQ size estimation.

BUCKET_FASTQS applies the barcode subsampling as given and the read sampling rate computed in SETUP_CHUNKS serially, and counts the number of reads that pass. It then computes a 'second-pass' subsampling
rate to acheive (within binomial sampling error) the requested number of reads.
'''



def main(args, outs):
    """Combine reads from multiple input FASTQ files, and potentially trim.
       Demultiplex outputs a series of FASTQ files with filenames of the form:
       read-[RA|I1|I2]_si-AGTAACGT_lane-001_chunk_001.fastq[.gz].
    """

    def check_key(n, dict_in, name, tys):
        if not dict_in.has_key(name):
            martian.exit("Entry %d in sample_def missing required field: %s" % (n, name))

        if not (type(dict_in[name]) in tys):
            martian.exit("Entry %d in sample_def for '%s' has incorrect type -- expecting %s, got %s" % (n, name, str(tys), type(dict_in[name])))


    if args.downsample is not None:
        if len(args.downsample.keys()) > 1:
            martian.exit("More than one downsampling mode requested. Please select a single downsampling mode")

        (k,v) = args.downsample.items()[0]
        if not k in ["gigabases", "subsample_rate", "target_reads"]:
            martian.exit("Unrecognized downsampling mode: %s" % k)

    # Check for self-consistent gem_group settings in the sample_def entries
    gem_groups = [x['gem_group'] for x in args.sample_def]
    all_null = all([x is None for x in gem_groups])
    all_int = all([type(x) is int for x in gem_groups])

    if not (all_null or all_int):
        martian.exit("Inconsistent gem_group tags. Please specify all gem_group tags as null, or all gem_group tags with an integer")

    # If all gem_groups are set to null, then set them all to 1
    if all_null:
        for sample_item in args.sample_def:
            sample_item['gem_group'] = 1

    # Predicted input bases
    total_seq_bases = 0
    # Predicted input reads
    total_input_reads = 0

    # verify input mode upfront
    if args.input_mode not in ["BCL_PROCESSOR", "ILMN_BCL2FASTQ"]:
        martian.throw("Unrecognized input_mode: %s" % args.input_mode)

    for (idx, sample_item) in enumerate(args.sample_def):
        # validate fields
        check_key(idx, sample_item, "read_path", [str, unicode])
        check_key(idx, sample_item, "lanes",  [list, type(None)])
        check_key(idx, sample_item, "gem_group", [int, type(None)])
        if args.input_mode == "BCL_PROCESSOR":
            check_key(idx, sample_item, "sample_indices", [list, type(None)])
        elif args.input_mode == "ILMN_BCL2FASTQ":
            check_key(idx, sample_item, "sample_names", [list, type(None)])

    interleaved_read_type = "RA"

    chunks = []
    read_groups = set()

    for read_chunk in args.sample_def:
        # Each sample_def entry can have a separate pre-applied downsampling rate
        # We adjust the estimated data in that chunk to account for this
        # subsampling
        chunk_subsample_rate = read_chunk.get('subsample_rate', 1.0)

        bc_in_read = {}
        if read_chunk.has_key('bc_in_read'):
            if read_chunk['bc_in_read'] is not None:
                bc_in_read['bc_in_read'] = read_chunk['bc_in_read']
                bc_in_read['bc_length'] = read_chunk['bc_length']

        path = read_chunk['read_path']
        lanes = read_chunk['lanes']
        gem_group = read_chunk['gem_group']
        unbarcoded = read_chunk.get('unbarcoded')
        sample_id = args.sample_id
        library_id = read_chunk.get('library_id', 'MissingLibrary')

        # split on BCL_PROCESSOR / ILMN_BCL2FASTQ
        # the main difference is that BCL_PROCESSOR uses interleaved reads and labels FASTQs by sample index;
        # whereas ILMN_BCL2FASTQ uses R1/R2 and labels by sample name

        if args.input_mode == "BCL_PROCESSOR":
            sample_index_strings, msg = tk_preflight.check_sample_indices(read_chunk)
            if sample_index_strings is None:
                martian.exit(msg)

            sample_seq_bases = 0
            read_length = 100 # Should be overwritten below

            find_func = tk_fasta.find_input_fastq_files_10x_preprocess
            for sample_index in sample_index_strings:
                # process interleaved reads
                reads = find_func(path, interleaved_read_type, sample_index, lanes)
                for read in reads:
                    _, predicted_seq_bases, read_length = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases

            sample_seq_bases = chunk_subsample_rate * sample_seq_bases
            bp_per_read_pair = 2*read_length

            martian.log_info("Input data: Predict %f GB from %s. (%d bp per read pair)" % (float(sample_seq_bases)/1e9, path, bp_per_read_pair))
            total_seq_bases += sample_seq_bases
            total_input_reads += float(sample_seq_bases)/read_length

            for sample_index in sample_index_strings:
                reads = find_func(path, interleaved_read_type, sample_index, lanes)
                # TODO confirm that this works with cellranger
                si_read, bc_read = ("I1", "I2")
                if 'barcode_read' in read_chunk and read_chunk['barcode_read'] == 'I1':
                    si_read, bc_read = ("I2", "I1")
                sis = find_func(path, si_read, sample_index, lanes)

                # allow empty sample index case if all reads in lane are same sample
                if sis is None or sis == []:
                    sis = [None] * len(reads)

                if not unbarcoded:
                    barcodes = find_func(path, bc_read, sample_index, lanes)
                    if len(barcodes) == 0:
                        barcodes = [None] * len(reads)
                else:
                    barcodes = [None] * len(reads)

                # calculate chunks
                for r,b,si in zip(reads, barcodes, sis):
                    (flowcell, lane) = get_run_data(r)
                    rg_string = tk_bam.pack_rg_string(sample_id, library_id, gem_group, flowcell, lane)
                    new_chunk = {
                        'read1': r, 'read2': None, 'reads_interleaved': True, 'barcode': b,
                        'sample_index': si, 'barcode_reverse_complement': False, 'gem_group': gem_group,
                        'subsample_rate': chunk_subsample_rate, 'read_group': rg_string
                    }
                    new_chunk.update(bc_in_read)
                    chunks.append(new_chunk)
                    read_groups.add(rg_string)

        elif args.input_mode == "ILMN_BCL2FASTQ":
            sample_names = read_chunk['sample_names']

            read_length1 = None
            read_length2 = None
            sample_seq_bases = 0
            find_func = tk_fasta.find_input_fastq_files_bcl2fastq_demult
            for sample_name in sample_names:
                # process read 1
                reads = find_func(path, "R1", sample_name, lanes)
                for read in reads:
                    _, predicted_seq_bases, read_length1 = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases
                # process read 2
                reads = find_func(path, "R2", sample_name, lanes)
                for read in reads:
                    _, predicted_seq_bases, read_length2 = fastq_data_estimate(read)
                    sample_seq_bases += predicted_seq_bases

            if read_length1 is None and read_length2 is None:
                martian.exit("No input FASTQs were found for the requested parameters.")
            elif read_length1 is None:
                martian.exit("No input FASTQs were found for Read1.")
            elif read_length2 is None:
                martian.exit("No input FASTQs were found for Read2.")

            sample_seq_bases = chunk_subsample_rate * sample_seq_bases
            bp_per_read_pair = read_length1 + read_length2

            martian.log_info("Input data: Predict %f GB from %s. (%d bp per read pair)" % (float(sample_seq_bases)/1e9, path, bp_per_read_pair))
            total_seq_bases += sample_seq_bases
            total_input_reads += float(sample_seq_bases)*2/(read_length1 +
                read_length2)

            for sample_name in sample_names:
                r1_reads = find_func(path, "R1", sample_name, lanes)
                r2_reads = find_func(path, "R2", sample_name, lanes)

                # TODO confirm that this works with cellranger
                si_read, bc_read = ("I1", "I2")
                if 'barcode_read' in read_chunk and read_chunk['barcode_read'] == 'I1':
                    si_read, bc_read = ("I2", "I1")
                sis = find_func(path, si_read, sample_name, lanes)

                # allow empty sample index case if all reads in lane are same sample
                if sis is None or sis == []:
                    sis = [None] * len(r1_reads)

                # in Chromium chemistry... there shouldn't be separate barcode reads...
                if not unbarcoded:
                    barcodes = find_func(path, bc_read, sample_name, lanes)
                    if len(barcodes) == 0:
                        barcodes = [None] * len(r1_reads)
                else:
                    barcodes = [None] * len(r1_reads)

                # again, with Chromium, the barcodes should be an array of Nones, but
                # just in case...
                if not (len(r1_reads) == len(r2_reads) == len(barcodes)):
                    martian.log_info("Read 1 files: %s" % str(r1_reads))
                    martian.log_info("Read 2 files: %s" % str(r2_reads))
                    martian.log_info("Barcode files: %s" % str(barcodes))
                    martian.exit("Read1, Read2, and Barcode files are mismatched. Exiting pipline")

                # calculate chunks
                for r1,r2,b,si in zip(r1_reads, r2_reads, barcodes, sis):
                    (flowcell, lane) = get_run_data(r1)
                    rg_string = tk_bam.pack_rg_string(sample_id, library_id, gem_group, flowcell, lane)
                    new_chunk = {
                        'read1': r1, 'read2': r2, 'reads_interleaved': False, 'barcode': b,
                        'sample_index': si, 'barcode_reverse_complement': False, 'gem_group': gem_group,
                        'subsample_rate': chunk_subsample_rate, 'read_group': rg_string
                    }
                    new_chunk.update(bc_in_read)
                    chunks.append(new_chunk)
                    read_groups.add(rg_string)

    martian.log_info("Input data: Predict %f total GB" % (float(total_seq_bases)/1e9))

    if len(chunks) == 0:
        martian.exit("No input FASTQs were found for the requested parameters.")


    #
    # Downsampling setup
    #

    # The total available input raw gigabases of input data (est_gb), and the base pairs per read pair (bp_per_read_pair)
    # are estimated above.
    (est_gb, bp_per_read_pair) = (float(total_seq_bases)/1e9, bp_per_read_pair)

    downsample = args.downsample if args.downsample is not None else {}

    # Possible BC subsampling -- try to get the requested amount of data _after_ bc subsampling
    est_gb_post_bc = est_gb * downsample.get("bc_subsample_rate", 1.0)

    # Aim high to ensure that we won't be left with too few reads
    fudge_factor = 1.00

    downsample_succeeded = True

    if downsample.has_key("gigabases"):
        read_sample_rate = min(1.0, fudge_factor * downsample['gigabases'] / est_gb_post_bc)
        requested_read_pairs = int(1e9 * downsample['gigabases'] / bp_per_read_pair)
        downsample_succeeded = downsample['gigabases'] > est_gb_post_bc

    elif downsample.has_key("target_reads"):
        requested_read_pairs = int(downsample['target_reads'] / 2)
        est_read_pair_post_bc = 1e9 * est_gb_post_bc / bp_per_read_pair
        read_sample_rate = min(1.0, fudge_factor * requested_read_pairs / est_read_pair_post_bc)
        downsample_succeeded = requested_read_pairs > est_read_pair_post_bc

    elif downsample.has_key("subsample_rate"):
        read_sample_rate = min(1.0, downsample["subsample_rate"] / downsample.get("bc_subsample_rate", 1.0))
        requested_read_pairs = None
    else:
        if len(downsample.keys()) > 0:
            martian.exit("Unrecognized downsample request: %s.\n Please use 'gigabases', 'target_reads', or 'subsample_rate'" % str(downsample))
        read_sample_rate = 1.0
        requested_read_pairs = None

    ## Alert if user requests analysis on too many reads
    ## Three CS scenarios:
    ## no downsampling
    ## "gigabases" downsampling
    ## "target_reads" downsampling
    READ_THRESHOLD = 5*1000*1000*1000
    est_reads_post_ds = (requested_read_pairs*2 if requested_read_pairs is 
        not None else total_input_reads)
    martian.log_info("Estimate %.3f M reads entering pipeline" % 
        (est_reads_post_ds/1e6))
    if est_reads_post_ds > READ_THRESHOLD:
        martian.alarm("We will be processing data from %.1f billion reads "\
            "and the pipeline run time will likely exceed 24 hours. Please "\
            "consult the 10x support website for guidance on run times. You "\
            "can reduce the number of reads using the downsample/maxreads "\
            "command-line option." % (est_reads_post_ds/1e9))

    martian.log_info("Downsampling request: %s" % str(downsample))
    martian.log_info("Base pairs per read pair: %s" % bp_per_read_pair)
    martian.log_info("Estimated Input: %.2f GB, Initial Downsample Rate: %.3f. Requested total reads: %s" % (est_gb, read_sample_rate, str(requested_read_pairs)))

    # Copy over the per-chunk subsample rates
    if read_sample_rate is not None:
        for chunk in chunks:
            chunk['subsample_rate'] = chunk.get('subsample_rate', 1.0) * read_sample_rate
            if downsample.has_key("bc_subsample_rate"):
                chunk["bc_subsample_rate"] = downsample["bc_subsample_rate"]

    outs.requested_read_pairs = requested_read_pairs

    martian.log_info("Input reads: %s" % str(chunks))
    outs.chunks = chunks
    outs.read_groups = [rg for rg in read_groups]

    downsample_info = {}
    downsample_info['available_gb'] = est_gb
    downsample_info['requested_gb'] = downsample.get('gigabases', None)
    downsample_info['requested_rate'] = read_sample_rate
    downsample_info['post_downsample_gb'] = float(requested_read_pairs * bp_per_read_pair) / 1e9 if requested_read_pairs is not None else None
    downsample_info['downsample_succeeded'] = downsample_succeeded

    with open(outs.downsample_info, 'w') as downsample_out:
        tenkit.safe_json.dump_numpy(downsample_info, downsample_out)

    check_fastqs(outs.chunks)

    # Give out full path to BC whitelist
    if args.barcode_whitelist:
        outs.barcode_whitelist_path = bc_utils.barcode_whitelist_path(args.barcode_whitelist)
    else:
        outs.barcode_whitelist_path = None


def check_fastqs(chunks):
    keys = ['read1', 'read2', 'barcode', 'sample_index']

    def check_fastq(fastq):
        # Check if fastq is readable
        if not os.access(fastq, os.R_OK):
            martian.exit("Do not have file read permission for FASTQ file: %s" % fastq)

        # Check if fastq is gzipped
        gzip_suffix = '.gz'
        is_gzip_fastq = True
        try:
            with gzip.open(fastq) as f:
                f.read(1)
        except:
            is_gzip_fastq = False

        if is_gzip_fastq and not fastq.endswith(gzip_suffix):
            martian.exit("Input FASTQ file is gzipped but filename does not have %s suffix: %s" % (fastq, gzip_suffix))
        if not is_gzip_fastq and fastq.endswith(gzip_suffix):
            martian.exit("Input FASTQ file is not gzipped but filename has %s suffix: %s" % (fastq, gzip_suffix))

    for chunk in chunks:
        for key in keys:
            fastq = chunk.get(key)
            if fastq is not None:
                check_fastq(fastq)

def infer_barcode_reverse_complement(barcode_whitelist, barcode_files):
    rc_valid_count = 0
    reg_valid_count = 0
    if barcode_whitelist:
        barcode_rc = []
        for barcode_file in barcode_files:
            read_num = 0

            if barcode_file[-3:] == ".gz":
                barcode_open_file = gzip.open(barcode_file)
            else:
                barcode_open_file = open(barcode_file, 'r')
            read_iter = tk_fasta.read_generator_fastq(barcode_open_file)
            for (name, seq, qual) in read_iter:
                if seq in barcode_whitelist:
                    reg_valid_count += 1
                if tk_seq.get_rev_comp(seq) in barcode_whitelist:
                    rc_valid_count += 1
                if read_num > 1000:
                    break
                read_num += 1

            if tk_stats.robust_divide(float(rc_valid_count), float(rc_valid_count + reg_valid_count)) > 0.75:
                barcode_rc.append(True)
            else:
                barcode_rc.append(False)
            barcode_open_file.close()
        return barcode_rc
    else:
        return [False] * len(barcode_files)

def fastq_data_estimate(fn, num_reads = 4000000):
    # Open reader
    if fn[-2:] == 'gz':
        reader = gzip.open(fn)
        is_gz = True
    else:
        reader = open(fn, 'r')
        is_gz = False

    gen = tk_fasta.read_generator_fastq(reader)
    rds = itertools.islice(gen, num_reads)

    input_lens = [(len(header) + len(r) + len(qual) + 5, len(r)) for (header, r, qual) in rds]
    total_seq_len = sum(x[1] for x in input_lens)
    total_data_len = sum(x[0] for x in input_lens)
    file_sz = os.path.getsize(fn)
    read_length = total_seq_len / len(input_lens)

    # NOTE: do not try and use the gzip footer containing the length of the compressed data
    # that only reflects the length of the final gzip block. A valid gzip file may have
    # many blocks, so that field cannot be relied upon.

    if is_gz:
        compressed_sz = reader.myfileobj.tell()
        predicted_sz = total_data_len / compressed_sz * file_sz
    else:
        predicted_sz = file_sz

    read_yield = float(len(input_lens)) / total_data_len
    seq_yield = float(total_seq_len) / total_data_len
    predicted_reads = read_yield * predicted_sz
    predicted_seq = seq_yield * predicted_sz

    # Log estimate of downsampling
    martian.log_info("Estimates for: %s" % fn)
    dbg_str =  "compressed_size: %.2f, predicted_size: %.2f" % \
               (file_sz / 1e9, predicted_sz / 1e9)
    martian.log_info(dbg_str)

    return (predicted_reads, predicted_seq, read_length)


def get_run_data(fn):
    """ Parse flowcell + lane from the first FASTQ record.
    NOTE: we don't check whether there are multiple FC / lanes in this file.
    """
    if fn[-2:] == 'gz':
        reader = gzip.open(fn)
    else:
        reader = open(fn, 'r')

    gen = tk_fasta.read_generator_fastq(reader)

    try:
        (name, seq, qual) = gen.next()
        read_name_split = re.split(':', name)
        if len(read_name_split) >= 4:
            (flowcell, lane) = read_name_split[2:4]
        else:
            (flowcell, lane) = ("unknown_fc", 0)
        return (flowcell, lane)
    except StopIteration:
        # empty fastq
        martian.exit("FASTQ is empty: %s" % fn)
