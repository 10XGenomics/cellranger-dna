#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import os
import re
import socket

import martian
import tenkit.fasta as tk_fasta
import tenkit.preflight as tk_preflight
from longranger.cnv import contig_manager

__MRO__ = """
stage CNV_CALLER_PREFLIGHT(
    in string sample_id,
    in string fastq_mode,
    in map[]  sample_def,
    in int    force_cells,
    in map    cnv_params,
    in map    downsample,
    in string reference_path,
    in bool   check_executables,
    src py     "stages/preflight/cnv_caller",
)
"""

def main(args, outs):
    hostname = socket.gethostname()

    # Sample ID / pipestance name
    if args.sample_id is not None:
        if not re.match("^[\w-]+$", args.sample_id):
            martian.exit("Sample name may only contain letters, numbers, underscores, and dashes: " + args.sample_id)

    # Check numerical options
    # types are already checked by mrp so only need to check ranges
    if args.force_cells is not None and (args.force_cells < 1 or
        args.force_cells > 20000):
        martian.exit("MRO parameter force_cells must be a positive integer"\
            " <= 20000.")

    # check min_ploidy, max_ploidy
    if args.cnv_params is not None:
        min_ploidy = args.cnv_params.get("min_ploidy", None)
        max_ploidy = args.cnv_params.get("max_ploidy", None)
        if min_ploidy is not None and min_ploidy <= 0:
            martian.exit("Command line argument soft-min-avg-ploidy must be a "\
                "positive real number.")
        if max_ploidy is not None and (max_ploidy <= 0 or max_ploidy > 8.0):
            martian.exit("Command line argument soft-max-avg-ploidy must be a "\
                "positive real number <= 8.")
        if (min_ploidy is not None and max_ploidy is not None and 
            max_ploidy <= min_ploidy):
            martian.exit("Command line arguments must satisfy "\
                "soft-min-avg-ploidy < soft-max-avg-ploidy.")

    # check downsample options
    if args.downsample is not None and len(args.downsample.keys()) > 0:
        keys = args.downsample.keys()
        if len(keys) > 1:
            martian.exit("Please supply either maxreads or downsample but not "\
                "both.")
        key = keys[0]
        value = args.downsample[key]
        param_map = {"target_reads" : "maxreads", "gigabases" : "downsample"}
        bad_value = False
        try:
            float(value)
            bad_value = value < 1e-12
        except ValueError:
            bad_value = True
        if bad_value:
            cs_key = param_map[key]
            martian.exit("Command line argument %s must be a positive number" %
                cs_key)

    # FASTQ input
    for idx, sample_def in enumerate(args.sample_def):
        read_path = sample_def["read_path"]
        if not read_path:
            martian.exit("Must specify a read_path containing FASTQs.")
        if not read_path.startswith('/'):
            martian.exit("Specified FASTQ folder must be an absolute path: %s" % read_path)
        if not os.path.exists(read_path):
            martian.exit("On machine: %s, specified FASTQ folder does not exist: %s" % (hostname, read_path))
        if not os.access(read_path, os.X_OK):
            martian.exit("On machine: %s, longranger does not have permission to open FASTQ folder: %s" % (hostname, read_path))
        if not os.listdir(read_path):
            martian.exit("Specified FASTQ folder is empty: " + read_path)

        library_id = sample_def.get("library_id")
        if library_id is not None:
            if not re.match("^[\w-]+$", library_id):
                martian.exit("Library name may only contain letters, numbers, underscores, and dashes: " + library_id)

        lanes = sample_def["lanes"]
        if lanes is not None:
            for lane in lanes:
                if not tk_preflight.is_int(lane):
                    martian.exit("Lanes must be a comma-separated list of numbers.")

        if args.fastq_mode == "BCL_PROCESSOR":
            sample_indices, msg = tk_preflight.check_sample_indices(sample_def)
            if sample_indices is None:
                martian.exit(msg)

            find_func = tk_fasta.find_input_fastq_files_10x_preprocess
            reads = []
            for sample_index in sample_indices:
                # process interleaved reads
                reads.extend(find_func(read_path, "RA", sample_index, lanes))
            if len(reads) == 0:
                martian.exit("No input FASTQs were found for the requested parameters.")
        elif args.fastq_mode == "ILMN_BCL2FASTQ":
            sample_names = sample_def.get("sample_names", None)
            if sample_names is None:
                martian.exit("Entry {} in sample_def missing required field: sample_names".format(idx))
            find_func = tk_fasta.find_input_fastq_files_bcl2fastq_demult
            reads1 = []
            reads2 = []
            for sample_name in sample_names:
                r1 = find_func(read_path, "R1", sample_name, lanes)
                r2 = find_func(read_path, "R2", sample_name, lanes)
                if len(r1) != len(r2):
                    martian.exit("Entry {} in sample_defs are missing input FASTQs.".format(idx))
                reads1.extend(r1)
                reads2.extend(r2)
            if len(reads1) == 0 and len(reads2) == 0:
                martian.exit("No input FASTQs were found for the requested parameters.")
        else:
            martian.exit("Unrecognized fastq_mode: {}".format(args.fastq_mode))

    # Reference
    ok, msg = tk_preflight.check_refdata(args.reference_path, max_contigs=None)
    if ok:
        martian.log_info(msg)
    else:
        martian.exit(msg)
    contig_defs_json_path = os.path.join(args.reference_path, "fasta", 
        "contig-defs.json")
    faidx_path = os.path.join(args.reference_path, "fasta", 
        "genome.fa.fai")
    error_msg = contig_manager.verify_contig_defs(contig_defs_json_path,
        faidx_path)
    if error_msg is not None:
        martian.exit(error_msg)

    try:
        ref = contig_manager.contig_manager(args.reference_path)
    except Exception as e:
        martian.exit("Unexpected error occurred.\n%s"%str(e))

    # too many contigs
    primary = ref.primary_contigs(allow_sex_chromosomes=True)
    num_primary_contigs = len(primary)
    if num_primary_contigs > 100:
        martian.exit("There can be at most 100 primary contigs.")

    # contig length checks
    chrom_length_dict = ref.get_contig_lengths()

    contig_length_exit = 500 * 1000
    contig_length_warn = 10 ** 7
    offending_contigs_warn = []
    offending_contigs_exit = []
    for c in primary:
        clen = chrom_length_dict[c]
        if clen < contig_length_exit:
            offending_contigs_exit.append(c)
        elif clen < contig_length_warn:
            offending_contigs_warn.append(c)
    if len(offending_contigs_exit) > 0:
        martian.exit("Primary contig(s) \"%s\" are shorter than %d bases. "\
            "Every primary contig must be at least %d bases "\
            "in length."%(",".join(offending_contigs_exit), contig_length_exit,
                          contig_length_exit))
    elif (not args.check_executables) and len(offending_contigs_warn) > 0:
        martian.alarm("Primary contig(s) \"%s\" are shorter than %d bases. "\
            "Every primary contig is recommended to be at least %d bases "\
            "in length."%(",".join(offending_contigs_warn), contig_length_warn,
                          contig_length_warn))

    # Open file handles limit 
    if args.check_executables:
        ok, msg = tk_preflight.check_open_fh()
        if not ok:
            martian.exit(msg)

    martian.log_info(tk_preflight.record_package_versions())

