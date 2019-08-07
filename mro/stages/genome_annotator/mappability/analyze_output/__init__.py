#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
from longranger.cnv import coverage_matrix
import numpy as np
import tenkit.bam as tk_bam
import martian
import pandas as pd
import crdna.read_filter

__MRO__="""
stage ANALYZE_OUTPUT(
    in bam[] bams,
    in string reference_path,
    in int samples_per_bin,
    in int window_size,
    out h5 mappability,
    src py "stages/genome_annotator/mappability/analyze_output",
) split using (
    in bam bam,
)"""

# Split: input is an array of bam paths from ALIGN. Just run each chunk
# on that array.
def split(args):

    chunk_defs = [ {'bam':x} for x in args.bams]

    return {'chunks':chunk_defs, 'join':{'__mem_gb':2}}


# Common preparation function to main,join. Open the input bam file
# and create mappability arrays for each chromosome.
def prepare(bam_path, window_size):

    bam_in = tk_bam.create_bam_infile(bam_path)

    chroms = bam_in.references

    chrom_sizes = dict(zip(bam_in.references, bam_in.lengths))

    chrom_store = {}
    martian.log_info(chroms)
    martian.log_info(chrom_sizes)

    for k in chroms: 
        chrom_store[k]=np.zeros(int(np.ceil(1.*chrom_sizes[k]/window_size)))
    
    return bam_in, chrom_store

# Compute the partial mappability track for a single BAM file
def main(args, outs):

    bam_in, chrom_store = prepare(args.bam, args.window_size)

    # Iterate of the BAM file. Add a point for every read that landed within
    # its window. bump_success_count scales these so that we'll end up with
    # 1.0 if a window is perfectly mapped.
    for read in bam_in.fetch(until_eof=True):
        if mapped_right(read, args.window_size):
            bump_success_count(chrom_store, read, args.window_size, args.samples_per_bin)
    
    # Dump the results to an HDF5 matrix
    store = pd.HDFStore(outs.mappability, "w")
    for c in chrom_store.keys():
        #XXX Contigs isn't really the right name here (we change it to map, later)
        # but it means the coverage_matrix API works on the H5 file.
        store["/contigs/"+c]=pd.Series(chrom_store[c])
    store.close()

def join(args, outs, chunk_defs, chunk_outs):
    """Combine a bunch of partial mappability vectors. This is easy:
    they're already normalized so that we can just add them up"""

    bam_in, chrom_store = prepare(args.bams[0], args.window_size)

    store = pd.HDFStore(outs.mappability, "w")

    # Loop over every contig of every chunk
    for chunk in chunk_outs:
        in_store = pd.HDFStore(chunk.mappability, "r")
        for c in coverage_matrix.list_all_contigs(in_store):
            chrom_store[c] += in_store[coverage_matrix.contig_to_coverage_path(c)].values
        in_store.close()

    # Format the results
    print "# mappability by bin"
    print "contig, nbins, mean, min, p25, p50, p50, max"
    for c in chrom_store.keys():
        # because of fuzzy accounting when paired reads split bins we can have
        # mappability scores > 1.0
        # it's innocent to fix that here
        np.clip(chrom_store[c],0.0,1.0,out=chrom_store[c])
        store["/map/"+c]=pd.Series(chrom_store[c])
        l = len(chrom_store[c])
        u = np.nanmean(chrom_store[c])
        p = np.nanpercentile(chrom_store[c], [0, 25, 50, 75, 100])
        print ("{}, {}, " + ", ".join(["{:.3f}"]*6)).format(c, l, u, *p)
    #
    # aggregate tracks will combine this with the constants section from other .h5 files
    #
    store["constants"] = pd.Series({"mappability_samples_per_bin": args.samples_per_bin})
    store.close()


def mapped_right(read, window_size):
    """Did a read map to the right place?"""
    # Parse the query name
    _, real_chr, real_r1_pos, real_r2_pos = read.query_name.split(":")
    
    score = crdna.read_filter.inserts_per_alignment(read)
    if (score > 0.0):

        if (read.is_read1):
            real_pos = real_r1_pos
        else:
            real_pos = real_r2_pos

        # If you're not even on the right contig....
        if (real_chr != read.reference_name):
            return False
       

        # TODO BAM file are "1" biased and fasta indexing is "0" biased. Need to account
        # for that.
        read_bin_idx = int(int(real_pos)/window_size)
        emperical_bin_idx = int(read.reference_start/window_size)

        # TODO Need to consider R2. 

        if read_bin_idx == emperical_bin_idx:
            return True
        else:
            return False
    else:
        return False

def bump_success_count(chrom_store, read, window_size, samples):
    c = chrom_store[read.reference_name]
    c[read.reference_start/window_size]+=0.5/float(samples)
    #martian.log_info("B: %s %i" % (read.reference_name, read.reference_start/window_size))

    
