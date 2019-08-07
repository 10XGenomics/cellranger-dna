#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import itertools
import numpy as np
import scipy.special
import math
import longranger.cnv.contig_manager as contig_manager
# from crdna.read_filter import stringent_read_filter
import tenkit.bam as tk_bam
import crdna.bio_io as crdna_io
import tenkit.safe_json
import tenkit.stats as tk_stats
import pandas as pd
import json
from collections import defaultdict
from tools.io_util import combine_csv
import barcodes.utils as bc_utils
from crdna.constants import RAW_BARCODE_TAG, PROCESSED_BARCODE_TAG, PROFILE_MAPQ_THRESHOLD,\
    DEFAULT_AMPLICON_LENGTH, DUPLICATE_COUNT_TAG

__MRO__ = """
stage REPORT_SINGLECELL(
    in  bam      input,
    in  string   barcode_whitelist,
    in  string   reference_path,
    in  map      cell_barcodes,
    in  h5       profiles,
    in  json     duplicate_summary,
    out json     summary,
    out csv      barnyard,
    out csv      barnyard_hits,
    out json     barcode_histogram,
    src py       "stages/reporter/report_singlecell",
) split using (
    in  string   chunk_start,
    in  string   chunk_end,
)
"""

############################### SPLIT #######################################

def split(args):
    if args.input is None or args.barcode_whitelist is None:
        chunk_defs = [{'chunk_start':"0", 'chunk_end':"0"}]
        return {'chunks': chunk_defs}

    # Some R&D bc sets have very small diversity -- don't run on them
    barcode_whitelist = bc_utils.load_barcode_whitelist(args.barcode_whitelist)
    if len(barcode_whitelist) < 100:
        chunk_defs = [{'chunk_start':"0", 'chunk_end':"0"}]
        return {'chunks': chunk_defs}

    min_chunks = 20
    if len(barcode_whitelist) > 1e6:
        min_chunks = 100

    bam_in = tk_bam.create_bam_infile(args.input)
    chunks = tk_bam.chunk_bam_records(bam_in, groupbybarcode, 
                chunk_size_gb=8.0, min_chunks=min_chunks)
    for c in chunks:
        c['__mem_gb'] = 3

    return {'chunks': chunks, 'join': {'__mem_gb': 6}}

############################### CHUNK #######################################
       
def main(args, outs):
    #min_insert_size = 0
    #max_insert_size = 1e4
    
    ## sc purity threshold: what fraction of contamination by another species
    ## will we tolerate
    SC_PURITY_THRESHOLD = 0.95

    args.coerce_strings()
    outs.coerce_strings()

    # Bail out if there no valid barcodes
    if args.barcode_whitelist is None or args.input is None:
        outs.summary = None
        return

    ## group bam records by barcode NO_BARCODE/raw barcode tag/processed barcode tag
    bam_in = tk_bam.create_bam_infile(args.input)
    bam_chunk = tk_bam.read_bam_chunk(bam_in, (args.chunk_start, args.chunk_end))
    bc_read_iter = itertools.groupby(bam_chunk, groupbybarcode)

    ## compute species_list
    refs = bam_in.references
    ref = contig_manager.contig_manager(args.reference_path)
    species_list = ref.list_species()
    has_species_info = (species_list != [""])
    species_list.sort()
    genome_size = sum(ref.get_contig_lengths().values())

    ## index cells of each species
    cell_index = {}
    for sp in species_list:
        bc_list = args.cell_barcodes.get(sp, {}).keys()
        bc_list.sort( )
        for i, b in enumerate(bc_list):
            y = cell_index.get(b, "")
            if len(y) == 0:
                cell_index[b] = "%s_cell_%d"%(sp, i)
            else:
                cell_index[b] = y + "_" + "%s_cell_%d"%(sp, i)

    ## construct and write header for barnyard file
    barnyard_file = open(outs.barnyard, 'w')
    barnyard_header = (['BC'] + ["cell_id"] +
        [s+("_" if has_species_info else "")+"reads_mapq_60" for s in species_list] +
        [s+("_" if has_species_info else "")+"contigs" for s in species_list] +
        ['mapped',
        'num_mapped_bases',
        'soft_clip_frac',
        'insert_p50',
        'num_mapped_pos',
        'mapped_frac',
        'amp_rate',
        'library_complexity',
        'dup_ratio',
        'num_pairs'] +
        ["is_%s_cell_barcode"%s for s in species_list])
    waste_keys = ["no_barcode", "non_cell_barcode", "unmapped",
                  "low_mapq_lt_%d"%PROFILE_MAPQ_THRESHOLD,
                  "dups", "denominator", "unusable_read"]
    fractional_waste_keys = [
                  "no_barcode_frac", "non_cell_barcode_frac", "unmapped_frac",
                  "low_mapq_lt_%d_frac"%PROFILE_MAPQ_THRESHOLD, "dups_frac"]

    barnyard_header.extend(waste_keys)
    barnyard_header.extend(fractional_waste_keys)
    barnyard_file.write( ",".join(barnyard_header) + "\n" )

    ## wasted data categories

    ## construct and write header for barnyard_hits file
    barnyard_hits_file = open( outs.barnyard_hits, "w" )
    bh_header = ["barcode", "is_whitelisted"]
    bh_header.extend(["is_%s_cell_barcode"%s for s in species_list])
    bh_header.extend([refname for refname in bam_in.references])
    barnyard_hits_file.write( ",".join(bh_header) + "\n" )

    # For each barocode, count # per each contig, number per each window (for each window size)
    # number per species (if available in contig), number per species
    # TODO: Add detailed matrix by contigs, windows output
    num_sc_bcs = 0
    num_qual_reads = 0
    num_sc_reads = 0

    ploidy = 2
    bc_hist = {}

    ## count number of raw barcodes that exactly match whitelist
    ## without any error correction
    raw_bc_on_whitelist = 0
    # dup_summary = json.load(open(args.duplicate_summary))
    # pcr_dup_fraction = dup_summary['dup_fraction']['pcr']
    #barcode_whitelist = bc_utils.load_barcode_whitelist(args.barcode_whitelist)
    for bc, reads in bc_read_iter:
        ## collect various forms of wasted data here per barcode
        wastebin = defaultdict(int)

        bh_hits = [0 for _ in bam_in.references]
        dup_count = 1
        non_dup = 1
        bc_count = 0
        num_per_species = defaultdict(int)
        contigs_per_species = defaultdict(set)

        total_reads_by_clip = np.zeros(2, dtype=float)

        insert_length = []
        num_pairs = 0
        num_mapped = 0
        num_mapped_bases = 0
        pos_set = set([])
        for r in reads:
            ## secondary/supplementary are never counted towards anything
            if r.is_secondary or r.is_supplementary:
                continue

            ## include everything in the denominator
            wastebin["denominator"] += 1

            ## how many reads have >= 10 soft clipped bases
            if r.cigartuples is not None:
                cigar_dict = dict(r.cigartuples)
                soft_clip_index = int(cigar_dict.get(4, 0) >= 10)
                total_reads_by_clip[soft_clip_index] += 1

            if barnyard_hits_include(r):
                bh_hits[r.tid] += 1
            ## non-whitelisted barcodes count as wasted data
            if not "-" in bc:
                wastebin["no_barcode"] += 1
                continue

            if bc[:-2] == r.get_tag(RAW_BARCODE_TAG):
                raw_bc_on_whitelist += 1

            is_cell_bc_read = True

            ## waste hierarchy
            ## if not a cell or if read doesn't belong to species, then waste
            ## else if not mapped, then waste
            ## else if mapq< 30, then waste
            ## else if dup, then waste

            ## is this is a contaminant read from a different species
            ## it is wasted
            contig = refs[r.tid]
            read_species = ref.species_from_contig(contig)
            if ( not(read_species in args.cell_barcodes) or
                 not(bc in args.cell_barcodes[read_species]) ):
                wastebin["non_cell_barcode"] += 1
                is_cell_bc_read = False
            elif r.is_unmapped:
                wastebin["unmapped"] += 1
            elif r.mapq < PROFILE_MAPQ_THRESHOLD:
                wastebin["low_mapq_lt_%d"%PROFILE_MAPQ_THRESHOLD] += 1
            elif r.is_duplicate:
                wastebin["dups"] += 1
            bad_map_or_dup = (r.is_unmapped or
                              (r.mapq < PROFILE_MAPQ_THRESHOLD) or
                              r.is_duplicate)

            if is_cell_bc_read:
                bc_count += 1
                # if (stringent_read_filter(r, True) and
                #         not(r.is_unmapped) and not(r.mate_is_unmapped)):
                #     if r.is_duplicate:
                #         dup_count += 1
                #     else:
                #         non_dup += 1
                if r.has_tag(DUPLICATE_COUNT_TAG):
                    dup_count += r.get_tag(DUPLICATE_COUNT_TAG)
                    non_dup += 1
            elif bad_map_or_dup:
                # unusable reads are those that are non-cell barcodes that are
                # also any of unmapped, low mapq, nor dups
                wastebin['unusable_read'] += 1

            ## whether we have a cell barcode or not, count these stats
            if not bad_map_or_dup:
                num_mapped += 1
                num_mapped_bases += r.reference_length

                pos_set.add((r.reference_name, r.reference_start/1000))

                ## if read is part of a proper pair, only count read or its pair
                if r.is_proper_pair:
                    if r.is_read1:
                        insert_length.append( r.template_length )
                        num_pairs += 1
                    else:
                        continue

                ## Use MAPQ >= 60 to get accurate mappings only for barnyard stuff
                if r.mapq < 60:
                    continue
                num_qual_reads += 1
                if has_species_info:
                    num_per_species[read_species] += 1
                    contigs_per_species[read_species].add(contig)
            ## end of loop over reads in this barcode
            assert wastebin['denominator'] - wastebin['no_barcode'] - wastebin['unusable_read'] == num_mapped + \
                   wastebin["low_mapq_lt_%d" % PROFILE_MAPQ_THRESHOLD] + wastebin['unmapped'] + wastebin['dups']

        ## compute the library complexity and amp
        ## NOTE: insert length is hardcoded as 250, so the amp rate is really the
        ## library complexity in different units 
        num_amplicons = num_mapped - num_pairs
        dup_ratio = tk_stats.robust_divide(float(dup_count + non_dup), float(non_dup))
        
        library_complexity = tk_stats.robust_divide(num_amplicons, (dup_ratio-1.0)*2)

        amp_rate = tk_stats.robust_divide(float(library_complexity * DEFAULT_AMPLICON_LENGTH) ,
            float(ploidy * genome_size))

        bc_hist[bc] = bc_count
        map_rate = tk_stats.robust_divide(float(num_mapped), wastebin["denominator"])
        
        ## write row to barnyard_hits file
        bh_row = [ bc, int("-" in bc)]
        for s in species_list:
            bh_row.append( int(s in args.cell_barcodes and bc in args.cell_barcodes[s]) )
        bh_row.extend( bh_hits )
        barnyard_hits_file.write(",".join(map(str, bh_row)) + "\n" )

        ## write row to barnyard file
        barnyard_row = ([bc, cell_index.get(bc, "None")] +
            [num_per_species[s] for s in species_list] +
            [len(contigs_per_species[s]) for s in species_list] +
            [num_mapped, num_mapped_bases] +
            [tk_stats.robust_divide(total_reads_by_clip[1], sum(total_reads_by_clip)),
            np.median(insert_length) if len(insert_length) else np.nan,
            len(pos_set),
            map_rate,
            amp_rate,
            library_complexity,
            dup_ratio,
            num_pairs])
        for speci in species_list:
            barnyard_row.append( int((speci in args.cell_barcodes) and 
                (bc in args.cell_barcodes[speci])) )

        for key in waste_keys:
            fkey = key + "_frac"
            if (fkey in fractional_waste_keys):
                wastebin[fkey] = tk_stats.robust_divide(float(wastebin[key]), float(wastebin["denominator"]))
        barnyard_row.extend( [ wastebin[x] for x in waste_keys ] )
        barnyard_row.extend( [ wastebin[x] for x in fractional_waste_keys ] )

        barnyard_file.write( ",".join(map(str, barnyard_row)) + "\n")
        
        ## metrics relating to purity - only for multi species
        if has_species_info and len(species_list) >= 2:
            counts_by_species = [float(num_per_species[s]) for s in species_list]

            major_species_index = np.argmax( counts_by_species )
            major_species = species_list[major_species_index]
            species_purity = tk_stats.robust_divide( counts_by_species[major_species_index],
                np.sum(counts_by_species) )

            if species_purity >= SC_PURITY_THRESHOLD:
                num_sc_bcs += 1
                num_sc_reads += num_per_species[major_species]
        ## END of loop over barcodes

    summary_info = {}
    summary_info['num_sc_bcs'] = num_sc_bcs
    summary_info['num_sc_qual_reads'] = num_qual_reads
    summary_info['num_sc_reads'] = num_sc_reads
    summary_info['raw_bc_on_whitelist'] = raw_bc_on_whitelist

    barnyard_file.close()
    barnyard_hits_file.close()
    
    with open(outs.summary, 'w') as summary_file:
        summary_file.write(tenkit.safe_json.safe_jsonify(summary_info))

    with open(outs.barcode_histogram, 'w') as bc_hist_file:
        bc_hist_file.write(tenkit.safe_json.safe_jsonify(bc_hist))

############################### JOIN #######################################
def join(args, outs, chunk_defs, chunk_outs):
    num_sc_bcs = 0
    num_qual_reads = 0
    num_sc_reads = 0
    bc_counts = {}
    
    ## compute species_list
    ref = contig_manager.contig_manager(args.reference_path)
    species_list = ref.list_species()
    species_list.sort()

    ## doublet rate estimation
    total_unique_cell_barcodes = set()
    total_cell_barcodes = []
    species_counts = {}
    for (species, species_barcodes) in args.cell_barcodes.iteritems():
        species_counts[species] = 0
        for bc in species_barcodes.iterkeys():
            total_cell_barcodes.append(bc)
            total_unique_cell_barcodes.add(bc)
            species_counts[species] += 1
    counts = species_counts.values()

    observed_doublets = len(total_cell_barcodes) - len(total_unique_cell_barcodes)
    observed_doublet_rate = tk_stats.robust_divide(observed_doublets,
        float(len(total_cell_barcodes)))
    
    inferred_doublets = float('NaN')
    inferred_doublet_rate = float('NaN')
    if len(species_counts) > 1:
        inferred_doublets = _infer_multiplets_from_observed(observed_doublets,
            counts[0], counts[1])
        inferred_doublet_rate = tk_stats.robust_divide(float(inferred_doublets),
            float(len(total_cell_barcodes)))
    
    ## combine barnyard_hits chunks
    combine_csv([c.barnyard_hits for c in chunk_outs], outs.barnyard_hits,
                 header_lines=1)
    
    ## aggregate summary.json from chunks
    raw_bc_on_whitelist = 0
    for j,chunk_out in enumerate(chunk_outs):
        if chunk_out.summary is None: continue
        chunk_summary = json.loads(open(chunk_out.summary).read())
        num_sc_bcs += chunk_summary['num_sc_bcs']
        num_qual_reads += chunk_summary['num_sc_qual_reads']
        num_sc_reads += chunk_summary['num_sc_reads']
        raw_bc_on_whitelist += chunk_summary['raw_bc_on_whitelist']

        chunk_bc_counts_file = open(chunk_out.barcode_histogram)
        chunk_bc_counts = json.loads(chunk_bc_counts_file.read())
        bc_counts.update(chunk_bc_counts)

    ## combine barnyard chunks
    combine_csv([c.barnyard for c in chunk_outs], outs.barnyard,
                 header_lines=1)

    n_reads = np.array(bc_counts.values())
    max_val = np.percentile(n_reads, 99.99) * 1.3
    min_val = n_reads.min()
    num_bins = 400
    step = math.ceil((max_val - min_val)/num_bins)
    if max_val - min_val < 1e-6:
        bins = np.array([min_val, min_val+1])
    else:
        bins = np.arange(min_val, max_val, step)
    (hist, edges) = np.histogram(n_reads, bins=bins)
    bc_hist = {int(edges[i]):hist[i] for i in range(len(bins)-1)}

    cells = 0
    for (speci, cell_list) in args.cell_barcodes.iteritems():
        cells += len(cell_list)
    summary_info = {}
    summary_info['cells_detected'] = cells
    summary_info['num_sc_bcs'] = num_sc_bcs
    summary_info['num_sc_qual_reads'] = num_qual_reads
    summary_info['num_sc_reads'] = num_sc_reads
    summary_info['fract_sc_reads'] = tk_stats.robust_divide(num_sc_reads, num_qual_reads)
    summary_info['observed_doublets'] = observed_doublets
    summary_info['obserbed_doublet_rate'] = observed_doublet_rate
    summary_info['inferred_doublets'] = inferred_doublets
    summary_info['inferred_doublet_rate'] = inferred_doublet_rate
    
    ## compute stats from barnyard file
    barnyard_df = pd.read_csv( outs.barnyard )
    bkeys = ["amp_rate", "library_complexity", "dup_ratio", "mapped", "mapped_frac"]
    for species in species_list:
        if len(species_list) == 1:
            key_suffix = ""
        else:
            key_suffix = "_" + species

        is_cell_filter = barnyard_df["is_%s_cell_barcode"%species] == 1
        species_barcodes = args.cell_barcodes.get(species, {} )
        
        ## compute quartiles, min, mean, max and CV
        for bkey in bkeys:
            vals = barnyard_df[bkey][is_cell_filter]
            for pct in [25, 50, 75]:
                summary_key = bkey + key_suffix + ("_cells_p%d"%pct)
                summary_info[summary_key] = tk_stats.robust_percentile(vals, pct)
            summary_key = bkey + key_suffix + "_cells_cv"
            summary_info[summary_key] = tk_stats.robust_divide(vals.std(), vals.mean())
            summary_key = bkey + key_suffix + "_cells_mean"
            summary_info[summary_key] = vals.mean()
            summary_key = bkey + key_suffix + "_cells_min"
            summary_info[summary_key] = vals.min()
            summary_key = bkey + key_suffix + "_cells_max"
            summary_info[summary_key] = vals.max()

    ## tabulate waste metrics from barnyard_hits file
    waste_keys = ["no_barcode", "non_cell_barcode", "unmapped", 
                  "low_mapq_lt_%d"%PROFILE_MAPQ_THRESHOLD, 
                  "dups", "denominator", "unusable_read"]
    bh_df = pd.read_csv( outs.barnyard)

    # calculate median percent unmapped (defined as unmapped / (denominator - non_cell_barcode - no_barcode)
    banyard_cell_df = bh_df[~(bh_df.cell_id == 'None')]
    unmapped_frac = 1.0 * banyard_cell_df['unmapped'] / banyard_cell_df['denominator']
    unmapped_frac = unmapped_frac.fillna(0)
    median_unmapped_frac = unmapped_frac.median()

    waste_totals = {}
    sum_waste_keys = 0.0
    for key in waste_keys:
        waste_totals[key] = float(bh_df[key].sum( ))
        if key != "denominator":
            sum_waste_keys += waste_totals[key]
    for level, key in enumerate(waste_keys):
        if key == "denominator":
            continue
        summary_info["waste_%s_reads"%key] = waste_totals[key]
        summary_info["frac_waste_%s"%(key)] = tk_stats.robust_divide(
            waste_totals[key], waste_totals["denominator"] )
    summary_info["waste_total_reads"] = sum_waste_keys
    summary_info["frac_waste_total"] = tk_stats.robust_divide(
        sum_waste_keys, waste_totals["denominator"] )
    summary_info['frac_raw_bc_on_whitelist'] = float(raw_bc_on_whitelist)/waste_totals["denominator"]
    summary_info['median_unmapped_frac'] = median_unmapped_frac

    ## compute leakage metric and add to summary_info
    if len(species_list) == 2:
        compute_leakage( outs.barnyard_hits, ref, summary_info )
    
    with open(outs.summary, 'w') as summary_file:
        summary_file.write(tenkit.safe_json.safe_jsonify(summary_info,pretty=True))

    with open(outs.barcode_histogram, 'w') as bc_hist_file:
        bc_hist_file.write(tenkit.safe_json.safe_jsonify(bc_hist))

    # logging
    print tenkit.safe_json.safe_jsonify(summary_info, pretty=True)

########################### SUPPORT FUNCTIONS #################################

## TODO: this function belongs in tenkit, move it
def groupbybarcode( rec ):
    if rec.has_tag(PROCESSED_BARCODE_TAG):
        return rec.get_tag(PROCESSED_BARCODE_TAG)
    else:
        return "NO_BARCODE"
 
def _infer_multiplets_from_observed(n_obs_multiplets, n_cells0, n_cells1):
    if n_cells0 == 0 or n_cells1 == 0:
        return 0

    # Prior probability of a doublet given counts for each cell type (ignore N_cells > 2)
    p_obs_multiplet = ( 2*(float(n_cells0)/float(n_cells0+n_cells1))*
                        (float(n_cells1)/float(n_cells0+n_cells1)) )

    # Brute force MLE of binomial n
    n_mle = 0
    if n_obs_multiplets > 0:
        likelihood = scipy.stats.binom.pmf(n_obs_multiplets,
            xrange(0, n_cells0 + n_cells1), p_obs_multiplet)
        n_mle = np.argmax(likelihood)
    return n_mle


def read_has_barcode(r):
    bc = crdna_io.get_read_barcode(r)
    if bc is None:
        return False
    else:
        return True

## when to include a read in the barnyard hits file
def barnyard_hits_include( r ):
    if r.mapq < 60:
        return False
    if r.is_duplicate:
        return False
    if r.is_secondary:
        return False
    return True

## compute leakage from barnyard hits file only for barnyard samples
def compute_leakage( bh_file, ref, summary_dict ):
    ## for a given chromosome, fraction of that chromosome seen in
    ## non-species cell barcodes, compared to the mean number for a 
    ## species cell barcode
    def read_ratio2( df, species, other_species, chrom ):
        is_sp = "is_%s_cell_barcode" % species
        not_sp = "is_%s_cell_barcode" % other_species
        is_sp_filter = (df[is_sp] == 1)&(df[not_sp] == 0)
        is_not_filter = (df[is_sp] == 0)&(df[not_sp] == 1)
        if ( is_sp_filter.sum() == 0 or is_not_filter.sum() == 0 ):
            return []
        Rsp = float(df[is_sp_filter][chrom].mean( ))
        Rn  = df[is_not_filter][chrom].values
        ratio = Rn/max( Rsp, 1.0 )
        return ratio

    species_list = ref.list_species( )

    ## restrict to cell barcodes that are not doublets or empty.
    df = pd.read_csv( bh_file )
    df = df[(df.is_whitelisted == 1)]

    pctiles = [25, 50, 75]
    for spi, sp in enumerate(species_list):
        autosomes = ref.primary_contigs(species=sp, allow_sex_chromosomes=False)
        avals = np.array([read_ratio2( df, sp, species_list[1-spi], chrom )
                          for chrom in autosomes])
        autosome_mean = avals.mean(axis=0)
        
        mito_info = ref.non_nuclear_contigs() is not None
        if mito_info:
            mitosomes = [chrom for chrom in ref.non_nuclear_contigs()
                         if ref.species_from_contig(chrom) == sp]
            mito  = np.array([read_ratio2(df, sp, species_list[1-spi], mito_chrom )
                              for mito_chrom in mitosomes])
            mito_mean = mito.mean(axis=0)
        for i, pct in enumerate(pctiles):
            summary_dict["%s_autosome_leakage_frac_p%d"%(sp, pct)] = tk_stats.robust_percentile(autosome_mean, pct)
            summary_dict["%s_autosome_leakage_frac_mean"%(sp)] = np.mean(autosome_mean)           
            if mito_info:
                summary_dict["%s_mito_leakage_frac_p%d"%(sp, pct)] = tk_stats.robust_percentile(mito_mean, pct)
                summary_dict["%s_mito_leakage_frac_mean"%(sp)] = np.mean(mito_mean)  

