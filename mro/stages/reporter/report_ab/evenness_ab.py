# evenness_ab.py

import numpy
import random
import statsmodels.api as sm
import pandas

################################################################################
## Auxiliary function
which = lambda lst:list(numpy.where(lst)[0])

################################################################################
## Merge bins: from 20 kb to 100 kb resolution
def merge_bins(profiles, n_merge):
    # TODO: handle exceptions
    if n_merge < 2:
        # TODO: warning: illegal n_merge
        return(profiles)
    # if n_merge
    n_chrom = len(profiles)
    if n_chrom < 1:
        # TODO: warning: no chromosomes!
        return(profiles)
    # if n_chrom
    n_cells = len(profiles[0])
    if n_cells < 1:
        # TODO: warning: no cells!
        return(profiles)
    # if n_cells
    merged_profiles = [None] * n_chrom
    for chrom in range(n_chrom):
        n_bins = len(profiles[chrom][0])
        if n_bins < n_merge:
            # TODO: warning
            pass
        # if n_bins
        #n_merged_bins = int(math.ceil(n_bins / float(n_merge)))
        tmp = numpy.cumsum(profiles[chrom], axis=1, dtype=float)
        tmp[:, n_merge:] = tmp[:, n_merge:] - tmp[:, :-n_merge]
        merged_profiles[chrom] = tmp[:, n_merge - 1::n_merge] # This throws away any leftover bins!
    # for chrom
    return(merged_profiles)
# merge_bins

################################################################################
##
def merge_mask(mask, n_merge):
    if n_merge < 2:
        return(mask)
    # if n_merge
    n_chrom = len(mask)
    if n_chrom < 1:
        # TODO: warning: no chromosomes!
        return(mask)
    # if n_chrom
    merged_mask = [None] * n_chrom
    for chrom in range(n_chrom):
        n_bins = len(mask[chrom])
        if n_bins < n_merge:
            # TODO: warning
            pass
        # if n_bins
        merged_mask[chrom] = numpy.array([
            numpy.prod(mask[chrom][i:(i + n_merge)])
            for i in range(0, n_bins + 1, n_merge)
            ][:-1])
    # for chrom
    return(merged_mask)
# merge_mask

################################################################################
##
def get_total_read_counts(profiles, mask):
    total_read_counts = [None]
    n_chrom = len(profiles)  
    if n_chrom < 1:  
        # TODO: warning: no chromosomes!  
        return(total_read_counts)  
    # if n_chrom  
    n_cells = len(profiles[0])  
    total_read_counts = [0] * n_cells
    if n_cells < 1:  
        # TODO: warning: no cells!  
        return(total_read_counts)  
    # if n_cells  
    for cell in range(n_cells):
        total_read_counts[cell] = count_reads(profiles, cell, mask)
    # for cell
#    for chrom in range(n_chrom):
#        chr_read_counts = numpy.sum(numpy.transpose(profiles[chrom]) * mask[chrom], axis=0)
#        total_read_counts[i] = [total_read_counts[i] + chr_read_counts[i] 
#            for i in range(n_cells)]
#    # for chrom
    return(total_read_counts)  
# get_total_read_couns

################################################################################
## The same functionality is implemented in count_bins, probably faster
def get_total_n_bins(mask):
    total_n_bins = 0
    n_chrom = len(mask)  
    if n_chrom < 1:
        # TODO: warning: no chromosomes!  
        return(total_n_bins)
    # if n_chrom  
    for chrom in range(n_chrom):
        total_n_bins += sum(mask[chrom])
    # for chrom
    return(total_n_bins)
# get_total_n_bins

################################################################################
## Scale with respect to total reads
def scale(profiles, mask):
    n_chrom = len(profiles)  
    if n_chrom < 1:
        # TODO: warning: no chromosomes!  
        return(profiles)
    # if n_chrom  
    n_cells = len(profiles[0])
    if n_cells < 1:
        # TODO: warning: no cells!  
        return(profiles)
    # if n_cells  
    total_n_bins = get_total_n_bins(mask)
    total_read_counts = get_total_read_counts(profiles, mask)
    scaled_profiles = [None] * n_chrom
    mean_reads_per_bin = numpy.array([total_read_counts[i] / float(total_n_bins)
        for i in range(n_cells)])
    for chrom in range(n_chrom):
        scaled_profiles[chrom] = profiles[chrom] / mean_reads_per_bin[:, None]
    # for chr
    return(scaled_profiles)
# scale

################################################################################
## Preprocess profiles
def process_profiles(raw_profiles, mask, n_merge):
    # Merge bins: from 20 kb to 100 kb resolution
    profiles = merge_bins(raw_profiles, n_merge)
    merged_mask = merge_mask(mask, n_merge)
    # Don't scale with respect to total reads!
    # profiles = scale(profiles, merged_mask)
    return((profiles, merged_mask))
# process_profiles

################################################################################
## The same functionality implemented in g_total_n_bins
def count_bins(mask):
    n_bins = sum([sum(mask_chrom) for mask_chrom in mask])
    return(n_bins)
# count_bins

################################################################################
##
def count_reads(profiles, cell, mask):
    n_reads = sum([sum(profiles[chrom][cell] * mask[chrom])
         for chrom in range(len(mask))])
    return(n_reads)
# count_reads

################################################################################
##
def get_subcluster_read_count(profiles, subcluster, mask):
    total_read_count = 0
    bin_count = count_bins(mask)
    for cell in subcluster:
        read_count = count_reads(
            profiles, cell, mask)
        total_read_count += read_count
    # for profile_id
    # TODO: what if bin_count == 0
    return(total_read_count / float(bin_count))
# get_subcluster_read_count

################################################################################
#
def aggregate_cells(profiles, cluster):
    n_chrom = len(profiles)
    if n_chrom < 1:
        # TODO: Warning: no chromosomes
        return([[]])
    # if n_chrom
    chromosomes = range(n_chrom)
    combined_profile = [None] * n_chrom
    initialize = True
    count = 0
    for cell in cluster:
        if 0 <= cell < len(profiles[0]):
            count = count + 1
            if initialize:
                for chrom in chromosomes:
                    combined_profile[chrom] = numpy.array([profiles[chrom][cell]])
                # for chrom
                initialize = False
            else:
                for chrom in chromosomes:
                    # This needs to be cleaned up:
                    combined_profile[chrom] += profiles[chrom][cell]
                # for chromosome
            # if initialize else
        # if profiles
    # for cell
    # Don't even bother scaling these profiles!
    return(combined_profile)
# aggregate_cells

################################################################################
def ab_linear_model(n_reads_per_bin, sd):
    # TODO: check inputs!
    n_chrom = len(sd)
    if n_chrom < 1:
        # TODO: Warning: no chromosomes
        return((None, None))
    # if n_chrom
    y = [sd[i] * sd[i] * n_reads_per_bin[i] for i in range(len(n_reads_per_bin))]
    x = sm.add_constant(n_reads_per_bin)
    model = sm.OLS(y, x)
    results = model.fit()
    (b, a) = results.params
    #results.tvalues
    #results.t_test([1, 0]))
    #results.f_test(numpy.identity(2)))
    #
    #
    return((a, b))
# ab_linear_model

################################################################################
################################################################################
################################################################################

################################################################################
def fit_ab(merged_profiles, cluster, merged_mask, n_trials=200, seed=0):
    def generate_pair(n_cells):
        for i in range(n_cells - 1):
            for j in range(i + 1, n_cells):
                yield [i, j]
            # for j
        # for i
    # generate_pair
    n_chrom = len(merged_mask)
    results = [None] * n_chrom
    if n_chrom < 1:
        # TODO: Warning: no chromosomes
        return(results)
    # if n_chrom
    chromosomes = range(n_chrom)
    n_cells = len(merged_profiles[0])
    if n_cells < 1:
        # TODO: warning: no cells!
        return(results)
    # if n_cells
    pair = generate_pair(n_cells)
    n_singlets_and_pairs = 0.5 * n_cells * (n_cells + 1)
    n_trials += n_singlets_and_pairs
    n_trials = int(n_trials)
    random.seed(seed)
    cluster_size = len(cluster)
    sd = [None] * n_chrom
    for chrom in chromosomes:
        sd[chrom] = [None] * n_trials
    # for chrom
    n_reads_per_bin = [None] * n_trials
    for trial in range(n_trials):
        if trial < n_cells:
            subcluster = [trial]
        elif trial < n_singlets_and_pairs:
            subcluster = next(pair)
        else:
            subcluster_size = random.randint(1, cluster_size)
            subcluster = random.sample(cluster, subcluster_size)
        # if trial else
        n_reads_per_bin[trial] = get_subcluster_read_count(
            merged_profiles, subcluster, merged_mask)
        combined_profile = aggregate_cells(
            merged_profiles, cluster=subcluster)
        combined_profile = scale(combined_profile, merged_mask)
        for chrom in chromosomes:
            tmp = combined_profile[chrom][0][merged_mask[chrom] == 1]
            # Remove spikes:
            mean_value = 1;#numpy.mean(tmp)
            sd_value = numpy.nanstd(tmp) # Note: this is not robust when biological variability is present, needs to be accompanied by differential variance ratio
            z = (tmp - mean_value) / sd_value
            tmp = tmp[abs(z) < 6]
            tmp = tmp / numpy.mean(tmp)
            sd[chrom][trial] = numpy.nanstd(tmp)
        # for chrom
#        if trial % 10 == 0:
#            print("Completed %d/%d trials" % (trial, n_trials))
#        # if trial
    # for trial
    for chrom in chromosomes:
        results[chrom] = ab_linear_model(n_reads_per_bin, sd[chrom])
    # for chrom
    return(results)
# fit_ab

################################################################################
# Integrate all pieces
def evaluate_ab(raw_profiles, cluster, mask, n_trials, seed, n_merge):
    (merged_profiles, merged_mask) = process_profiles(raw_profiles, mask, n_merge)
    results = fit_ab(merged_profiles, cluster, merged_mask, n_trials, seed)
    return(results)
# evaluate_ab

