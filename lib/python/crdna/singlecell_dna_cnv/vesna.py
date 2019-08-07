import os
import math
import numpy as np
import pandas as pd
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from scipy.stats import poisson
import crdna.constants
import martian
#from sklearn import linear_model

TRAINED_PARAMETERS = "TRAINED_PARAMETERS"
REFERENCE_ANNOTATION = "REFERENCE_ANNOTATION"
NONE = "NONE"
GC_0 = 0.4083597

################################################################################
def get_default_bin_parameters(chroms, chrom_sizes, bin_size_normalized):
    bin_parameters = []
    for index, (chrom_size, chrom) in enumerate(zip(chrom_sizes, chroms)):
        # Get bin size, generate Start and End columns
        n_bins = chrom_size // bin_size_normalized
        start = np.arange(n_bins) * bin_size_normalized + 1
        end = start + bin_size_normalized - 1
        # Order columns as follows: ["Chr", "Start", "End", "GCFraction", 
        #     "NFraction", "InterceptQ", "LinearQ", "QuadraticQ", "Confidence"]
        bin_parameters.append(np.array([
            np.array([index + 1] * n_bins),
            start,
            end,
            np.zeros(n_bins), # GC content
            np.zeros(n_bins), # N-base content
            np.ones(n_bins), # Mappability
            np.zeros(n_bins), # Linear GC parameters
            np.zeros(n_bins), # Quadratic GC parameters
            np.ones(n_bins)])) # Confidence
    # for chrom
    return(bin_parameters)
# get_default_bin_parameters

################################################################################
def load_track_parameters(tracks, bin_parameters, ref, bin_size=20000):
    if tracks == None or not os.path.exists(tracks):
        martian.log_warn("Could not find %s" % tracks)
    # if tracks
    #
    # Sort chromosomes in ascending order, 
    # only use standard chromosomes, discard unplaced contigs:
    standard_chroms = ref.primary_contigs(allow_sex_chromosomes = True)
    n_standard_chroms = len(standard_chroms)
    #
    # Initialize mappability, N-base, and GC tracks:
    mappability = [None] * n_standard_chroms
    gc = [None] * n_standard_chroms
    n = [None] * n_standard_chroms
    confidence = [None] * n_standard_chroms
    #
    # Get mappability, N-base, and GC tracks:
    track_store = pd.HDFStore(tracks, 'r')
    for key in track_store.keys():
        tokens = key.split("/")
        if len(tokens) < 3:
            continue
        # if tokens
        (_, track, chrom) = tokens[:3]
        if chrom in standard_chroms:
            index = standard_chroms.index(chrom)
            if track == "map":
                mappability[index] = np.array(track_store[key])
            elif track == "GC":
                gc[index] = np.array(track_store[key])
            elif track == "N":
                n[index] = np.array(track_store[key])
            elif track == "CONF":
                confidence[index] = np.array(track_store[key])
            else:
                pass
            # if track else
        # if chrom
    # for key
    track_store.close()
    #
    # Evaluate average mappability (needed to scale the intercept):
    all_map = []
    for index in range(n_standard_chroms):
        all_map += mappability[index][which(n[index] <= crdna.constants.N_BASE_THRESHOLD)].tolist()
    # for index
    mean_map = np.nanmean(all_map)
    #
    for chrom in standard_chroms:
        index = standard_chroms.index(chrom)
        # Get bin size, generate Start and End columns
        n_bins = len(n[index])
        start = np.arange(n_bins) * bin_size + 1
        end = start + bin_size - 1
        # Generate QuadraticQ column
        g = gc[index] - GC_0
        g2 = g * g
        # Order columns as follows: ["Chr", "Start", "End", "GCFraction", 
        #     "NFraction", "InterceptQ", "LinearQ", "QuadraticQ", "Confidence"]
        bin_parameters.append(np.array([
            np.array([standard_chroms.index(chrom) + 1] * n_bins),
            start,
            end,
            gc[index],
            n[index],
            mappability[index] / mean_map, 
            g,
            g2,
            confidence[index]]))
    # for chrom
    return(bin_parameters)
# load_track_parameters

################################################################################
def running_average(array, n_merge):
    tmp = np.cumsum(array, dtype=float)
    tmp = (tmp[n_merge:] - tmp[:-n_merge]) / n_merge
    merged_array = np.array([np.nanmean(array[:n_merge])] + 
        tmp[n_merge - 1::n_merge].tolist()) 
    return(merged_array)
# running_average

################################################################################
def merge_track_bins(bin_parameters, n_merge):
    n_chrom = len(bin_parameters)
    merged_bin_parameters = [None] * n_chrom
    old_start = bin_parameters[0][1]
    old_bin_size = old_start[1] - old_start[0] + 1
    merged_bin_size = old_bin_size * n_merge
    gc_column = 3
    n_column = 4
    map_column = 5
    conf_column = 8
    for chrom in range(n_chrom):
        n_bins = bin_parameters[chrom][0].shape[0]
        n_merged_bins = int(n_bins / n_merge)
        start = np.arange(n_merged_bins) * merged_bin_size + 1
        end = start + merged_bin_size - 1
        merged_g = running_average(bin_parameters[chrom][gc_column], n_merge)
        merged_g_g0 = merged_g - GC_0
        merged_n = running_average(bin_parameters[chrom][n_column], n_merge)
        merged_mappability = running_average(bin_parameters[chrom][map_column], n_merge)
        merged_confidence = running_average(bin_parameters[chrom][conf_column], n_merge)
        g2 = merged_g_g0 * merged_g_g0
        merged_bin_parameters[chrom] = [
            np.array([chrom] * n_merged_bins),
            start,
            end,
            merged_g,
            merged_n,
            merged_mappability,
            merged_g_g0,
            g2,
            merged_confidence]
    # for chrom
    return(merged_bin_parameters)
# merge_track_bins

################################################################################
def get_profile(profiles, cell):
    n_chrom = get_n_chrom(profiles)
    return([profiles[chrom][cell, :] for chrom in range(n_chrom)]);
# get_profile

################################################################################
## Auxiliary function
which = lambda lst:list(np.where(lst)[0])

################################################################################
## Load bin parameters
def load_bin_parameters(file_name, bin_parameters, chroms):
    store = pd.HDFStore(file_name, "r")
    for chrom in chroms:
        bin_parameters.append(np.array(store.get_node(chrom)))
    # for chrom
    store.close()
# load_bin_parameters

################################################################################
## Load data and masks
def load_data(file_name, raw_profiles, n_mask, chroms):
    store = pd.HDFStore(file_name, "r")
    for chrom in chroms:
        raw_profiles.append(np.array(store[chrom]))
        tmp = store["mask_%s"%chrom]
        n_bins = len(tmp)
        tmp = np.array([True] * n_bins)
        n_mask.append(tmp)
    # for chrom
    store.close()
# load_data

################################################################################
def get_n_chrom(profiles):
    return(len(profiles))
# get_n_chrom

################################################################################
def merge_bins(profiles, n_merge, average=False):
    # TODO: handle exceptions
    if n_merge < 2:
        # TODO: warning: illegal n_merge
        return(profiles)
    # if n_merge
    n_chrom = get_n_chrom(profiles)
    if n_chrom < 1:
        # TODO: warning: no chromosomes!
        return(profiles)
    # if n_chrom
    n_cells = profiles[0].shape[0]
    if n_cells < 1:
        # TODO: warning: no cells!
        return(profiles)
    # if n_cells
    merged_profiles = [None] * n_chrom
    for chrom in range(n_chrom):
        n_bins = profiles[chrom].shape[1]
        if n_bins < n_merge:
            # TODO: warning
            pass
        # if n_bins
        #n_merged_bins = int(math.ceil(n_bins / float(n_merge)))

        for cell in range(n_cells):
            profile_selector = np.where(np.isnan(profiles[chrom][cell, :]))[0]
            profile_selector = profile_selector.tolist()
            profiles[chrom][cell, profile_selector] = 0
        # for cell

        tmp = np.cumsum(profiles[chrom], axis=1, dtype=float)
        tmp[:, n_merge:] = tmp[:, n_merge:] - tmp[:, :-n_merge]
        merged_profiles[chrom] = tmp[:, n_merge - 1::n_merge] # This throws away any leftover bins!
        if average:
            merged_profiles[chrom] /= n_merge
        # if average
    # for chrom
    return(merged_profiles)
# merge_bins

################################################################################
##
def merge_mask(mask, n_merge):
    if n_merge < 2:
        return(mask)
    # if n_merge
    n_chrom = get_n_chrom(mask)
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
        merged_mask[chrom] = np.array([
            np.prod(mask[chrom][i:(i + n_merge)])
            for i in range(0, n_bins + 1, n_merge)
            ][:-1])
    # for chrom
    return(merged_mask)
# merge_mask

################################################################################
##
def get_total_read_counts(profiles, mask):
    total_read_counts = [None]
    n_chrom = get_n_chrom(profiles)
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
    for chrom in range(n_chrom):
        chr_read_counts = np.sum(np.transpose(profiles[chrom])[which(mask[chrom])], axis=0)
        total_read_counts = [total_read_counts[cell] + chr_read_counts[cell]
            for cell in range(n_cells)]
    # for chrom
    return(total_read_counts)
# get_total_read_couns
        
################################################################################
##          
def get_total_n_bins(mask):
    total_n_bins = 0
    n_chrom = get_n_chrom(mask)
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
    n_chrom = get_n_chrom(profiles)
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
    mean_reads_per_bin = np.array([total_read_counts[cell] / total_n_bins
        for cell in range(n_cells)])
    for chrom in range(n_chrom):
        scaled_profiles[chrom] = profiles[chrom] / mean_reads_per_bin[:, None]
    # for chr
    return((total_read_counts, mean_reads_per_bin, scaled_profiles))
# scale

################################################################################
def get_common_n_bins_per_chrom(a, b): 
    n_chrom = get_n_chrom(a)
    assert(n_chrom == get_n_chrom(b))
    return([np.min([len(a[chrom]), len(b[chrom])])
                     for chrom in range(n_chrom)])
# get_common_n_bins_per_chrom

################################################################################
## 
def create_confident_mask(conf_df, n_chrom):
    confident_mask = []
    for chrom in range(n_chrom):
        selector = which(conf_df["chrom"] == "hg19_chr" + str(chrom + 1))
        tmp = (conf_df["perc"].values[selector] > 0.70)
        confident_mask.append(tmp)
    # for chrom
    return(confident_mask)
# create_confident_mask

################################################################################
def combine_masks(mask_A, mask_B):
    n_chrom = get_n_chrom(mask_A)
    assert(n_chrom > 0)
    assert(n_chrom == get_n_chrom(mask_B))
    n_bins = get_common_n_bins_per_chrom(mask_A, mask_B)
    combined_mask = [[mask_A[chrom][bin_index] and 
                      mask_B[chrom][bin_index] 
                     for bin_index in range(n_bins[chrom])]
                     for chrom in range(n_chrom)]
    return(combined_mask)
# combine_masks

################################################################################
## Preprocess profiles
def process_profiles(raw_profiles, n_mask, confident_mask, n_merge):
    # Merge bins: from 20 kb to 100 kb resolution
    merged_profiles = merge_bins(raw_profiles, n_merge)
    merged_n_mask = merge_mask(n_mask, n_merge)
    combined_mask = combine_masks(merged_n_mask, confident_mask)
    # Scale with respect to total reads
    (total_read_counts, mean_reads_per_bin, scaled_profiles) = scale(
        merged_profiles, combined_mask)
    return(total_read_counts, mean_reads_per_bin, 
           scaled_profiles, merged_profiles, combined_mask)
# process_profiles

################################################################################
def get_gc_bias_coefficients(
        profiles, combined_mask, bin_parameters, 
        low_gc=0.3, high_gc=0.7):
    # Get GC content per bin
    gc_column = 3
    n_chrom = get_n_chrom(profiles)
    n_bins = [profiles[chrom].shape[1] for chrom in range(n_chrom)]
    n_cells = profiles[0].shape[0]
    gc = np.concatenate([bin_parameters[chrom][gc_column,:n_bins[chrom]] 
                         for chrom in range(n_chrom)], axis=0)
    m = np.concatenate([combined_mask[chrom][:n_bins[chrom]]
                        for chrom in range(n_chrom)], axis=0)
    flag = [((low_gc <= g) and (g <= high_gc)) for g in gc]
    flag = flag and m
    gc = gc[which(flag)]
    gc0 = 0.4083597 # np.mean(gc)
    gc = gc - gc0
    poly = PolynomialFeatures(degree=2, include_bias=False)
    X = gc[:,np.newaxis]
    poly.fit(X)
    X_poly = poly.transform(X)
    gc_bias_coefficients = []
    for cell in range(n_cells):
        counts = np.concatenate([profiles[chrom][cell,:] 
                                 for chrom in range(n_chrom)], axis=0)
        y = counts[which(flag)]
        reg = LinearRegression().fit(X_poly, y)
        gc_bias_coefficients.append(
            [reg.intercept_, reg.coef_[0], reg.coef_[1]])
    # for cell
    return(gc_bias_coefficients)
# get_gc_bias_coefficients

################################################################################
def get_column(column_name):
    column_names = ["Chr", "Start", "End", "GCFraction", "NFraction", 
                    "InterceptQ", "LinearQ", "QuadraticQ"]
    assert(column_name in column_names)
    return(column_names.index(column_name))
# get_column

################################################################################
#def get_gc_bias_coefficients(
#        profiles, combined_mask, bin_parameters, 
#        low_gc=0.3, high_gc=0.7):
#    # Get GC content per bin
#    gc_column = get_column("GCFraction")
#    intercept_column = get_column("InterceptQ")
#    n_chrom = get_n_chrom(profiles)
#    n_bins = [profiles[chrom].shape[1] for chrom in range(n_chrom)]
#    n_cells = profiles[0].shape[0]
#    gc = np.concatenate([bin_parameters[chrom][gc_column, :n_bins[chrom]] 
#                         for chrom in range(n_chrom)], axis=0)
#    intercept = np.concatenate([bin_parameters[chrom][intercept_column, :n_bins[chrom]]
#                        for chrom in range(n_chrom)], axis=0)
#    m = np.concatenate([combined_mask[chrom][:n_bins[chrom]]
#                        for chrom in range(n_chrom)], axis=0)
#    flag = [((low_gc <= g) and (g <= high_gc)) for g in gc]
#    flag = flag and m
#    gc = gc[which(flag)]
#    intercept = intercept[which(flag)]
#    gc0 = 0.409519 # np.mean(gc) # Using 0.409519 instead of mean GC because I forgot to mask GC when training bin parameters
#    gc = gc - gc0
#    X1 = (gc / intercept)[:,np.newaxis]
#    X2 = (gc * gc / intercept)[:,np.newaxis]
#    X = np.hstack((X1, X2))
#    gc_bias_coefficients = []
#    for cell in range(n_cells):
#        counts = np.concatenate([profiles[chrom][cell,:] 
#                                 for chrom in range(n_chrom)], axis=0)
#        y = counts[which(flag)] / intercept
#        regression = LinearRegression().fit(X, y)
#        print(regression.intercept_)
#        print(regression.coef_)
#        gc_bias_coefficients.append(
#            [regression.intercept_, 
#             regression.coef_[0], 
#             regression.coef_[1]])
#    # for cell
#    return(gc_bias_coefficients)
## get_gc_bias_coefficients

################################################################################
## Integrate all pieces
def get_posterior_ploidy_probability(
        raw_profiles, n_mask, n_merge, conf_df, bin_parameters, max_ploidy=6):
    n_chrom = get_n_chrom(raw_profiles)
    confident_mask = create_confident_mask(conf_df, n_chrom)
    (total_read_counts, mean_reads_per_bin, 
     scaled_profiles, merged_profiles, combined_mask) = process_profiles(
        raw_profiles, n_mask, confident_mask, n_merge)
    gc_bias_coefficients = get_gc_bias_coefficients(
        scaled_profiles, combined_mask, bin_parameters)
    ploidy_probability = []
    n_cells = scaled_profiles[0].shape[0]
    dummy_cell = None
    for cell in range(n_cells):
        try:
            scaled_profile = get_profile(scaled_profiles, cell)
            merged_profile = get_profile(merged_profiles, cell)
            p_matrix = create_p_matrix(scaled_profile, max_ploidy)    
            prior = initialize_prior(scaled_profile, combined_mask, p_matrix)
            p_matrix = initialize_posterior(
                merged_profile, 
                combined_mask, 
                mean_reads_per_bin[cell], 
                gc_bias_coefficients[cell],
                bin_parameters, p_matrix, prior)
            ploidy_probability.append(p_matrix)
        except Exception as error:
            print('[get_posterior_ploidy_probability] caught this error: ' + repr(error))
            ploidy_probability.append(dummy_cell)
        # try except
        #print("Completed cell %d/%d" % (cell, n_cells))
    # for cell
    return(ploidy_probability)
# get_posterior_ploidy_probability

################################################################################
def get_posterior_ploidy_probability_from_intermediate_data(max_ploidy,
            merged_profiles, scaled_profiles, poisson_expectations, combined_mask):
    ploidy_probability = []
    n_cells = scaled_profiles[0].shape[0]
    dummy_cell = None
    for cell in range(n_cells):
        try:
            if poisson_expectations[cell] == None:
                pass
            # if None
            scaled_profile = get_profile(scaled_profiles, cell)
            merged_profile = get_profile(merged_profiles, cell)
            poisson_expectations_single_cell = poisson_expectations[cell]
            p_matrix = create_p_matrix(scaled_profile, max_ploidy)
            prior = initialize_prior(scaled_profile, combined_mask, p_matrix)
            p_matrix = initialize_posterior_from_poisson_expectations_single_cell(
                p_matrix, merged_profile, poisson_expectations_single_cell, combined_mask, prior);
            ploidy_probability.append(p_matrix)
        except Exception as error:
            print('[get_posterior_ploidy_probability_from_intermediate_data] caught this error: ' + repr(error))
            ploidy_probability.append(dummy_cell)
        # try except
        #print("Completed cell %d/%d" % (cell, n_cells))
    # for cell
    return(ploidy_probability)
# get_posterior_ploidy_probability_from_intermediate_data

################################################################################
def create_p_matrix(profile, max_ploidy=6):
    # Set up probability matrix: for each bin, ploidies 0 to maxPloidy
    n_chrom = get_n_chrom(profile)
    p_matrix = []
    for chrom in range(n_chrom):
        n_bins = profile[chrom].shape[0]
        p_matrix.append(np.zeros((n_bins, max_ploidy + 1)))
    # for chr
    return(p_matrix)
# create_p_matrix

################################################################################
def get_max_ploidy(p_matrix):
    max_ploidy = p_matrix[0].shape[1] - 1
    return(max_ploidy)
# get_max_ploidy

################################################################################
def initialize_prior(scaled_profile, combined_mask, p_matrix):
    # Set up prior distribution based on histogram of scaled counts
    eps = 1e-6
    n_chrom = get_n_chrom(scaled_profile)
    max_ploidy = get_max_ploidy(p_matrix)
    prior = []
    for chrom in range(n_chrom):
        chr_prior = np.zeros(max_ploidy + 1)
        lower_bound = 0
        upper_bound = lower_bound + 0.25
        n_bins = scaled_profile[chrom].shape[0]
        n_unmasked_bins = float(len(which(combined_mask[chrom]))) + eps
        for ploidy in range(max_ploidy + 1):
            n_inside = len(which(
                [(lower_bound <= scaled_profile[chrom][bin_index]) and
                 (scaled_profile[chrom][bin_index] <= upper_bound) and
                 combined_mask[chrom][bin_index]
                 for bin_index in range(n_bins)]))
            lower_bound = upper_bound
            upper_bound += 0.5
            chr_prior[ploidy] = n_inside / n_unmasked_bins 
        # for ploidy
        prior.append(chr_prior)
    # for chrom
    return(prior)
# initialize_prior
# 
################################################################################
def initialize_posterior(
            merged_profile, combined_mask, 
            average_reads_per_bin, gc_bias_coefficients,
            bin_parameters, p_matrix, prior):
    # Initialize posterior probabilities (each bin, ploidies 0 to maxPloidy):
    n_chrom = get_n_chrom(p_matrix)
    n_bins = [merged_profile[chrom].shape[0] for chrom in range(n_chrom)]
    max_ploidy = get_max_ploidy(p_matrix)
    dummy = [None] * (max_ploidy + 1)
    intercept_column = get_column("InterceptQ")
    linear_column = get_column("LinearQ")
    quadratic_column = get_column("QuadraticQ")
    for chrom in range(n_chrom):
        n_chr_bins = n_bins[chrom]
        q0 = bin_parameters[chrom][intercept_column, :n_chr_bins]
        q1 = bin_parameters[chrom][linear_column, :n_chr_bins]
        q2 = bin_parameters[chrom][quadratic_column, :n_chr_bins]
        for bin_index in range(n_chr_bins):
            if combined_mask[chrom][bin_index]:
                mu_diploid = average_reads_per_bin * (
                    q0[bin_index] + 
                    q1[bin_index] * gc_bias_coefficients[1] + 
                    q2[bin_index] * gc_bias_coefficients[2])
                mu = mu_diploid * 0.5 * np.arange(max_ploidy + 1)
                p_counts = [prior[chrom][ploidy] *
                            poisson.pmf(merged_profile[chrom][bin_index], mu[ploidy])
                            for ploidy in range(max_ploidy + 1)]
                for ploidy in range(max_ploidy + 1):
                    if math.isnan(p_counts[ploidy]):
                        p_counts[ploidy] = 0.0
                    # if p_counts
                # for ploidy
                z = sum(p_counts)
                p_matrix[chrom][bin_index, :] = [
                    p_counts[ploidy] / z for ploidy in range(max_ploidy + 1)]
            else:
                p_matrix[chrom][bin_index, :] = dummy
            # if combined_mask else
        # for bin
    # for chrom
    return(p_matrix)
# initialize_posterior

################################################################################
def get_poisson_expectations_single_cell(
        p_matrix_single_cell, average_reads_per_bin_single_cell, 
        gc_bias_coefficients_single_cell, combined_mask, bin_parameters):
    #print('Entering get_poisson_expectations_single_cell()')
    #print('p_matrix_single_cell.shape=%d,%d,%d'%(
    #      len(p_matrix_single_cell),
    #      len(p_matrix_single_cell[0]),
    #      len(p_matrix_single_cell[0][0])))
    n_chrom = get_n_chrom(p_matrix_single_cell)
    n_bins = [p_matrix_single_cell[chrom].shape[0] for chrom in range(n_chrom)]
    max_ploidy = get_max_ploidy(p_matrix_single_cell)
    intercept_column = get_column("InterceptQ")
    linear_column = get_column("LinearQ")
    quadratic_column = get_column("QuadraticQ")
    poisson_expectations = np.zeros_like(p_matrix_single_cell)
    for chrom in range(n_chrom):
        n_chr_bins = n_bins[chrom]
        q0 = bin_parameters[chrom][intercept_column, :n_chr_bins]
        q1 = bin_parameters[chrom][linear_column, :n_chr_bins]
        q2 = bin_parameters[chrom][quadratic_column, :n_chr_bins]
        poisson_expectations[chrom] = [np.array([None] * (max_ploidy + 1))] * n_chr_bins
        for bin_index in range(n_chr_bins):
            if combined_mask[chrom][bin_index]:
                mu_diploid = average_reads_per_bin_single_cell * (
                    q0[bin_index] +
                    q1[bin_index] * gc_bias_coefficients_single_cell[1] +
                    q2[bin_index] * gc_bias_coefficients_single_cell[2])
                mu = mu_diploid * 0.5 * np.arange(max_ploidy + 1)
                poisson_expectations[chrom][bin_index] = mu
            # if combined_mask
        # for bin_index
    # for chrom
    #print('Leaving get_poisson_expectations_single_cell()')
    return(poisson_expectations)
# get_poisson_expectations_single_cell

################################################################################
def get_poisson_expectations(average_reads_per_bin,
        merged_profiles, gc_bias_coefficients, combined_mask, bin_parameters, max_ploidy):
    #print('Entering get_poisson_expectations()')
    n_cells = len(average_reads_per_bin)
    poisson_expectations = [None] * n_cells
    dummy_cell = None
    for cell in range(n_cells):
        try:
            merged_profile = get_profile(merged_profiles, cell)
            p_matrix_single_cell = create_p_matrix(merged_profile, max_ploidy)
            poisson_expectations_single_cell = get_poisson_expectations_single_cell(
                p_matrix_single_cell, average_reads_per_bin[cell],
                gc_bias_coefficients[cell], combined_mask, bin_parameters)
            poisson_expectations[cell] = poisson_expectations_single_cell
        except Exception as error:
            print('[get_poisson_expectations] caught this error: ' + repr(error))
            poisson_expectations[cell] = dummy_cell
        # try except
        print('Completed cell %d/%d' % (cell, n_cells))
    # for cell
    #print('Leaving get_poisson_expectations()')
    return(poisson_expectations)
# get_poisson_expectations

################################################################################
def initialize_posterior_from_poisson_expectations_single_cell(
        p_matrix, merged_profile, poisson_expectations_single_cell, combined_mask, prior):
    eps = 1e-6 # Serves to prevent division by zero
    # Initialize posterior probabilities (each bin, ploidies 0 to maxPloidy):
    n_chrom = get_n_chrom(p_matrix)
    n_bins = [merged_profile[chrom].shape[0] for chrom in range(n_chrom)]
    max_ploidy = get_max_ploidy(p_matrix)
    dummy = [None] * (max_ploidy + 1)
    for chrom in range(n_chrom):
        n_chr_bins = n_bins[chrom]
        for bin_index in range(n_chr_bins):
            if combined_mask[chrom][bin_index]:
                mu = poisson_expectations_single_cell[chrom][bin_index]
                p_counts = [prior[chrom][ploidy] *
                            poisson.pmf(merged_profile[chrom][bin_index], mu[ploidy])
                            for ploidy in range(max_ploidy + 1)]
                for ploidy in range(max_ploidy + 1):
                    if math.isnan(p_counts[ploidy]):
                        p_counts[ploidy] = 0.0
                    # if p_counts
                # for ploidy
                z = sum(p_counts) + eps 
                p_matrix[chrom][bin_index, :] = [
                    p_counts[ploidy] / z for ploidy in range(max_ploidy + 1)]
            else:
                p_matrix[chrom][bin_index, :] = dummy
            # if combined_mask else
        # for bin
    # for chrom
    return(p_matrix)
# initialize_posterior_from_poisson_expectations_single_cell

################################################################################
def get_max_likelihood_ploidy_per_cell(p_matrix):
    n_chrom = get_n_chrom(p_matrix)
    n_bins = [p_matrix[chrom].shape[0] for chrom in range(n_chrom)]
    max_likelihood_ploidy_per_cell = [None] * n_chrom
    for chrom in range(n_chrom):
        tmp = [float("nan")] * n_bins[chrom]
        for bin_index in range(n_bins[chrom]):
            if math.isnan(p_matrix[chrom][bin_index][0]):
                pass
            else:
                tmp[bin_index] = np.argmax(p_matrix[chrom][bin_index, :])
            # if NaN else
        # for bin_index
        max_likelihood_ploidy_per_cell[chrom] = tmp
    # for chrom
    return(max_likelihood_ploidy_per_cell)
# get_max_likelihood_ploidy_per_cell

################################################################################
def get_max_likelihood_ploidies(ploidy_probability):
    n_cells = len(ploidy_probability)
    max_likelihood_ploidies = [None] * n_cells
    for cell in range(n_cells):
        try:
            if ploidy_probability[cell] == None:
                print("No data for cell %d" % cell)
                continue
            # if None
            p_matrix = ploidy_probability[cell]
            max_likelihood_ploidies[cell] = get_max_likelihood_ploidy_per_cell(p_matrix)
        except Exception as error:
            print('get_max_likelihood_ploidies] caught this error: ' + repr(error))
        # try except
    # for cell
    return(max_likelihood_ploidies)
#  get_max_likelihood_ploidies

# ################################################################################
# getMeanVarPloidy <- function( pMatrix ) {
#     meanVarPloidy <- list();
#     meanVarPloidy$Mean <- list();
#     meanVarPloidy$Var <- list();
#     nChr <- length( pMatrix );
#     maxPloidy <- getMaxPloidy( pMatrix );
#     ploidy <- 0:maxPloidy;
#     for ( chr in 1:nChr ) {
#         chrMean <- apply( pMatrix[[ chr ]], 1, function( x ) sum( x * ploidy ) );
#         chrSecondMoment <- apply( pMatrix[[ chr ]], 1, function( x ) sum( x * ploidy * ploidy ) );
#         chrVar <- chrSecondMoment - chrMean * chrMean;
#         meanVarPloidy$Mean[[ chr ]] <- chrMean;
#         meanVarPloidy$Var[[ chr ]] <- chrVar;
#     } # for chr
#     return( meanVarPloidy );
# } # getMeanVarPloidy
def get_mean_var_ploidy_single_cell(p_matrix):
    mean_ploidy = []
    var_ploidy = []
    max_ploidy = get_max_ploidy(p_matrix)
    n_chrom = get_n_chrom(p_matrix)
    n_bins = [p_matrix[chrom].shape[0] for chrom in range(n_chrom)]
    for chrom in range(n_chrom):
        chr_mean = np.array([sum([ploidy * p_matrix[chrom][bin_index][ploidy] 
                                  for ploidy in range(max_ploidy)]) 
                             for bin_index in range(n_bins[chrom])])
        chr_second_moment = np.array([sum([ploidy * ploidy * p_matrix[chrom][bin_index][ploidy] 
                                  for ploidy in range(max_ploidy)]) 
                             for bin_index in range(n_bins[chrom])])
        chr_var = chr_second_moment - chr_mean * chr_mean
        mean_ploidy.append(chr_mean)
        var_ploidy.append(chr_var)
    # for chrom
    return((mean_ploidy, var_ploidy))
# get_mean_var_ploidy_single_cell

################################################################################
def get_mean_var_ploidy(ploidy_probability):
    n_cells = len(ploidy_probability)
    mean_var_ploidy = [(None, None)] * n_cells
    for cell in range(n_cells):
        if ploidy_probability[cell] == None:
            print("No data for cell %d" % cell)
            continue
        # if None
        mean_var_ploidy[cell] = get_mean_var_ploidy_single_cell(
            ploidy_probability[cell])
    # for cell
    return(mean_var_ploidy)
# get_mean_var_ploidy

################################################################################
def get_ratio(merged_profiles, poisson_expectations, combined_mask, mean_reads_per_bin):
    ratio = merged_profiles
    n_chrom = merged_profiles.shape[0]
    n_cells = merged_profiles[0].shape[0]
    for chrom in range(n_chrom):
        for cell in range(n_cells):
            n_bins = len(merged_profiles[chrom][cell])
            tmp = [poisson_expectations[cell][chrom][bin_index][2] 
                for bin_index in range(n_bins)]
            tmp = [float('NaN') if tmp[bin_index] == None else tmp[bin_index] 
                for bin_index in range(n_bins)]
            diploid_expectation = np.array(tmp)
            ratio[chrom][cell] = merged_profiles[chrom][cell] / diploid_expectation
        # for cell
    # for chrom
    return(np.array(ratio))
# get_ratio

# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# # Iteratively refine posterior matrix.
# refinePMatrix <- function( pMatrix, maskedProfile, averageReadsPerBin, gcCoefficients, 
#                            binParameters, characteristicLength=5, extent=10 ) {
#     tmp <- dexp( 0:extent, 1 / characteristicLength );
#     to <- tmp[ 1 ];
#     tmp <- tail( tmp, -1 );
#     p <- c( rev( tmp ), to, tmp );
#     p <- p / sum( p );
#     newPMatrix <- pMatrix;
#     maxPloidy <- getMaxPloidy( pMatrix );
#     nChr <- length( pMatrix );
#     for ( chr in 1:nChr ) {
#         q0 <- binParameters[[ chr ]]$InterceptQ;
#         q1 <- binParameters[[ chr ]]$LinearQ;
#         q2 <- binParameters[[ chr ]]$QuadraticQ;
#         nBins <- nrow( pMatrix[[ chr ]] );
#         for ( bin in 1:nBins ) {
#             # First reevaluate prior:
#             start <- max( bin - extent, 1 );
#             if ( start == 1 ) {
#                 trimLeft <- extent - bin + 1;
#             } else {
#                 trimLeft <- 0;
#             } # if start else
#             end <- min( bin + extent, nBins );
#             if ( end == nBins ) {
#                 trimRight <- extent - ( nBins - bin );
#             } else {
#                 trimRight <- 0;
#             } # if start else
#             weights <- p[ ( 1 + trimLeft ):( 2 * extent + 1 - trimRight ) ];
#             prior <- rep( NA, maxPloidy + 1 );
#             for ( ploidy in 0:maxPloidy ) {
#                 prior[[ ploidy + 1 ]] <- sum( pMatrix[[ chr ]][ start:end, ploidy + 1 ] * weights, na.rm=TRUE );
#             } # for ploidy
#             prior <- prior / sum( prior, na.rm=TRUE );
#             #
#             muDiploid <- averageReadsPerBin * 
#                 ( q0[[ bin ]] + q1[[ bin ]] * gcCoefficients$IQ1 + q2[[ bin ]] * gcCoefficients$IQ2 );
#             muPloidy <- muDiploid * 0.5 * ( 0:maxPloidy );
#             pCountsConditionalOnPloidy <- sapply( muPloidy, function( x ) {
#                 dpois( maskedProfile[[ chr ]][[ bin ]], x )
#             } );
#             tmp <- pCountsConditionalOnPloidy * prior;
#             newPMatrix[[ chr ]][ bin, ] <- tmp / sum( tmp );
#         } # for bin
#     } # for chr
#     return( newPMatrix );
# } # refinePMatrix
# 
# ################################################################################
# getScorePerBin <- function( 
#     maskedProfile, averageReadsPerBin, pMatrix,
#     binParameters, gcCoefficients, chr, bin, 
#     testLevels, nullLevels ) {
#     maxPloidy <- length( pMatrix );
#     q0 <- binParameters[[ chr ]]$InterceptQ[[ bin ]];
#     q1 <- binParameters[[ chr ]]$LinearQ[[ bin ]];
#     q2 <- binParameters[[ chr ]]$QuadraticQ[[ bin ]];
#     muDiploid <- averageReadsPerBin * 
#         ( q0 + q1 * gcCoefficients$IQ1 + q2 * gcCoefficients$IQ2 );
#     muPloidy <- muDiploid * 0.5 * ( 0:maxPloidy );
#     pCountsConditionalOnPloidy <- sapply( muPloidy, function( x ) {
#         dpois( maskedProfile[[ chr ]][[ bin ]], x )
#     } );
#     pTest <- sum( pCountsConditionalOnPloidy[ testLevels ] );
#     pNull <- sum( pCountsConditionalOnPloidy[ nullLevels ] );
#     score <- log10( pTest / pNull );
#     return( score );
# } # getScorePerBin
# 
# ################################################################################
# getDeletionScorePerBin <- function( 
#     maskedProfile, averageReadsPerBin, pMatrix,
#     binParameters, gcCoefficients, chr, bin ) {
#     maxPloidy <- length( pMatrix );
#     homozygousDeletion <- 1;
#     heterozygousDeletion <- 2;
#     diploid <- 3;
#     return( getScorePerBin( 
#         maskedProfile, averageReadsPerBin, pMatrix,
#         binParameters, gcCoefficients, chr, bin, 
#         testLevels=homozygousDeletion:heterozygousDeletion, 
#         nullLevels=diploid:maxPloidy ) );
# } # getDeletionScorePerBin
# 
# ################################################################################
# getDuplicationScorePerBin <- function( 
#     maskedProfile, averageReadsPerBin, pMatrix, 
#     binParameters, gcCoefficients, chr, bin ) {
#     maxPloidy <- length( pMatrix );
#     homozygousDeletion <- 1;
#     heterozygousDuplication <- 4;
#     diploid <- 3;
#     return( getScorePerBin( 
#         maskedProfile, averageReadsPerBin, pMatrix,
#         binParameters, gcCoefficients, chr, bin, 
#         testLevels=heterozygousDuplication:maxPloidy, 
#         nullLevels=homozygousDeletion:diploid ) );
# } # getDuplicationScorePerBin
# 
# ################################################################################
# visualizePosterior <- function( results, chr, what="MaxLikelihoodPloidy", how="Line" ) {
#     jet.colors <- colorRampPalette( 
#             c( "#00007F", "blue", "#007FFF", "cyan",
#                "#7FFF7F", "yellow", "#FF7F00", "red", 
#                "#7F0000" ) );
#     titleText <- paste0( "Chr", chr );
#     maxPloidy <- getMaxPloidy( results$PMatrix );
#     if ( as.character( what ) == "MaxLikelihoodPloidy" ) {
#         if ( as.character( how ) == "Line" ) {
#             plot( results$MaxLikelihoodPloidy[[ chr ]], type="l", ylim=c( 0, maxPloidy ),
#                   xlab="Bins", ylab="Maximum Likelihood Ploidy", 
#                   main=titleText );
#         } else if ( as.character( how ) == "Contour" ) {
#             nBins <- nrow( results$PMatrix[[ chr ]] );
#             filled.contour( 1:nBins, 0:maxPloidy, as.matrix( results$PMatrix[[ chr ]] ), 
#                             nlevels=20,
#                             xlab="Bins", ylab="Ploidy",
#                             col=jet.colors( 20 ) );
#             title( titleText );
#         } else if ( as.character( how ) == "Perspective" ) {
#             nBins <- nrow( results$PMatrix[[ chr ]] );
#             persp( 1:nBins, 0:maxPloidy, as.matrix( results$PMatrix[[ chr ]] ),
#                    ticktype="detailed", theta=30, phi=30,
#                    expand=0.5, shade=0.5, col="cyan", ltheta=-30,
#                    xlab="Bins", ylab="Ploidy", zlab="Probability" );  
#             title( titleText );
#         } # if how else
#     } else if ( as.character( what ) == "DeletionScore" ) {
#         plot( results$DeletionScore[[ chr ]], type="l" );
#     } else if ( as.character( what ) == "DuplicationScore" ) {
#         plot( results$DuplicationScore[[ chr ]], type="l" );
#     } # if what else
# } # visualizePosterior
# 
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# aggregateCells <- function( pMatrices ) {
#     nCells <- length( pMatrices );
#     if ( nCells < 1 ) {
#         return( NULL );
#     } # if nCells
#     pMatrix <- pMatrices[[ 1 ]];
#     if ( nCells == 1 ) {
#         return( pMatrix );
#     } # if nCells
#     nChr <- length( pMatrix );
#     maxPloidy <- getMaxPloidy( pMatrix );
#     for ( cell in 2:nCells ) {
#         cat( paste0( "Cell: ", cell, "\n" ) );
#         pMatrix <- lapply( 1:nChr, function( chr ) {
#             cat( paste0( nrow(pMatrix[[chr]]), ", ", ncol(pMatrix[[chr]]), ", ", 
#                          nrow(pMatrices[[cell]][[chr]]), ", ", ncol(pMatrices[[cell]][[chr]]), "\n"))
#             return( as.data.frame( as.matrix( pMatrix[[ chr ]] ) * 
#                                    as.matrix( pMatrices[[ cell ]][[ chr ]] ) ) );
#         } );
#     } # for cell
#     for ( chr in 1:nChr ) {
#         pMatrix[[ chr ]] <- pMatrix[[ chr ]] / rowSums( pMatrix[[ chr ]] );
#     } # for chr
#     result <- list( 
#         "PMatrix"=pMatrix,
#         "MaxLikelihoodPloidy"=getMaxLikelihoodPloidy( pMatrix ) );
#     return( result );
# } # aggregateCells
# 
# ################################################################################
# kullbackLeiblerDivergence <- function( p, q ) {
#     eps <- 1e-2;
#     p <- p + eps;
#     q <- q + eps;
#     return( sum( p * log( p / q ) ) );
# } # kullbackLeiblerDivergence
# 
# ################################################################################
# jensenShannonDivergence <- function( p, q ) {
#     return( 0.5 * ( 
#         kullbackLeiblerDivergence( p, q ) + 
#         kullbackLeiblerDivergence( q, p ) ) );
# } # jensenShannonDivergence
#
# ################################################################################
# getDissimilarity <- function( pMatrixA, pMatrixB, q=0.95 ) {
#     nChr <- length( pMatrixA );
#     d <- rep( NA, nChr );
#     maxPloidy <- getMaxPloidy( pMatrixA );
#     for ( chr in 1:nChr ) {
#         selector <- which( 
#             !is.na( pMatrixA[[ chr ]][ , 1 ] ) &
#             !is.na( pMatrixB[[ chr ]][ , 1 ] ) );
#         dChr <- sapply( selector, function( bin ) {
#             return( jensenShannonDivergence( 
#                 pMatrixA[[ chr ]][ bin, ], 
#                 pMatrixB[[ chr ]][ bin, ] ) );
#         } );
#         # hist( dChr, breaks=50 );
#         d[[ chr ]] <- quantile( dChr, q, na.rm=TRUE );
#     } # for chr
#     return( d );
# } # getDissimilarity
# 
# ################################################################################
# d <- list();
# nClusterSelection <- nrow( clusterSelection );
# for ( cellIndexA in 10:( nClusterSelection - 1 ) ) {
# #for ( cellIndexA in 9 ) {
#     cellA <- clusterSelection$Cell[[ cellIndexA ]];
#     load( paste0( projectFolder, "results_", cellA, ".RData" ) );
#     pMatrixA <- results$PMatrix;
#     for ( cellIndexB in ( cellIndexA + 1 ):nClusterSelection ) {
#     #for ( cellIndexB in 14:29 ) {
#         cellB <- clusterSelection$Cell[[ cellIndexB ]];
#         load( paste0( projectFolder, "results_", cellB, ".RData" ) );
#         pMatrixB <- results$PMatrix;
#         d[[ paste0( cellA, "_", cellB ) ]] <- getDissimilarity( pMatrixA, pMatrixB, q=0.95 );
#         plot( sapply( d, max ), type="b" );
#         cat( paste0( "Completed ", cellA, ",", cellB, "\n" ) );
#     } # for cellB
#     save( list="d", file=paste0( projectFolder, "d.RData" ) );
# } # for cellA
# save( list="d", file=paste0( projectFolder, "d.RData" ) );
#
# ################################################################################
# makeDistanceMatrix <- function( d ) {
#     n <- length( d );
#     nPoints <- 0.5 * ( 1 + sqrt( 8 * n + 1 ) );
#     dMatrix <- matrix( 0, nrow=nPoints, ncol=nPoints );
#     count <- 0;
#     for ( i in 1:( nPoints - 1 ) ) {
#         for ( j in ( i + 1 ):nPoints ) {
#             count <- count + 1;
#             tmp <- sqrt( max( d[[ count ]] ) ); # Shannon-Jensen Divergence is squared distance
#             dMatrix[ i, j ] <- tmp;
#             dMatrix[ j, i ] <- tmp;
#         } # for j
#     } # for i
#     return( dMatrix );
# } # makeDistanceMatrix
# 
# ################################################################################
# makeOldScale <- function( dMatrix ) {
#     # Old scaling factor
#     nPoints <- nrow( dMatrix );
#     oldScale <- 0;
#     for ( i in 1:( nPoints - 1 ) ) {
#         for ( j in ( i + 1 ):nPoints ) {
#             D <- dMatrix[ i, j ];
#             oldScale <- oldScale + D * D;
#         } # for j
#     } # for i
#     oldScale <- oldScale / ( nPoints * nPoints );
#     return( oldScale );
# } # makeOldScale
# 
# ################################################################################
# makeD02 <- function( dMatrix, oldScale ) {
#     # Squared distances from the center of mass
#     nPoints <- nrow( dMatrix );
#     d02 <- rep( 0, nPoints );
#     for ( i in 1:nPoints ) {
#         for ( j in 1:nPoints ) {
#             D <- dMatrix[ i, j ];
#             d02[[ j ]] <- d02[[ j ]] + D * D;
#         } # for j
#     } # for i
#     d02 <- d02 / nPoints - oldScale;
#     return( d02 );
# } # makeD02
# 
# ################################################################################
# makeMetricMatrix <- function( dMatrix, d02 ) {
#     nPoints <- nrow( dMatrix );
#     m <- matrix( 0, nrow=nPoints, ncol=nPoints );
#     for ( i in 1:( nPoints - 1 ) ) {
#         m[ i, i ] <- 2 * d02[[ i ]];
#         for ( j in ( i + 1 ):nPoints ) {
#             D <- dMatrix[ i, j ];
#             tmp <- d02[[ i ]] + d02[[ j ]] - D * D;
#             m[ i, j ] <- tmp;
#             m[ j, i ] <- tmp;
#         } # for j
#     } # for i
#     m <- 0.5 * m;
#     return( m );
# } # makeMetricMatrix
# 
# ################################################################################
# distanceGeometry <- function( d ) {
#     dMatrix <- makeDistanceMatrix( d );
#     # heatmap( dMatrix )
#     oldScale <- makeOldScale( dMatrix );
#     d02 <- makeD02( dMatrix, oldScale );
#     m <- makeMetricMatrix( dMatrix, d02 );
#     # heatmap( m )
#     lambdaU <- eigen( m, symmetric=TRUE, only.values=FALSE );
#     # Eigenalues are already sorted, 
#     # Eigenvectors are already normalized! 
#     x <- sqrt( lambdaU$values[[ 1 ]] ) * lambdaU$vectors[ , 1 ]; # / norm( lambdaU$vectors )
#     y <- sqrt( lambdaU$values[[ 2 ]] ) * lambdaU$vectors[ , 2 ]; # / norm( lambdaU$vectors )
#     # plot( x, rep( 0, length( x ) ) );
#     return( data.frame( x, y ) );
# } # distanceGeometry
# 
################################################################################
################################################################################
################################################################################
################################################################################
