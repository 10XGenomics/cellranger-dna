# differential_variance_ratio.py
#
import numpy as np
import tenkit.stats as tk_stats
import pandas as pd
from statsmodels import robust
import martian
# 
################################################################################
## Auxiliary function
which = lambda lst:list(np.where(lst)[0])

# ################################################################################
##
def get_n_bins(values):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.get_n_bins()')
    ##print('Input: values=')
    ##print(values)
    n_bins = values.shape[0]
    #print('Output: n_bins=%d' % n_bins)
    #print('Leaving differential_variance_ratio.get_n_bins()')
    #print('.' * 80)
    return(n_bins)
# get_n_bins

# ################################################################################
##
def generate_diff_values(values):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.generate_diff_values()')
    #print('Input: values=')
    #print(values)
    diff_values = np.zeros(1)
    if values.shape[0] > 1:
        diff_values = np.hstack((diff_values, np.diff(values)))
    # if values 
    ##print('Output: diff_values=')
    ##print(diff_values)
    #print('Leaving differential_variance_ratio.generate_diff_values()')
    #print('.' * 80)
    return(diff_values)
# generate_diff_values

# ################################################################################
##
def trim(values, trim_level=0.01):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.trim()')
    ##print('Input: values=')
    ##print(values)
    #print('Input: trim_level=')
    #print(trim_level)
    n = get_n_bins(values)
    trimmed_values = np.sort(values)
    n_throw = int(0.5 * trim_level * n)
    #print('n=%d' % n)
    #print('n_throw=%d' % n_throw)
    if (n_throw > 0) and (2 * n_throw < n):
        trimmed_values = trimmed_values[n_throw:-n_throw]
    else:
        trimmed_values = values
    # if n_throw else
    ##print('Output: trimmed_values=')
    ##print(trimmed_values)
    #print('Leaving differential_variance_ratio.trim()')
    #print('.' * 80)
    return(trimmed_values)
# trim

# ################################################################################
## 
def differential_variance_ratio(values, diff_values, trim_level=0.01):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.differential_variance_ratio()')
    ##print('Input: values=')
    ##print(values)
    ##print('Input: diff_values=')
    ##print(diff_values)
    #print('Input: trim_level=%f' % trim_level)
    eps = 1e-6
    n = get_n_bins(values)
    var = np.nanvar(values)
    #print('var=%f' % var)
    # Trim diff_values:
    trimmed_diff_values = trim(diff_values, trim_level)
    ##print('trimmed_diff_values=')
    ##print(trimmed_diff_values)
    diff_var = np.nanvar(trimmed_diff_values)
    #print('diff_var=%f' % diff_var)
    diff_var_ratio = 2.0 * tk_stats.robust_divide(var, diff_var + eps)
    #print('Output: diff_var_ratio=%f' % diff_var_ratio)
    #print('Leaving differential_variance_ratio.differential_variance_ratio()')
    #print('.' * 80)
    return(diff_var_ratio)
# differential_variance_ratio

# ################################################################################
##
def get_z_score(values, diff_values, trim_level=0.01):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.get_z_score()')
    ##print('Input: values')
    ##print(values)
    ##print('Input: diff_values')
    ##print(diff_values)
    ##print('Input: trim_level=')
    ##print(trim_level)
    n = get_n_bins(values)
    assert(n > 0)
    diff_var_ratio = differential_variance_ratio(values, diff_values, trim_level)
    z = (diff_var_ratio - 1) * np.sqrt(n)
    #print('Output: z=%f' % z)
    #print('Leaving differential_variance_ratio.get_z_score()')
    #print('.' * 80)
    return(z)
# get_z_score

# ################################################################################
def assign_blocks(copy_number):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.assign_blocks()')
    ##print('Input: copy_number')
    ##print(copy_number)
    blocks = copy_number
    na_selector = which(np.isnan(blocks))
    blocks = np.diff(np.array([copy_number[0]] + copy_number.tolist()))
    n_blocks = blocks.shape[0]
    selector = which([np.isnan(blocks[i]) and not (i in na_selector)
                     for i in range(n_blocks)])
    blocks[selector] = 1
    blocks = np.abs(blocks)
    blocks[na_selector] = 0
    blocks = np.cumsum(blocks)
    blocks[na_selector] = np.nan
    ##print('Output: blocks=')
    ##print(blocks)
    #print('Leaving differential_variance_ratio.assign_blocks()')
    #print('.' * 80)
    return(blocks)
# assign_blocks
# 
# ################################################################################
##
def initialize_segments(n_segments):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.initialize_segments()')
    #print('Input: n_segments=%d' % n_segments)
    segments = pd.DataFrame({
        'Start':       np.array([np.nan] * n_segments),
        'End':         np.array([np.nan] * n_segments),
        'Length':      np.array([np.nan] * n_segments),
        'Median':      np.array([np.nan] * n_segments),
        'MAD':         np.array([np.nan] * n_segments),
        'Mean':        np.array([np.nan] * n_segments),
        'SD':          np.array([np.nan] * n_segments),
        'Ratio':       np.array([np.nan] * n_segments),
        'F':           np.array([np.nan] * n_segments),
        'CopyNumber':  np.array([np.nan] * n_segments),
        'Homogeneous': np.array([np.nan] * n_segments),
        'Block':       np.array([np.nan] * n_segments)})
    ##print('Output: segments')
    ##print(segments)
    #print('Leaving differential_variance_ratio.initialize_segments()')
    #print('.' * 80)
    return(segments)
# initialize_segments

# ################################################################################
##
def merge_integer_blocks(segments):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.merge_integer_blocks()')
    ##print('Input: segments')
    ##print(segments)
    blocks = np.unique(segments['Block'])
    n_blocks = blocks.shape[0]
    blocks_not_nan = [not np.isnan(blocks[i]) for i in range(n_blocks)]
    blocks = blocks[which(blocks_not_nan)]
    n_homogeneous = blocks.shape[0]
    n_na_blocks = n_blocks - sum(blocks_not_nan)
    n_segments = n_homogeneous + n_na_blocks
    new_segments = initialize_segments(n_segments)
    #
    # Merge neighboring homogeneous blocks: 
    for index in range(n_homogeneous):
        block = blocks[index]
        selector = which(segments['Block'] == block)
        start = segments.loc[selector, 'Start'].min()
        end = segments.loc[selector, 'End'].max()
        new_segments.loc[index, 'Block'] = index
        new_segments.loc[index, 'Start'] = start
        new_segments.loc[index, 'End'] = end
        new_segments.loc[index, 'Length'] = end - start + 1
        new_segments.loc[index, 'Median'] = segments.loc[selector, 'Median'].mean()
        new_segments.loc[index, 'MAD'] = segments.loc[selector, 'MAD'].mean()
        new_segments.loc[index, 'Mean'] = segments.loc[selector, 'Mean'].mean()
        new_segments.loc[index, 'SD'] = segments.loc[selector, 'SD'].mean()
        new_segments.loc[index, 'Ratio'] = segments.loc[selector, 'Ratio'].mean()
        new_segments.loc[index, 'F'] = segments.loc[selector, 'F'].mean()
        new_segments.loc[index, 'CopyNumber'] = segments.loc[min(selector), 'CopyNumber']
        new_segments.loc[index, 'Homogeneous'] = True
    # for index
    #
    # Record nonhomogeneous blocks (to be further subdivided):
    counter = n_homogeneous
    for index in which(np.isnan(segments.loc[:, 'Block'])):
        new_segments.loc[counter, :] = segments.loc[index, :] # Changing things in place?
        counter += 1
    # for index
    new_segments.index = pd.Index(range(new_segments.shape[0]))
    ##print('Output: new_segments')
    ##print(new_segments)
    #print('Leaving differential_variance_ratio.merge_integer_blocks()')
    #print('.' * 80)
    return(new_segments)
# merge_integer_blocks
# 
# ################################################################################
##
def merge_fractional_blocks(segments, values, trim_level=0.01, f_cutoff=3.0, z_cutoff=3.0):
    #print('-' * 80)
    #print('Entering differential_variance_ratio.merge_integer_blocks()')
    segments.sort_values(
        by='Start', axis=0, ascending=True, inplace=True, 
        kind='quicksort', na_position='last')
    #print('Input: sorted segments')
    #print(segments)
    block = 0
    homogeneous_selector = np.where(segments['Homogeneous'] & (segments['Length'] > 1))[0]
    if len(homogeneous_selector) == 0:
        sd_homogeneous = 1e-6
    else:
        sd_homogeneous = np.nanmean(segments.loc[homogeneous_selector, 'SD'])
    # if homogeneous_selector else
    n_segments = segments.shape[0]
    for segment in range(n_segments - 1):
        segments.loc[segment, 'Block'] = np.nan
        if segments.loc[segment, 'Homogeneous']:
            segments.loc[segment, 'Block'] = block
            next_segment = segment + 1
            if segments.loc[next_segment, 'Homogeneous']:
                delta = np.abs(segments.loc[segment, 'Mean'] - segments.loc[next_segment, 'Mean'])
                n_A = segments.loc[segment, 'Length']
                if n_A == 1:
                    sd_A = sd_homogeneous
                else:
                    sd_A = segments.loc[segment, 'SD']
                # if n_A else
                n_B = segments.loc[next_segment, 'Length']
                if n_B == 1:
                    sd_B = sd_homogeneous
                else:
                    sd_B = segments.loc[next_segment, 'SD']
                # if n_B else
                sd = np.sqrt(sd_A * sd_A / n_A + sd_B * sd_B / n_B)
                z = delta / sd
                if z > z_cutoff:
                    block += 1
                # if z
                #if z <= z_cutoff:
                #    segments.loc[next_segment, 'Block'] = block
                #else:
                #    block += 1
                #    segments.loc[next_segment, 'Block'] = block
                ## if z else
            # if next_segment homogeneous
        else:
            block += 1
        # if segment homogeneous else
    # for segment
    #
    blocks = np.unique(segments['Block'])
    n_blocks = blocks.shape[0]
    blocks_not_nan = [not np.isnan(blocks[i]) for i in range(n_blocks)]
    blocks = blocks[which(blocks_not_nan)]
    n_homogeneous = blocks.shape[0]
    n_na_blocks = n_blocks - sum(blocks_not_nan)
    n_segments = n_homogeneous + n_na_blocks
    new_segments = initialize_segments(n_segments)
    #
    # Merge neighboring homogeneous blocks: 
    for index in range(n_homogeneous):
        block = blocks[index]
        selector = which(segments['Block'] == block)
        #print('block=%d' % block)
        #print('selector:')
        #print(selector)
        start = segments.loc[selector, 'Start'].min()
        end = segments.loc[selector, 'End'].max()
        tmp = np.array(values[start:(end + 1)])

        #print('DEBUG: len(tmp)=%d, should be the same as len(diff_tmp)!' % tmp.shape[0])
        #print('tmp:')
        #print(tmp.shape)
        #print(tmp)

        if tmp.shape[0] > 0:
            diff_tmp = generate_diff_values(tmp)
        
            #print('DEBUG: len(diff_tmp)=%d, should be the same as len(tmp)=%d!' % (diff_tmp.shape[0], tmp.shape[0]))
            #print('Median=%f' % float(np.nanmedian(tmp)))
            #print('MAD=%f' % float(robust.mad(tmp)))
            #print('Mean=%f' % np.nanmean(tmp))
            #print('SD=%f' % np.nanstd(tmp))
            #print('Ratio=%f' % differential_variance_ratio(tmp, diff_tmp, trim_level))
            #print('F=%f' % get_z_score(tmp, diff_tmp, trim_level))
        
            new_segments.loc[index, 'Median'] = np.nanmedian(tmp)
            new_segments.loc[index, 'MAD'] = robust.mad(tmp)
            new_segments.loc[index, 'Mean'] = np.nanmean(tmp)
            new_segments.loc[index, 'SD'] = np.nanstd(tmp)
            new_segments.loc[index, 'Ratio'] = differential_variance_ratio(tmp, diff_tmp, trim_level)
            new_segments.loc[index, 'F'] = get_z_score(tmp, diff_tmp, trim_level)
        else:
            print('WARNING: 0-length segment. Start: %d, End: %d' % (int(start), int(end)))
        # if tmp else
        new_segments.loc[index, 'Block'] = index
        new_segments.loc[index, 'Start'] = start
        new_segments.loc[index, 'End'] = end
        new_segments.loc[index, 'Length'] = end - start + 1
        #new_segments.loc[index, 'Median'] = segments.loc[selector, 'Median'].mean()
        #new_segments.loc[index, 'MAD'] = segments.loc[selector, 'MAD'].mean()
        #new_segments.loc[index, 'Mean'] = segments.loc[selector, 'Mean'].mean()
        #new_segments.loc[index, 'SD'] = segments.loc[selector, 'SD'].mean()
        #new_segments.loc[index, 'Ratio'] = segments.loc[selector, 'Ratio'].mean()
        #new_segments.loc[index, 'F'] = segments.loc[selector, 'F'].mean()
        #new_segments.loc[index, 'CopyNumber'] = np.nanmean(segments.loc[selector, 'Mean'])
        new_segments.loc[index, 'CopyNumber'] = new_segments.loc[index, 'Mean']
        new_segments.loc[index, 'Homogeneous'] = new_segments.loc[index, 'F'] < f_cutoff
    # for index
    #
    # Record nonhomogeneous blocks (to be further subdivided):
    counter = n_homogeneous
    for index in which(np.isnan(segments.loc[:, 'Block'])):
        new_segments.loc[counter, :] = segments.loc[index, :] # Changing things in place?
        counter += 1
    # for index
    new_segments.index = pd.Index(range(new_segments.shape[0]))
    #print('Output: new_segments')
    #print(new_segments)
    #print('Leaving differential_variance_ratio.merge_integer_blocks()')
    #print('.' * 80)
    return(new_segments)
# merge_fractional_blocks
# 
# ################################################################################
## 
def remove_NAs(values):
    '''Profile segmentation does not like NaNs'''
    n_bins = values.shape[0]
    na_indices = np.where(np.isnan(values))[0].tolist()
    non_na_indices = set(range(n_bins)).difference(set(na_indices))
    non_na_indices = sorted(list(non_na_indices))
    values = values[non_na_indices]
    return({'NonNAIndices':non_na_indices, 
            'NAIndices':na_indices, 
            'NonNAValues':values})
# remove_NAs
#
# ################################################################################
## 
def reintroduce_NA_bins(segments, non_na_indices):
    '''Purpose: reintroduce NA bins and shift CNV breakpoints to their correct locations.'''
    segments['Start'] = np.array(non_na_indices)[segments['Start']]
    segments['End'] = np.array(non_na_indices)[segments['End']]
    segments['Length'] = (segments['End'] - segments['Start'] + 1).tolist()
    segments['Start'] = segments['Start'].tolist()
    segments['End'] = segments['End'].tolist()
    return(segments)
# reintroduce_NA_bins
#
# ################################################################################
## 
def partition_profile(
    values, segment_length=18, f_cutoff=3.0, trim_level=0.01, min_length=1, 
    enable_fractional_elevations=False, z_cutoff=3.0):
    ##print('-' * 80)
    #print('Entering differential_variance_ratio.partition_profile()')
    #
    clean_results = remove_NAs(values)
    non_na_indices = clean_results['NonNAIndices']
    na_indices = clean_results['NAIndices']
    values = clean_results['NonNAValues']
    #
    diff_values = generate_diff_values(values)
    n_bins = get_n_bins(values)
    n_segments = int(np.ceil(tk_stats.robust_divide(n_bins, float(segment_length))))
    diff_indices = np.argsort(np.abs(diff_values))
    breakpoints = diff_indices[-n_segments:]
 
    ##print('diff_values:')
    ##print(diff_values)
    #print('n_segments: %d' % n_segments)
    ##print('diff_indices:')
    ##print(diff_indices)
    ##print('breakpoints:')
    ##print(breakpoints)

    start_points = np.unique(np.concatenate((np.zeros(1), breakpoints)))
    start_points = np.unique(start_points)
    end_points = np.unique(np.concatenate((start_points[1:] - 1, np.array([n_bins - 1]))))
    n_segments = start_points.shape[0]
    segments = initialize_segments(n_segments)
    if n_segments == 0:
        martian.log_warn('WARNING: n_segments = 0, returning empty DataFrame')
        return(segments)
    # if n_segments
    #assert(n_segments == end_points.shape[0])
    # Debug:
    if n_segments != end_points.shape[0]:
        martian.log_warn('WARNING: n_segments=%d, length(end_points)=%d' % (n_segments, end_points.shape[0]))
        return(segments)
    # if n_segments
    #
    segments['Start'] = start_points
    segments['End'] = end_points
    segments['Length'] = end_points - start_points + 1
    indices = np.arange(n_segments)
    for index in indices:
        start = segments.loc[index, 'Start']
        end = segments.loc[index, 'End'] + 1
        tmp = np.array(values[start:end])

        #print('DEBUG: len(tmp)=%d, should be the same as len(diff_tmp)!' % tmp.shape[0])

        if tmp.shape[0] > 0:
            diff_tmp = generate_diff_values(tmp)

            #print('DEBUG: len(diff_tmp)=%d, should be the same as len(tmp)=%d!' % (diff_tmp.shape[0], tmp.shape[0]))

            segments.loc[index, 'Median'] = np.nanmedian(tmp)
            segments.loc[index, 'MAD'] = robust.mad(tmp)
            segments.loc[index, 'Mean'] = np.nanmean(tmp)
            segments.loc[index, 'SD'] = np.nanstd(tmp)
            segments.loc[index, 'Ratio'] = differential_variance_ratio(tmp, diff_tmp, trim_level)
            segments.loc[index, 'F'] = get_z_score(tmp, diff_tmp, trim_level)
        else:
            martian.log_warn('WARNING: 0-length segment. Start: %d, End: %d' %
                  (int(start), int(end)))
        # if tmp else
    # for index    
    #
    # Identify homogeneous/inhomogeneous segments
    segments['Homogeneous'] = (segments['F'] < f_cutoff)
    homogeneous_selector = which(segments['Homogeneous'])
    #
    # Compare homogeneous segments
    if enable_fractional_elevations:
        segments.loc[index, 'CopyNumber'] = segments.loc[index, 'Mean']
        #
        # Merge neighboring homogeneous segments
        segments = merge_fractional_blocks(segments, values, trim_level, f_cutoff, z_cutoff)
    else:
        segments.loc[homogeneous_selector, 'CopyNumber'] = np.round(
            segments.loc[homogeneous_selector, 'Median'])
        #segments.loc[homogeneous_selector, 'CopyNumber'] = np.round(
        #    segments.loc[homogeneous_selector, 'Median'] * 2) / 2
        segments['Block'] = assign_blocks(segments['CopyNumber'])
        #
        # Merge neighboring homogeneous segments
        segments = merge_integer_blocks(segments)
    # if enable_fractional_elevations else
    #
    # Recursion: partition inhomogeneous segments
    if segment_length > min_length:
        try:
            new_length = max(min_length, segment_length / 5)
            unfinished_selector = which(np.isnan( segments['Block']))
            while len(unfinished_selector) > 0:
                index = unfinished_selector[0]
                try:
                    selector = np.arange(int(segments.loc[index, 'Start']), 
                                         int(segments.loc[index, 'End'] + 1))
                    new_values = values[selector]
                    new_segments = partition_profile(
                        new_values, new_length,
                        f_cutoff, trim_level)
                    prefix = segments.loc[index, 'Start']
                    new_segments['Start'] += prefix
                    new_segments['End'] += prefix
                    segments.drop(index, inplace=True)
                    segments = pd.concat([segments, new_segments])
                    segments.index = pd.Index(range(segments.shape[0]))
                    segments['Block'] = assign_blocks(segments['CopyNumber'])
                    unfinished_selector = which(np.isnan( segments['Block']))
                except Exception as error:
                    martian.log_warn('[partition_profile:(for index)] caught this error: ' 
                          + repr(error))
                # try except
            # for index
        except Exception as error:
            martian.log_warn('[partition_profile:(if segment_length)] caught this error: ' 
                  + repr(error))
        # try except
    # if segment_length
    #
    indices = np.argsort(segments['Start'])
    segments.index = pd.Index(range(segments.shape[0]))
    segments = segments.loc[indices, :].dropna(how='all')
    segments.index = pd.Index(range(segments.shape[0]))
    unfinished_selector = which(np.isnan(segments['Block']))
    for index in unfinished_selector:
        selector = np.arange(int(segments.loc[index, 'Start']), 
                             int(segments.loc[index, 'End'] + 1))
        if enable_fractional_elevations:
            segments.loc[index, 'CopyNumber'] = segments.loc[index, 'Mean']
        else:
            segments.loc[index, 'CopyNumber'] = np.round(
                np.nanmean(values[selector]))
            #segments.loc[index, 'CopyNumber'] = np.round(
            #    np.nanmean(values[selector] * 2)) / 2
        # if !enable_fractional_elevations
    # for index
    if enable_fractional_elevations:
        segments = merge_fractional_blocks(segments, values, trim_level, f_cutoff, z_cutoff);
    else:
        segments['Block'] = assign_blocks(segments['CopyNumber']);
        segments = merge_integer_blocks(segments);
    # if enable_fractional_elevations else
    #
    indices = np.argsort(segments['Start'])
    segments = segments.loc[indices, :]
    #indices = segments.index
    #
    segments = reintroduce_NA_bins(segments, non_na_indices)
    #
    #print('Leaving differential_variance_ratio.partition_profile()')
    #print('.' * 80)
    return(segments);
# partition_profile


