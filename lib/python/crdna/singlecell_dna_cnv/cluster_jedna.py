# cluster_jedna.py

#import math
import numpy as np
#import pandas as pd
from differential_variance_ratio import generate_diff_values 
from differential_variance_ratio import differential_variance_ratio 
from differential_variance_ratio import get_z_score 

################################################################################
## Auxiliary function
which = lambda lst:list(np.where(lst)[0])

################################################################################
## Merge bins: from 20 kb to 100 kb resolution
def merge_bins(profiles, n_merge):
    print('DEBUG: Entering cluster_jedna.merge_bins()')

    #print('DEBUG input profiles[0][0]:') 
    #print(profiles[0][0].tolist())

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
        # n_bins = len(profiles[chrom][0])
        n_bins = profiles[chrom].shape[1]
        print('chrom=%d, n_bins=%d' % (chrom, n_bins))
        if n_bins < n_merge:
            # TODO: warning
            pass
        # if n_bins
        #n_merged_bins = int(math.ceil(n_bins / float(n_merge)))
        tmp = np.cumsum(profiles[chrom], axis=1, dtype=float) 
        tmp[:, n_merge:] = tmp[:, n_merge:] - tmp[:, :-n_merge]
        merged_profiles[chrom] = tmp[:, n_merge - 1::n_merge] # This throws away any leftover bins!
    # for chrom

    #print('DEBUG: merged_profiles[0][0]:')
    #print(merged_profiles[0][0])

    print('DEBUG: Leaving cluster_jedna.merge_bins()')
    return(merged_profiles)
# merge_bins

################################################################################
##
def merge_mask(mask, n_merge):
    print('DEBUG: Entering cluster_jedna.merge_mask()')
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
        merged_mask[chrom] = np.array([
            bool(np.prod(mask[chrom][i:(i + n_merge)]))
            for i in range(0, n_bins + 1, n_merge)
            ][:-1])
    # for chrom
    print('DEBUG: Leaving cluster_jedna.merge_mask()')
    return(merged_mask)
# merge_mask

################################################################################
##
def get_total_read_counts(profiles, mask):
    print('DEBUG: Entering cluster_jedna.get_total_read_counts()')
    total_read_counts = [None]
    n_chrom = len(profiles)  
    if n_chrom < 1:  
        # TODO: warning: no chromosomes!  
        return(total_read_counts)  
    # if n_chrom  
    n_cells = profiles[0].shape[0]
    total_read_counts = [0] * n_cells
    if n_cells < 1:  
        # TODO: warning: no cells!  
        return(total_read_counts)  
    # if n_cells  
    for chrom in range(n_chrom):
        chr_read_counts = np.nansum(np.transpose(profiles[chrom])[mask[chrom]], axis=0)
        total_read_counts = [total_read_counts[i] + chr_read_counts[i] 
            for i in range(n_cells)]
    # for chrom
    print('DEBUG: Leaving cluster_jedna.get_total_read_counts()')
    return(total_read_counts)  
# get_total_read_counts

################################################################################
##
def get_total_n_bins(mask):
    print('DEBUG: Entering cluster_jedna.get_total_n_bins()')
    total_n_bins = 0
    n_chrom = len(mask)  
    if n_chrom < 1:
        # TODO: warning: no chromosomes!  
        return(total_n_bins)
    # if n_chrom  
    for chrom in range(n_chrom):
        total_n_bins += sum(mask[chrom])
    # for chrom
    print('DEBUG: Leaving cluster_jedna.get_total_n_bins()')
    return(total_n_bins)
# get_total_n_bins

################################################################################
## Scale with respect to total reads
def scale(profiles, mask):
    print('DEBUG: Entering cluster_jedna.scale()')
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

    #print('DEBUG profiles[0][0] before get_total_read_counts:')
    #print(profiles[0][0].tolist())

    total_read_counts = get_total_read_counts(profiles, mask)

    #print('DEBUG profiles[0][0] after get_total_read_counts:')
    #print(profiles[0][0].tolist())

    scaled_profiles = [None] * n_chrom
    mean_reads_per_bin = np.array([total_read_counts[i] / float(total_n_bins) 
        for i in range(n_cells)])
    print('DEBUG: mean_reads_per_bin')
    print(mean_reads_per_bin)
    print('total_n_bins=%d' % total_n_bins)

    for chrom in range(n_chrom):
        #print('DEBUG profiles[%d]' % chrom)
        #print(profiles[chrom].tolist())
        scaled_profiles[chrom] = profiles[chrom] / mean_reads_per_bin[:, None]
        #print('DEBUG scaled_profiles[%d]' % chrom)
        #print(scaled_profiles[chrom].tolist())
    # for chr
    print('DEBUG: Leaving cluster_jedna.scale()')
    return(scaled_profiles)
# scale

################################################################################
## Preprocess profiles
def process_profiles(raw_profiles, mask, n_merge):
    print('DEBUG: Entering cluster_jedna.process_profiles()')
    # Merge bins: from 20 kb to 100 kb resolution

    #print('DEBUG raw_profiles[0][0]:')
    #print(raw_profiles[0][0].tolist())
    #print('DEBUG mask[0]:')
    #print(mask[0].tolist())

    clean_profiles, clean_mask = remove_NAs(raw_profiles, mask)
    merged_profiles = merge_bins(clean_profiles, n_merge)

    #print('DEBUG profiles[0][0] after merging:')
    #print(merged_profiles[0][0].tolist())

    merged_mask = merge_mask(clean_mask, n_merge)

    #print('DEBUG mask[0] after merging:')
    #print(merged_mask[0].tolist())

    # Scale with respect to total reads
    profiles = scale(merged_profiles, merged_mask)

    print('DEBUG: Leaving cluster_jedna.process_profiles()')
    return(profiles, merged_mask)
# process_profiles

################################################################################
def evaluate_distances(profiles, mask, method='dvr_per_chrom'):
    print('DEBUG: Entering cluster_jedna.evaluate_distances()')
    if method == 'dvr_per_chrom':
        return(evaluate_dvr_per_chrom(profiles, mask, weight=100.0))
    elif method == 'dvr_per_genome':
        return(evaluate_dvr_per_genome(profiles, mask, 
            trim_level=0.01, tile_size=18, weight=1.0, dvr_z_cutoff=3.0, chrom_diff_cutoff=3.0))
    else:
        print('Unknown method: %s' % method)
        return({'DVR':[], 'ChromDelta':[], 'Z':[], 'Penalty':[]})
    # end if method else
    print('DEBUG: Leaving cluster_jedna.evaluate_distances()')
# evaluate_distances

################################################################################
def distance(cell_i, cell_j, profiles, mask, n_chrom, trim_level=0.01, weight=100):
    print('DEBUG: Entering cluster_jedna.distance()')
    dvr = [None] * n_chrom
    chrom_delta = [None] * n_chrom
    z = [None] * n_chrom
    penalty = [None] * n_chrom
    for chrom in range(n_chrom):

        #print('DEBUG profiles chrom %d cell %d:' %(chrom, cell_i))
        #print(profiles[chrom][cell_i, :].tolist())
        #print('DEBUG profiles chrom %d cell %d:' %(chrom, cell_j))
        #print(profiles[chrom][cell_j, :].tolist())

        delta = np.subtract(
            profiles[chrom][cell_i, :],
            profiles[chrom][cell_j, :])
        delta = delta[which(mask[chrom])]

        #print('DEBUG delta before cleanup:')
        #print(delta.shape)
        #print(delta.tolist())

        delta = delta[np.where(np.logical_not(np.isnan(delta)))[0]]

        #print('DEBUG delta after cleanup:')
        #print(delta.shape)
        #print(delta.tolist())

        diff_delta = generate_diff_values(delta)
        dvr[chrom] = differential_variance_ratio(delta, diff_delta, trim_level)
        tmp_z = get_z_score(delta, diff_delta, trim_level)
        tmp_chrom_delta = np.abs(np.nanmean(delta))
        z[chrom] = tmp_z
        chrom_delta[chrom] = tmp_chrom_delta
        penalty[chrom] = tmp_z + weight * tmp_chrom_delta

        print('DEBUG: z_score=%f' % tmp_z)
        print('DEBUG: mean(delta)=%f' % tmp_chrom_delta)
        print('DEBUG: chrom=%d' % chrom)
        print('DEBUG: dvr=%f' % dvr[chrom])
        print('DEBUG: final z=%f' % z[chrom])

    # for chr
    print('DEBUG: Leaving cluster_jedna.distance()')
    return({'DVR':dvr, 'ChromDelta':chrom_delta, 'Z':z, 'Penalty':penalty})
# distance

################################################################################
def get_n_tiles(sequence, tile_size):
    print('DEBUG: Entering cluster_jedna.get_n_tiles()')
    n_tiles = 0
    n_chrom = len(sequence)
    print('DEBUG: n_chrom=%d' % n_chrom)
    for chrom in range(n_chrom):
        n_bins = len(sequence[chrom])
        n_chrom_tiles = n_bins // tile_size
        n_chrom_tiles *= 2
        n_tiles += n_chrom_tiles
        print('DEBUG: n_chrom_tiles=%d' % n_chrom_tiles)
    # for chrom
    print('DEBUG: n_tiles=%d' % n_tiles)
    print('DEBUG: Leaving cluster_jedna.get_n_tiles()')
    return(n_tiles)
# get_n_tiles

################################################################################
def remove_NAs(profiles, mask):
    n_chrom = len(mask)
    clean_mask = [None] * n_chrom
    clean_profiles = [None] * n_chrom
    for chrom in range(n_chrom):
        n_bins = len(mask[chrom])
        profile_selector = [not np.any(np.isnan(profiles[chrom][:, bin_index]))
            for bin_index in range(n_bins)]
        selector = [mask[chrom][bin_index] and profile_selector[bin_index]
            for bin_index in range(n_bins)]
        selector = np.where(selector)[0].tolist()
        clean_mask[chrom] = np.array(mask[chrom])[selector]
        clean_profiles[chrom] = profiles[chrom][:, selector]
    # for chrom
    return((clean_profiles, clean_mask))
# remove_NAs

################################################################################
# Evaluate distances
def evaluate_dvr_per_genome(profiles, mask, 
    trim_level=0.01, tile_size=18, weight=1.0, dvr_z_cutoff=3.0, chrom_diff_cutoff=3.0, min_tile_length=12):
    print('DEBUG: Entering cluster_jedna.evaluate_dvr_per_genome()')
    d = None
    n_chrom = len(profiles)
    if n_chrom < 1:
        # TODO: warning - no chromosomes
        return(d)
    # if n_chrom
    if len(mask) != n_chrom:
        # TODO: warning - mask does not match profiles
        return(d)
    # if mask
    n_cells = profiles[0].shape[0]
    if n_cells < 1:
        # TODO: warning: no cells
        return(d)
    # if n_celss
    #
    n_pairs = int(0.5 * n_cells * ( n_cells - 1 ))
    d = {'I':([None] * n_pairs),
         'J':([None] * n_pairs),
         'DVR':([None] * n_pairs),
         'ChromDelta':([None] * n_pairs),
         'Z':([None] * n_pairs),
         'Penalty':([None] * n_pairs)}
    count = 0
    profiles, mask = remove_NAs(profiles, mask)
    n_tiles = get_n_tiles(mask, tile_size)
    #
    print('DEBUG: n_tiles=%d' % n_tiles)
    #
    sqrt_tile_size = np.sqrt(tile_size)
    for cell_i in range(0, n_cells - 1):
        for cell_j in range(cell_i + 1, n_cells):
            dvr_z_values = np.zeros(n_tiles)
            chrom_diff_values = np.zeros(n_tiles)
            tile_count = 0
            for chrom in range(n_chrom):
                n_bins = len(mask[chrom])
                n_chrom_tiles = (n_bins // tile_size) * 2
                starts = np.arange(0, n_bins, tile_size // 2)[:n_chrom_tiles]
                #
                print('DEBUG: chrom=%d, n_bins=%d, n_chrom_tiles=%d' % (chrom, n_bins, n_chrom_tiles))
                #print('DEBUG: starts:')
                #print(starts)
                #
                for start in starts:
                    tile_end = min(start + tile_size, n_bins - 1)
                    if tile_end > (start + min_tile_length):
                        tile_i = profiles[chrom][cell_i, start:tile_end]
                        tile_j = profiles[chrom][cell_j, start:tile_end]
                        delta = np.subtract(tile_i, tile_j)
                        diff_delta = generate_diff_values(delta)
                        # This is not used downstream, commenting it out:
                        #dvr = differential_variance_ratio(delta, diff_delta, trim_level)
                        dvr_z = get_z_score(delta, diff_delta, trim_level)
                        #
                        print('DEBUG n_tiles=%d, tile_count=%d, dvr_z=%f' % (n_tiles, tile_count, dvr_z))
                        #
                        dvr_z_values[tile_count] = dvr_z
                        #
                        #print('DEBUG dvr_z_values:')
                        #print(dvr_z_values)
                        #
                        chrom_diff = abs(np.nanmean(delta)) * sqrt_tile_size / np.sqrt(np.nanvar(tile_i) + np.nanvar(tile_j))
                        # 
                        print('DEBUG chrom_diff=%f' % chrom_diff)
                        #
                        chrom_diff_values[tile_count] = chrom_diff
                        #
                        print('DEBUG: dvr_z=%f, chrom_diff=%f' % (dvr_z, chrom_diff))
                        #
                        tile_count += 1
                    # if tile_end
                # for start
            # for chrom
            print('DEBUG: tile_count=%d, len(dvr_z_values)=%d' % (tile_count, len(dvr_z_values)))
            dvr_z_values = dvr_z_values[:tile_count] # dvr_z_values[np.where(~np.isnan(dvr_z_values))[0]]
            chrom_diff_values = chrom_diff_values[:tile_count] # chrom_diff_values[np.where(~np.isnan(chrom_diff_values))[0]]
            n_dvr_z = dvr_z_values.shape[0]
            n_chrom_diff = chrom_diff_values.shape[0]
            n_dvr_above_cutoff = len(np.where(dvr_z_values > dvr_z_cutoff)[0])
            n_chrom_diff_above_cutoff = len(np.where(chrom_diff_values > chrom_diff_cutoff)[0])
            d['I'][count] = cell_i
            d['J'][count] = cell_j
            dvr_penalty = n_dvr_above_cutoff / float(n_dvr_z)
            chrom_diff_penalty = n_chrom_diff_above_cutoff / float(n_chrom_diff)
            combined_penalty = len(np.where(np.logical_or(dvr_z_values > dvr_z_cutoff, chrom_diff_values > chrom_diff_cutoff))[0]) / float(n_chrom_diff)
            d['DVR'][count] = dvr_penalty
            d['ChromDelta'][count] = chrom_diff_penalty
            #d['Penalty'][count] = dvr_penalty + weight * chrom_diff_penalty
            d['Penalty'][count] = combined_penalty
            count += 1
        # for cell_j
    # for cell_i
    print('DEBUG: Leaving cluster_jedna.evaluate_dvr_per_genome()')
    return(d)
# evaluate_dvr_per_genome

################################################################################
## Evaluate distances
def evaluate_dvr_per_chrom(profiles, mask, weight=100):
    print('DEBUG: Entering cluster_jedna.evaluate_dvr_per_chrom()')
    d = None
    n_chrom = len(profiles)
    if n_chrom < 1:
        # TODO: warning - no chromosomes
        return(d)
    # if n_chrom
    if len(mask) != n_chrom:
        # TODO: warning - mask does not match profiles
        return(d)
    # if mask
    n_cells = profiles[0].shape[0]
    if n_cells < 1:
        # TODO: warning: no cells
        return(d)
    # if n_celss
    ##...........................................................................
    #def distance(cell_i, cell_j, trim_level=0.01):
    #    dvr = [None] * n_chrom
    #    chrom_delta = [None] * n_chrom
    #    z = [None] * n_chrom
    #    penalty = [None] * n_chrom
    #    for chrom in range(n_chrom):
    #        delta = np.subtract(
    #            profiles[chrom][cell_i, :],
    #            profiles[chrom][cell_j, :])
    #        delta = delta[which(mask[chrom])]
    #        diff_delta = generate_diff_values(delta)
    #        dvr[chrom] = differential_variance_ratio(delta, diff_delta, trim_level)
    #        tmp_z = get_z_score(delta, diff_delta, trim_level)
    #        tmp_chrom_delta = np.abs(np.nanmean(delta))
    #        z[chrom] = tmp_z
    #        chrom_delta[chrom] = tmp_chrom_delta
    #        penalty[chrom] = tmp_z + weight * tmp_chrom_delta
    #        print('DEBUG: z_score=%f' % tmp_z)
    #        print('DEBUG: mean(delta)=%f' % tmp_chrom_delta)
    #        print('DEBUG: chrom=%d' % chrom)
    #        print('DEBUG: dvr=%f' % dvr[chrom])
    #        print('DEBUG: final z=%f' % z[chrom])
    #    # for chr
    #    return({'DVR':dvr, 'ChromDelta':chrom_delta, 'Z':z, 'Penalty':penalty})
    ## distance
    #
    n_pairs = int(0.5 * n_cells * ( n_cells - 1 ))
    d = {'I':([None] * n_pairs),
         'J':([None] * n_pairs),
         'DVR':([None] * n_pairs),
         'ChromDelta':([None] * n_pairs),
         'Z':([None] * n_pairs),
         'Penalty':([None] * n_pairs)}
    count = 0
    for cell_i in range(0, n_cells - 1):
        for cell_j in range(cell_i + 1, n_cells):
            d['I'][count] = cell_i
            d['J'][count] = cell_j
            #tmp = distance(cell_i, cell_j)
            tmp = distance(cell_i, cell_j, profiles, mask, n_chrom)
            #print('DEBUG: tmp.keys()')
            #print(tmp)
            #print(tmp.keys())
            d['DVR'][count] = tmp['DVR']
            d['ChromDelta'][count] = tmp['ChromDelta']
            #print('Before getting Z')
            d['Z'][count] = tmp['Z']
            d['Penalty'][count] = tmp['Penalty']
            #print('After getting Z')
            count += 1
        # for cell_j
    # for cell_i
    print('DEBUG: Leaving cluster_jedna.evaluate_dvr_per_chrom()')
    return(d)
# evaluate_dvr_per_chrom

################################################################################
def get_z(d, chrom):
    print('DEBUG: Entering cluster_jedna.get_z()')
    n_pairs = len(d['Penalty'])
    z = {'I': d['I'], 
         'J': d['J'], 
         'DVR': [d['DVR'][index][chrom] for index in range(n_pairs)],
         'ChromDelta': [d['ChromDelta'][index][chrom] for index in range(n_pairs)],
         'Z': [d['Z'][index][chrom] for index in range(n_pairs)],
         'Penalty': [d['Penalty'][index][chrom] for index in range(n_pairs)]}
    return(z)
    print('DEBUG: Leaving cluster_jedna.get_z()')
# get_z

################################################################################
## Greedy association
def cluster_greedily(z, score_cutoff=15):
    print('DEBUG: Entering cluster_jedna.cluster_greedily()')
    #
    print('DEBUG: input z[I]:')
    print(z['I'])
    print('DEBUG: input z[J]:')
    print(z['J'])
    print('DEBUG: input z[ChromDelta]:')
    print(z['ChromDelta'])
    print('DEBUG: input z[Z]:')
    print(z['Z'])
    print('DEBUG: input z[Penalty]:')
    print(z['Penalty'])
    #
    #epsilon = 1e6
    profile_ids = range(z.get('J')[-1] + 1)
    clusters = [] 
    seed = np.argmin(z.get('Penalty')) 
    cluster = [z.get('I')[seed], z.get('J')[seed]] 
    profile_ids.remove(z.get('I')[seed])
    profile_ids.remove(z.get('J')[seed])
    print('=' * 80)
    #print('New cluster: [%d-$d].\nClustered: %d. Unclustered: %d. Score: [%f]' % 
    #    (cluster[0], cluster[1], len(cluster), len(profile_ids), z.get('Penalty')[seed])) 
    print('New cluster: [TODO].\nClustered: %d. Unclustered: %d. Score: [%f]' % 
        (len(cluster), len(profile_ids), z.get('Penalty')[seed])) 
    min_d = np.percentile(z.get('Penalty'), score_cutoff)
    #
    print('DEBUG: min_d=%f' % min_d)
    #
    while len(profile_ids) > 0:
         #
         # Proximity within cluster
         cluster_size = len(cluster)
         #
         print('DEBUG 1')
         print('DEBUG: cluster_size=%d' % cluster_size)
         #
         n_pairs = cluster_size * (cluster_size - 1) / 2
         #
         print('DEBUG 2')
         print('n_pairs=%d' % n_pairs)
         #
         tmp = [None] * n_pairs
         #
         print('DEBUG 3')
         #
         count = 0
         for k in range(0, cluster_size - 1):
             for s in range(k + 1, cluster_size):
                 #
                 print('DEBUG 4')
                 print('DEBUG: k=%d, s=%d, cluster[k]=%d, cluster[s]=%d' % (k, s, cluster[k], cluster[s]))
                 #
                 auxiliary = map(lambda x, y:((x == cluster[k] and y == cluster[s]) or
                                              (x == cluster[s] and y == cluster[k])), 
                                              z.get('I'), z.get('J'))
                 #
                 print('DEBUG 5')
                 #
                 selector = which(auxiliary)
                 #
                 print('DEBUG 6')
                 print('DEBUG selector:')
                 print(selector)
                 #
                 tmp[count] = z.get('Penalty')[selector[0]]
                 #
                 print('DEBUG 7')
                 print('count=%d, tmp[count]=%f' % (count, tmp[count]))
                 #
                 count += 1
             # for s
         # for k
         #d_within = np.mean(tmp)
         #s_within = np.std(tmp)
         #
         # Proximity between cluster members and outside indices
         #
         print('DEBUG 8')
         #
         d_between = [None] * len(profile_ids)
         #
         print('DEBUG 9')
         #
         #s_between = [None] * len(profile_ids)
         index = 0
         for out_id in profile_ids:
             #
             print('DEBUG 10')
             print('out_id=%d' % out_id)
             #
             tmp = [None] * cluster_size
             #
             print('DEBUG 11')
             #
             count = 0
             for cluster_id in cluster:
                 #
                 print('DEBUG 12')
                 print('cluster_id=%d' % cluster_id)
                 #
                 auxiliary = map(lambda x, y:((x == out_id and y == cluster_id) or
                                              (x == cluster_id and y == out_id)),
                                              z.get('I'), z.get('J'))
                 #
                 print('DEBUG 13')
                 #
                 selector = which(auxiliary)
                 #
                 print('DEBUG 14')
                 #print('auxiliary:')
                 #print(auxiliary)
                 print('selector:')
                 print(selector)
                 print('tmp:')
                 print(tmp)
                 print('count=%d' % count)
                 #print('I=%d, J=%d' % (I, J))
                 #
                 if len(selector) > 0: # THIS IS A HACK - need to revisit the logic
                     print('len(selector)>0')
                     tmp[count] = z.get('Penalty')[selector[0]]
                 else:
                     print('len(selector)=0')
                     tmp[count] = 1e3 # THIS IS A HACK - need to revisit the logic
                 # if selector
                 #
                 print('DEBUG 15')
                 #
                 count += 1
             # for cluster_id
             #
             print('DEBUG 16')
             #
             d_between[index] = np.mean(tmp)
             #
             print('DEBUG 17')
             print('d_between[index]=%f' % d_between[index])
             #
             #s_between[index] = np.std(tmp)
             index += 1
         # for out_id
         #
         # Find the best candidates
         auxiliary = map(lambda x:x <= min_d, d_between)
         #
         print('DEBUG 18')
         print('auxiliary:')
         print(auxiliary)
         #
         k_best = which(auxiliary)
         #
         print('DEBUG 19')
         print('k_best:')
         print(k_best)
         #
         # Decision:
         if len(k_best) > 0:
             #
             # If yes, augment the cluster
             new_ids = [profile_ids[i] for i in k_best]
             #
             print('DEBUG 20')
             print('new_ids:')
             print(new_ids)
             #
             cluster.extend(new_ids)
             #
             print('DEBUG 21')
             print('Extended cluster:')
             print(cluster)
             #
             [profile_ids.remove(out_id) for out_id in new_ids]
             #
             print('DEBUG 22')
             #
             print('Growing cluster: [%s].\nClustered: %d. Unclustered: %d. Score: [%f, %f]' % 
                  (','.join(map(str, cluster)) + '\n', len(cluster), len(profile_ids), min(d_between), min_d))
         # if k_best
         #
         # Close this cluster
         clusters.append(cluster)
         #
         # Seed another cluster
         auxiliary = map(lambda x, y:(x in profile_ids) and (y in profile_ids),
                                     z.get('I'), z.get('J'))
         selector = which(auxiliary)
         if len(selector) > 0:
             seed = np.argmin([z.get('Penalty')[i] for i in selector])
             seed_index = selector[seed]
             min_z = z.get('Penalty')[seed_index]
         else:
             if len(profile_ids) == 1:
                 clusters.append([profile_ids[0]])
                 profile_ids = []
                 continue
             else:
                 print('DEBUG profile_ids:')
                 print(profile_ids)
                 continue
             # if len
         # if selector else
         print('=' * 80)
         if min_z > min_d:
             # Further grouping is too loose - treat all remaining barcodes as single-cell clusters
             for singlet_id in profile_ids:
                 clusters.append([singlet_id])
                 print('Single cell cluster: [%d].' % singlet_id) # TODO: display cluster
             # for singlet_id
             profile_ids = []
         else:
             print('New cluster: [%d-%d]. Clustered: %s. Unclustered: %d. Score: [%f]' % 
                 (z.get('I')[seed_index], 
                  z.get('J')[seed_index],
                  'TODO', len(profile_ids), min_z)) # TODO: display cluster, #clustered cells
             cluster = [z.get('I')[seed_index], z.get('J')[seed_index]]
             profile_ids.remove(z.get('I')[seed_index])
             profile_ids.remove(z.get('J')[seed_index])
         # if z else
    # while profile_ids
    #
    print('DEBUG: Leaving cluster_jedna.cluster_greedily()')
    return(clusters)
# cluster_greedily

################################################################################
# findOutliers <- function( cluster, clusterProfiles, totalReadsPerBin, hg19NGapsMask, pValueCutoff=0.001 ) {
#     outliers <- list();
#     meanClusterProfiles <- list();
#     nChr <- 24;
#     outCounter <- 0;
#     pBonferroni <- pValueCutoff / nChr;
#     for ( chr in 1:nChr ) {
#         #
#         # Evaluate mean cluster profile
#         meanClusterProfile <- 0 * clusterProfiles[[ 1 ]][[ chr ]];
#         for ( index in seq_along( cluster ) ) {
#             meanClusterProfile <- meanClusterProfile + clusterProfiles[[ index ]][[ chr ]];
#         } # for index
#         meanClusterProfile <- meanClusterProfile / totalReadsPerBin;
#         n <- min( length( genomeProfile[[ chr ]] ), length( hg19NGapsMask[[ chr ]] ) );
#         meanClusterProfiles[[ chr ]] <- meanClusterProfile[ 1:n ] * hg19NGapsMask[[ chr ]][ 1:n ];
#         #
#         # Evaluate distances from the mean cluster profile
#         delta2 <- rep( NA, length( cluster ) );
#         count <- 0;
#         for( index in seq_along( cluster ) ) {
#             count <- count + 1;
#             tmp <- scaleProfile( clusterProfiles[[ index ]] )[[ chr ]] * hg19NGapsMask[[ chr ]][ 1:n ];
#             delta2[[ count ]] <- sum( ( meanClusterProfiles[[ chr ]] - tmp )^2, na.rm=TRUE );
#         } # for index
#         #
#         # Identify outliers
#         tValue <- scale( log10( delta2 ) );
#         pValue <- 1 - pt( tValue, df=length( cluster ) );
#         outsiders <- which( pValue < pBonferroni );
#         if ( length( outsiders ) > 0 ) {
#             outCounter <- outCounter + 1;
#             outliers[[ outCounter ]] <- list( 'Chr'=chr, 'Outliers'=outsiders );
#         } # if outsiders
#     } # for chr
#     #
#     # Remove outliers from the cluster
#     # TODO
#     #
#     result <- list( 'Cluster'=cluster, 'Outliers'=outliers );
#     return( result );
# } # findOutliers
# 
# Cluster cleanup
def cleanup_clusters(clusters):
    print('DEBUG: Entering cluster_jedna.cleanup_clusters()')
    # TODO: Implement
    print('DEBUG: Leaving cluster_jedna.cleanup_clusters()')
    return(clusters)
# cleanup_clusters

################################################################################
## Integrate all pieces
def cluster(raw_profiles, mask, n_merge, score_cutoff=5):
    print('DEBUG: Entering cluster_jedna.cluster()')
    (profiles, merged_mask) = process_profiles(raw_profiles, mask, n_merge)
    d = evaluate_distances(profiles, merged_mask, method='dvr_per_genome')
    clusters = cluster_greedily(d, score_cutoff)
    clusters = cleanup_clusters(clusters)
    print('DEBUG: Leaving cluster_jedna.cluster()')
    return(clusters)
# cluster

################################################################################
## Integrate all pieces
def cluster_per_chrom(raw_profiles, mask, n_merge, score_cutoff=5):
    print('DEBUG: Entering cluster_jedna.cluster_per_chrom()')
    (profiles, merged_mask) = process_profiles(raw_profiles, mask, n_merge)
    d = evaluate_distances(profiles, merged_mask, method='dvr_per_chr')
    clusters_per_chrom = []
    chromosomes = range(len(profiles))
    #
    load_z_flag = False
    for chrom in chromosomes:
        if load_z_flag:
            load_z(chrom)
        else:
            z = get_z(d, chrom)
            save_z(chrom, z)
        # if load_z_else
        clusters = cluster_greedily(z, score_cutoff)
        clusters = cleanup_clusters(clusters)
        clusters_per_chrom.append(clusters)
    # for chrom
    #
    print('DEBUG: Leaving cluster_jedna.cluster_per_chrom()')
    return(clusters_per_chrom)
# cluster_per_chroms

################################################################################
def load_z(chrom):
    print('DEBUG: Entering cluster_jedna.load_z()')
    pass
    print('DEBUG: Leaving cluster_jedna.load_z()')
# load_z
################################################################################
def save_z(chrom, z):
    print('DEBUG: Entering cluster_jedna.save_z()')
    z_file = open('z_chrom%d.py' % (chrom + 1), 'w')
    z_file.write('chrom = %d\n' % chrom)
    z_file.write('I = [%s]\n' % ','.join(map(str, z['I'])))
    z_file.write('J = [%s]\n' % ','.join(map(str, z['J'])))
    z_file.write('DVR = [%s]\n' % ','.join(map(str, z['DVR'])))
    z_file.write('ChromDelta = [%s]\n' % ','.join(map(str, z['ChromDelta'])))
    z_file.write('Z = [%s]\n' % ','.join(map(str, z['Z'])))
    z_file.write('Penalty = [%s]\n' % ','.join(map(str, z['Penalty'])))
    z_file.write('z = {"I": I, "J": J, "DVR": DVR, "ChromDelta": ChromDelta, "Z": Z, "Penalty": Penalty}\n')
    z_file.close()
    print('DEBUG: Leaving cluster_jedna.save_z()')
# load_z

