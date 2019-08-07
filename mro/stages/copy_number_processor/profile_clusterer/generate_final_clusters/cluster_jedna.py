# cluster_jedna.py

import math
import numpy
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
    for chrom in range(n_chrom):
        chr_read_counts = numpy.sum(numpy.transpose(profiles[chrom])[mask[chrom]], axis=0)
        total_read_counts = [total_read_counts[i] + chr_read_counts[i] 
            for i in range(n_cells)]
    # for chrom
    return(total_read_counts)  
# get_total_read_couns

################################################################################
##
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
    # Scale with respect to total reads
    profiles = scale(profiles, merged_mask)
    return(profiles, merged_mask)
# process_profiles

################################################################################
## Evaluate distances
def evaluate_distances(profiles, mask):
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
    n_cells = len(profiles[0])
    if n_cells < 1:
        # TODO: warning: no cells
        return(d)
    # if n_celss
    #...........................................................................
    def distance(cell_i, cell_j):
        delta2 = [None] * n_chrom
        for chrom in range(n_chrom):
            delta = numpy.subtract(profiles[chrom][cell_i, ], profiles[chrom][cell_j, ])
            delta2[chrom] = math.sqrt(numpy.mean(delta[mask[chrom]] * delta[mask[chrom]]))
        # for chr
        return(delta2)
    # distance
    #
    n_pairs = int(0.5 * n_cells * ( n_cells - 1 ))
    d = {"I":([None] * n_pairs),
         "J":([None] * n_pairs),
         "D":([None] * n_pairs)}
    count = 0
    for cell_i in range(0, n_cells - 1):
        for cell_j in range(cell_i + 1, n_cells):
            d["I"][count] = cell_i
            d["J"][count] = cell_j
            d["D"][count] = distance(cell_i, cell_j)
            count += 1
        # for j
    # for i
    return(d)
# evaluate_distances

################################################################################
## Standardize distances
def standardize_distances(d):
    #...........................................................................
    def get_distances_per_chrom(chrom, d):
        d_per_chrom = [x[chrom] for x in d.get("D")]
        return(d_per_chrom)
    # get_distances_per_chrom
    #
    #...........................................................................
    def get_z_max(z_per_chrom):
        chromosomes = range(len(z_per_chrom))
        n = len(z_per_chrom.get(chromosomes[0]))
        z_max = [None] * n
        for index in range(0, n):
            z_max[index] = max([z_per_chrom.get(chrom)[index] for chrom in chromosomes])
        # for index
        return(z_max)
    # get_z_max
    #
    chromosomes = range(len(d.get("D")[0]))
    z_per_chrom = {}
    for chrom in chromosomes:
        d_per_chrom = get_distances_per_chrom(chrom, d)
        mean_d_per_chrom = numpy.mean(d_per_chrom)
        sd_d_per_chrom = numpy.std(d_per_chrom)
        z_per_chrom[chrom] = (d_per_chrom - mean_d_per_chrom) / sd_d_per_chrom
    # for chr
    z_max = get_z_max(z_per_chrom)
    z = {"I":d.get("I"), "J":d.get("J"), "Z":z_max}
    return(z)
# standardize_distances

################################################################################
## Greedy association
def cluster_greedily(z, score_cutoff=10):
    #epsilon = 1e6
    profile_ids = range(z.get("J")[-1] + 1)
    clusters = [] 
    seed = numpy.argmin(z.get("Z")) 
    cluster = [z.get("I")[seed], z.get("J")[seed]] 
    profile_ids.remove(z.get("I")[seed])
    profile_ids.remove(z.get("J")[seed])
    print("=" * 80)
    print("New cluster: [TODO].\nClustered: %d. Unclustered: %d. Score: [%f]" % 
        (len(cluster), len(profile_ids), z.get("Z")[seed])) # TODO: display cluster
    min_d = numpy.percentile(z.get("Z"), score_cutoff)
    while len(profile_ids) > 0:
         #
         # Proximity within cluster
         cluster_size = len(cluster)
         n_pairs = cluster_size * (cluster_size - 1) / 2
         tmp = [None] * n_pairs
         count = 0
         for k in range(0, cluster_size - 1):
             for s in range(k + 1, cluster_size):
                 auxiliary = map(lambda x, y:((x == cluster[k] and y == cluster[s]) or
                                              (x == cluster[s] and y == cluster[k])), 
                                              z.get("I"), z.get("J"))
                 selector = which(auxiliary)
                 tmp[count] = z.get("Z")[selector[0]]
                 count += 1
             # for s
         # for k
         #d_within = numpy.mean(tmp)
         #s_within = numpy.std(tmp)
         #
         # Proximity between cluster members and outside indices
         d_between = [None] * len(profile_ids)
         #s_between = [None] * len(profile_ids)
         index = 0
         for out_id in profile_ids:
             tmp = [None] * cluster_size
             count = 0
             for cluster_id in cluster:
                 auxiliary = map(lambda x, y:((x == out_id and y == cluster_id) or
                                              (x == cluster_id and y == out_id)),
                                              z.get("I"), z.get("J"))
                 selector = which(auxiliary)
                 tmp[count] = z.get("Z")[selector[0]]
                 count += 1
             # for cluster_id
             d_between[index] = numpy.mean(tmp)
             #s_between[index] = numpy.std(tmp)
             index += 1
         # for out_id
         #
         # Find the best candidates
         auxiliary = map(lambda x:x <= min_d, d_between)
         k_best = which(auxiliary)
         #
         # Decision:
         if len(k_best) > 0:
             #
             # If yes, augment the cluster
             new_ids = [profile_ids[i] for i in k_best]
             cluster.extend(new_ids)
             [profile_ids.remove(out_id) for out_id in new_ids]
             print("Growing cluster: [TODO].\nClustered: %d. Unclustered: %d. Score: [%f, %f]" % # TODO: display cluster
                  (len(cluster), len(profile_ids), min(d_between), min_d))
         # if k_best
         #
         # Close this cluster
         clusters.append(cluster)
         #
         # Seed another cluster
         auxiliary = map(lambda x, y:(x in profile_ids) and (y in profile_ids),
                                     z.get("I"), z.get("J"))
         selector = which(auxiliary)
         seed = numpy.argmin([z.get("Z")[i] for i in selector])
         seed_index = selector[seed]
         min_z = z.get("Z")[seed_index]
         print("=" * 80)
         if min_z > min_d:
             # Further grouping is too loose - treat all remaining barcodes as single-cell clusters
             for singlet_id in profile_ids:
                 clusters.append([singlet_id])
                 print("Single cell cluster: [TODO].") # TODO: display cluster
             # for singlet_id
             profile_ids = []
         else:
             cluster = [z.get("I")[seed_index], z.get("J")[seed_index]]
             profile_ids.remove(z.get("I")[seed_index])
             profile_ids.remove(z.get("J")[seed_index])
             print("New cluster: [TODO]. Clustered: %s. Unclustered: %d. Score: [%f]" % 
                 ("TODO", len(profile_ids), min_z)) # TODO: display cluster, #clustered cells
         # if z else
    # while profile_ids
    #
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
#             outliers[[ outCounter ]] <- list( "Chr"=chr, "Outliers"=outsiders );
#         } # if outsiders
#     } # for chr
#     #
#     # Remove outliers from the cluster
#     # TODO
#     #
#     result <- list( "Cluster"=cluster, "Outliers"=outliers );
#     return( result );
# } # findOutliers
# 
# Cluster cleanup
def cleanup_clusters(clusters):
    # TODO: Implement
    return(clusters)
# cleanup_clusters

################################################################################
## Integrate all pieces
def cluster(raw_profiles, mask, n_merge, score_cutoff=10):
    (profiles, merged_mask) = process_profiles(raw_profiles, mask, n_merge)
    d = evaluate_distances(profiles, merged_mask)
    z = standardize_distances(d)
    clusters = cluster_greedily(z, score_cutoff)
    clusters = cleanup_clusters(clusters)
    return(clusters)
# cluster
