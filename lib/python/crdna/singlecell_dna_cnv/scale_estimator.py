import numpy as np
from scipy.fftpack import fft, ifft

################################################################################
def get_aggregated_profile(profiles, mask, cluster):
    print('DEBUG Entering scale_estimator.get_aggregated_profile()')
    print('DEBUG scale_estimator.get_aggregated_profile(): cluster size %d' % len(cluster))
    aggregated_profile = []
    n_chrom = len(profiles)
    #
    print('DEBUG scale_estimator.get_aggregated_profile() n_chrom=%d' % n_chrom)
    #
    for chrom in range(n_chrom):
        #
        print('DEBUG chrom=%d' % chrom)
        print('DEBUG cluster:')
        print(cluster)
        print('DEBUG len(cluster)=%d' % len(cluster))
        print('DEBUG len(mask[chrom])=%d' % len(mask[chrom]))
        print('DEBUG len(mask[chrom][mask[chrom]])=%d' % len(mask[chrom][mask[chrom]]))
        #
        tmp = profiles[chrom][np.ix_(cluster, np.where(mask[chrom])[0])].sum(axis=0)
        #
        print('DEBUG tmp.shape:')
        print(tmp.shape)
        print('DEBUG tmp:')
        print(tmp.tolist())
        #
        tmp[np.where(np.isnan(tmp))] = 0.0
        #
        print('DEBUG tmp:')
        print(tmp.tolist())
        #
        aggregated_profile.append(tmp)
        #
        print('chrom=%d' % chrom)
        #
    # for chrom
    #
    print('DEBUG scale_estimator.get_aggregated_profile() Before: len(aggregated_profile)=%d' % len(aggregated_profile))
    #
    aggregated_profile = np.concatenate(aggregated_profile)
    #
    print('DEBUG scale_estimator.get_aggregated_profile() After: len(aggregated_profile)=%d' % len(aggregated_profile))
    print('DEBUG Leaving scale_estimator.get_aggregated_profile()')
    #
    return(aggregated_profile)
# get_aggregated_profile

################################################################################
def get_coarse_profile(aggregated_profile, step=25):
    profile_low_res = []
    n_low_res = int(aggregated_profile.shape[0] / step)

    print('DEBUG scale_estimator.get_coarse_profile() n_low_res=%d' % n_low_res)

    for i in range(n_low_res):
        start = i * step
        end = (i + 1) * step
        profile_low_res.append(np.mean(aggregated_profile[start:end]))
    # for i
    return(profile_low_res)
# get_coarse_profile

################################################################################
def get_count_frequency(aggregated_profile):
    max_count = round(np.max(aggregated_profile)) + 1
    #frequency, counts, patches = plt.hist(aggregated_profile, max_count + 1, range=[0, max_count])
    frequency, counts = np.histogram(aggregated_profile, bins=max_count + 1, range=[0, max_count])
    #counts = np.linspace(start=0, stop=max_count, num=max_count + 1)
    #
    # Take bin midpoints:
    counts = 0.5 * (counts + np.roll(counts, -1))[:-1]
    #
    print('#' * 80)
    print('DEBUG scale_estimator.get_count_frequency() counts:')
    print(counts)
    print('DEBUG scale_estimator.get_count_frequency() frequency:')
    print('............')
    print(frequency)
    #
    return({'Count':counts, 'Frequency':frequency})
# get_count_frequency

################################################################################
def add_spike_at_origin(count_frequency):
    c = count_frequency['Count']
    f = count_frequency['Frequency']
    mass = f.sum()
    f[0] = mass * 2 // 30
    f[1] = mass // 30
    return({'Count':c, 'Frequency':f})
# add_spike_at_origin

################################################################################
#def get_autocorrelation(count_frequency):
#    frequency = count_frequency['Frequency']
#    n_points = frequency.shape[0]
#    #tiled_frequency = np.tile(frequency, 2)
#    auto_cor = np.zeros(n_points)
#    for point in range(n_points):
#        tmp = np.sum(frequency * np.roll(frequency, point))
#        auto_cor[point] = tmp
#    # for point
#    return(auto_cor)
## get_autocorrelation
#
# Better implementation:
def get_autocorrelation(count_frequency):
    fft_freq = fft(count_frequency['Frequency'])
    auto_cor = ifft(fft_freq * np.conj(fft_freq))
    return({'Count':count_frequency['Count'] - np.min(count_frequency['Count']), 
            'AutoCorr':auto_cor.real})
# get_autocorrelation

################################################################################
#def suppress_baseline(auto_corr):
#    c = auto_corr['Count']
#    c2 = -c * (c - max(c))
#    ac = auto_corr['AutoCorr']
#    #plt.plot(c, ac)
#    return({'Count':c, 'AutoCorr':ac * c2})
## supress_baseline

################################################################################
def take_first_half(auto_corr):
    c = auto_corr['Count']
    ac = auto_corr['AutoCorr']
    n_points = c.shape[0] // 2 + 1
    return({'Count':c[:n_points], 'AutoCorr':ac[:n_points]})
# take_first_half

################################################################################
def remove_initial_slope(auto_corr):
    c = auto_corr['Count']
    ac = auto_corr['AutoCorr']
    n_points = c.shape[0]
    da = np.diff(ac)
    for i in range(n_points - 1):
        if da[i] <= 0:
            ac[i] = np.nan
        else:
            break
        # if da else
    # for i
    return({'Count':c, 'AutoCorr':ac})
# remove_initial_slope

################################################################################
def remove_baseline(auto_corr):
    c = auto_corr['Count']
    ac = auto_corr['AutoCorr']
    n_points = c.shape[0]
    ac -= np.nanmean(ac)
    sd = np.nanstd(ac)
    ac[ac <= sd] = np.nan
    ac = ac - np.nanmin(ac)
    ac[np.isnan(ac)] = 0.0
    return({'Count':c, 'AutoCorr':ac})
# remove_baseline

################################################################################
# This assumes that autocorrelation is smooth. 
def find_extrema(auto_corr):
    c = auto_corr['Count']
    ac = auto_corr['AutoCorr']
    da = np.diff(ac)
    n_pts = da.shape[0]
    count = 0
    extrema = []
    for index in range(n_pts - 1):
        if (da[index] * da[index + 1]) < 0:
            zero = 0.5 * (c[index] + c[index + 1]) # One could do a better job interpolating this
            height = ac[index]
            if da[index] > da[index + 1]:
                extremum_type = 'Max'
            else:
                extremum_type = 'Min'
            # if differential else
            extrema.append((zero, extremum_type, height))
        # if product
    # for index
    return(extrema)
# find_extrema

################################################################################
"""
def find_first_extremum(auto_corr):
    c = auto_corr['Count']
    ac = auto_corr['AutoCorr']
    n_pts = ac.shape[0]
    search_start = int(n_pts / 4.0)
    search_end = int(n_pts / 2.0) + 1
    max_index = search_start + np.argmax(ac[search_start:search_end])
    extremum = c[max_index] - c[0]
    return(extremum)
# find_first_extremum
"""

################################################################################
# Running median (see https://bitbucket.org/janto/snippets/src/tip/running_median.py?fileviewer=file-view-default)
from collections import deque
from bisect import insort, bisect_left
from itertools import islice
def running_median_insort(seq, window_size):
    """Contributed by Peter Otten"""
    seq = iter(seq)
    d = deque()
    s = []
    result = []
    for item in islice(seq, window_size):
        d.append(item)
        insort(s, item)
        result.append(s[len(d)//2])
    m = window_size // 2
    for item in seq:
        old = d.popleft()
        d.append(item)
        del s[bisect_left(s, old)]
        insort(s, item)
        result.append(s[m])
    return result
# running_median_insort

################################################################################
def smooth(auto_corr, window_size=25):
    c = auto_corr['Count']
    ac = auto_corr['AutoCorr']
    ac = running_median_insort(ac, window_size)
    return({'Count':c, 'AutoCorr':ac})
# smooth

################################################################################
def get_scale(profiles, mask, cluster, add_spike=True):
    aggregated_profile = get_aggregated_profile(profiles, mask, cluster)
    coarse_profile = get_coarse_profile(aggregated_profile)
    count_frequency = get_count_frequency(coarse_profile)
    if add_spike:
        count_frequency = add_spike_at_origin(count_frequency)
        print('DEBUG scale_estimator.get_scale() Added spike!')
    # if add_spike
    auto_corr = get_autocorrelation(count_frequency)
    auto_corr = smooth(auto_corr)
    auto_corr = take_first_half(auto_corr)
    auto_corr = remove_initial_slope(auto_corr)
    auto_corr = remove_baseline(auto_corr)
    # scale = find_first_extremum(auto_corr)

    print('x' * 80)
    print('DEBUG scale_estimator.get_scale() Counts:')
    print(auto_corr['Count'].tolist())
    print('DEBUG scale_estimator.get_scale() AutoCorr:')
    print(auto_corr['AutoCorr'].tolist())

    extrema = find_extrema(auto_corr)

    print('DEBUG scale_estimator.get_scale() extrema:')
    print(extrema)

    n_extrema = len(extrema)
    if n_extrema < 1:
        scale = np.nanmean(coarse_profile)
    else:
        max_index = np.argmax([extrema[i][2] for i in range(n_extrema)])
        scale = extrema[max_index][0] # This corresponds to a single chromosome copy; need to multiply by 2 to handle diploid regions
    # if n_extrema else
    #
    print('DEBUG scale_estimator.get_scale() scale=%f' % scale)
    #
    return(2.0 * scale) # Multiplication by 2 scales diploid regions to 1
    #return(7.75 * 25)
# get_scale




