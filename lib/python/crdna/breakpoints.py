## Support functions for CNV calling through breakpoints
import numpy as np
import scipy.stats
import bisect
import sys
import scipy.optimize
from collections import defaultdict, Counter
from constants import MIN_CURVE, MAX_CURVE
from plot_utils import aggregate_counts
import time
import pandas as pd

def parabola(x, o, l, q):
    """ Parabola representing GC bias model. """
    return np.clip(1 + (x-o)*l + np.power(x-o, 2)*q, MIN_CURVE, MAX_CURVE)

def max_look(y):
    return int(np.ceil(1.0*len(y)/10))

def llscore6(y, xi, read_threshold=200, MIN_LOOK = 5,
    log_func=sys.stdout.write):
    """
    Computes a log likelihood ratio testing whether the bins to the left of
    a position x and to the right of x are drawn from the same poisson
    distribution. Note, there is no filtering of outliers. And this includes
    GC correction.
    """
    ## how far to look in each direction
    MAX_LOOK = max_look(y)
    mean_reads = y.mean( )
    look = int(np.round(read_threshold/mean_reads)) if mean_reads > 0 else (MAX_LOOK + 1)

    if look > MAX_LOOK:
        sys.stderr.write("look = %d, exceeds MAX_LOOK %d\n"%(look, MAX_LOOK))
    look = np.clip(look, MIN_LOOK, MAX_LOOK)

    log_func("look = %d\n"%look)

    nbins = len(y)
    ll = np.zeros(nbins, dtype=float)

    y2 = np.zeros(nbins+2*look, dtype=y.dtype)
    xi2 = np.zeros(nbins+2*look, dtype=xi.dtype)

    y2[look:nbins+look] = y
    y2[:look] = y[-look:]
    y2[nbins+look:] = y[:look]

    xi2[look:nbins+look] = xi
    xi2[:look] = xi[-look:]
    xi2[nbins+look:] = xi[:look]

    for pos in xrange(look, nbins+look):
        lindex = pos-look
        rindex = pos+look

        lxi = xi2[lindex:pos]
        rxi = xi2[pos:rindex]

        lysum = y2[lindex:pos].sum()
        rysum = y2[pos:rindex].sum()
        lxisum = lxi.sum()
        rxisum = rxi.sum()

        lmean = max(lysum/lxisum, 1e-6)
        rmean = max(rysum/rxisum, 1e-6)
        amean = max((lysum+rysum)/(lxisum+rxisum), 1e-6)

        loglmean = np.log(lmean)
        logrmean = np.log(rmean)
        logamean = np.log(amean)

        ll[pos-look] = ((amean-lmean)*lxisum + (amean-rmean)*rxisum +
                        (loglmean-logamean)*lysum +
                        (logrmean-logamean)*rysum)

    return ll

def add_missed_breakpoints2(y, xi, segment_bdy2, min_segment_size = 5,
        min_ll = 5.0, log_func = sys.stdout.write):
    """ Improve segmentation by adding missed breakpoints. For any given 
        segment, see if dividing into segments results in a better log
        likelihood score. If so, split them up. Continue iteratively until
        no more improvements can be made. FASTER version of add_missed_breakpoints.
    """
    itercount = 0
    seen = set([])
    while True:
        itercount += 1
        
        log_func("Iteration: %d, %d segments\n"%(itercount, len(segment_bdy2)))
        
        extras = []
        for lpos, rpos in segment_bdy2:
            if rpos - lpos < min_segment_size:
                continue
            if (lpos, rpos) in seen:
                continue
            seen.add((lpos, rpos))
           
            xi_seg = xi[lpos:rpos]
            y_seg = y[lpos:rpos]

            lysum  = np.cumsum(y_seg).astype(float)
            rysum  = lysum[-1] - lysum
            lxisum = np.cumsum(xi_seg)
            rxisum = lxisum[-1] - lxisum

            lmean = np.clip(lysum.astype(float)/lxisum, 1e-2, None)
            rmean = np.clip(rysum.astype(float)/rxisum, 1e-2, None)
            amean = np.clip(lysum[-1]/lxisum[-1], 1e-2, None) ## just a number

            loglmean = np.log(lmean)
            logrmean = np.log(rmean)
            logamean = np.log(amean)

            delta_lls = ((amean-lmean)*lxisum + (amean-rmean)*rxisum +
                            (loglmean-logamean)*lysum +
                            (logrmean-logamean)*rysum)[:-2]
            
            mi = np.argmax(delta_lls)
            delta_ll_max = delta_lls[mi]
            if delta_ll_max > min_ll:
                extras.append(lpos+1+mi)
        new_segment_bdy = break_segments_at_points(segment_bdy2, extras,
            verbose=False)
        if len(extras) == 0:
            break
        segment_bdy2 = list(map(tuple, new_segment_bdy))

    return new_segment_bdy

def delete_unneeded_breakpoints2(y, xi, segment_bdy1, ll, 
        delta_ll_threshold=2.0, log_func=sys.stdout.write):
    """
        Analyze the left and right segments adjacent to a breakpoint and
        decide if the breakpoint is needed there. Do this based on the
        log likelihood score.
    """
    bps = np.array([x[0] for x in segment_bdy1] + [segment_bdy1[-1][1]])
    keep = np.ones_like(bps, dtype=bool)
    keep_going = True
    order = list(np.argsort(ll[bps[:-1]]))
    while keep_going:
        keep_going = False
        for i in order:
            pos = bps[i]
            if pos == 0 or pos == len(bps)-1:
                continue
            if not keep[i]:
                continue
            left = i-1
            while not keep[left]:
                left -= 1
            assert left >=0
            right = i+1
            while not keep[right]:
                right += 1
            assert right < len(bps)

            lpos = bps[left]
            rpos = bps[right]
            lpts = y[lpos:pos]
            rpts = y[pos:rpos]
            allpts = y[lpos:rpos]

            lxi = xi[lpos:pos]
            rxi = xi[pos:rpos]
            allxi = xi[lpos:rpos]

            mu_l = np.clip(lpts.sum( )/lxi.sum( ), 1e-2, None)*lxi
            mu_r = np.clip(rpts.sum( )/rxi.sum( ), 1e-2, None)*rxi
            mu_all = np.clip(allpts.sum( )/allxi.sum( ), 1e-2, None)*allxi

            delta_ll = ((-mu_l + lpts*np.log(mu_l)).sum( ) +
                   (-mu_r + rpts*np.log(mu_r)).sum( ) -
                   (-mu_all + allpts*np.log(mu_all)).sum( ))
            if delta_ll < delta_ll_threshold:
                keep[i] = False
                keep_going = True
    bps2 = bps[keep]
    log_func("Retaining %.2f pct of breakpoints\n"%(keep.sum()*100.0/keep.shape[0]))
    segment_bdy2 = [(bps2[i-1], bps2[i]) for i in xrange(1, len(bps2))]
    return segment_bdy2

def wiggle_breakpoints(y, xi, segment_bdy3, wiggle_width=5, num_iterations=1,
        verbose=False, log_func=sys.stdout.write):
    """ Wiggle each breakpoint by +/- wiggle_width and see if it produces
        a better solution based on the log likelihood ratio test."""
    t0 = time.time()
    count = 0
    segment_bdy4 = map(tuple, segment_bdy3)
    while count < num_iterations:
        did_nothing = True
        bps = [x[0] for x in segment_bdy4] + [segment_bdy4[-1][1]]
        new_bps = [0]
        for index in xrange(1, len(bps)-1):
            b = bps[index]
            lpos = new_bps[index-1]
            rpos = bps[index+1]

            wiggle_width = 5
            if b - lpos == 1 or rpos - b == 1:
                new_bps.append(b)
                continue
            lwiggle = max(b-wiggle_width, lpos+1)
            rwiggle = min(b+wiggle_width, rpos-1)
            delta_lls = []
            for pos in xrange(lwiggle, rwiggle):
                lpts = y[lpos:pos]
                rpts = y[pos:rpos]
                allpts = y[lpos:rpos]

                lxi = xi[lpos:pos]
                rxi = xi[pos:rpos]
                allxi = xi[lpos:rpos]

                mu_l = np.clip(lpts.sum( )/lxi.sum( ), 1e-2, None)*lxi
                mu_r = np.clip(rpts.sum( )/rxi.sum( ), 1e-2, None)*rxi
                mu_all = np.clip(allpts.sum( )/allxi.sum( ), 1e-2, None)*allxi

                delta_ll = ((-mu_l + lpts*np.log(mu_l)).sum( ) +
                       (-mu_r + rpts*np.log(mu_r)).sum( ) -
                       (-mu_all + allpts*np.log(mu_all)).sum( ))
                delta_lls.append(delta_ll)
            new_b = lwiggle + np.argmax(delta_lls)
            try:
                gain = (delta_lls[new_b-lwiggle] - delta_lls[b-lwiggle])
            except:
                print "-"*40
                print lpos, b, rpos
                print "argmax", np.argmax(delta_lls), len(delta_lls)
                print "wiggle", lwiggle, rwiggle
                print new_b-lwiggle, b-lwiggle
                raise Exception("Wiggling breakpoints produced an invalid breakpoint configuration")
            new_bps.append(new_b)
            if new_b != b:
                did_nothing = False
                if verbose:
                    print "%6d -> %6d : (%2d) gain +%.2f"%(b, new_b, delta_lls[b-lwiggle], gain)
        new_bps.append(bps[-1])
        assert len(new_bps) == len(bps)
        segment_bdy4 = [(new_bps[i], new_bps[i+1]) for i in xrange(0, len(new_bps)-1)]
        count += 1
        if did_nothing:
            break
    log_func("%.2f s spent wiggling, %d iterations\n"%(time.time()-t0, count))
    return segment_bdy4

def get_breakpoint_positions( y, ll, xi, ll_threshold = 5):
    """ Determine breakpoint positions given the cell profile and the 
        log likelihood score. The GC correction is contained in xi.
    """

    EPS = 1e-12
    ## first get the local maxima of the log likelihoods
    bp_cands = [0]
    for i in xrange(1, ll.shape[0]-1):
        delta_ll_left = ll[i] - ll[i-1]
        delta_ll_right = ll[i] - ll[i+1]
        if (delta_ll_left > EPS and delta_ll_right > EPS and 
            ll[i] > ll_threshold):
            bp_cands.append(i)
    bp_cands.append(ll.shape[0])
    bp_cands = np.array(bp_cands)
    return bp_cands

def break_segments_at_points( intervals, positions, verbose=False ):
    """ Given a list of [start, end) intervals that completely cover a half-open interval,
        introduce breaks at positions specified by positions. 
    """
    starts, ends = map(list, zip(*intervals))
    for chrom_end_pos in positions:
        assert chrom_end_pos <= ends[-1], "segment end does not include end of genome"
        assert chrom_end_pos >= starts[0], "start of genome to the left of first segment"
        spos = bisect.bisect_left(starts, chrom_end_pos)

        if verbose:
            print "-"*80
            print "position %d"%chrom_end_pos
            print "spos %d"%spos
        if spos < len(starts) and starts[spos] == chrom_end_pos:
            if verbose:
                print "already a segment boundary"
        else:
            epos = bisect.bisect_left(ends, chrom_end_pos)
            if verbose:
                print "epos %d"%epos
            assert spos == epos+1
            if chrom_end_pos < ends[epos]:
                if verbose:
                        print "Break [%d, %d) into [%d, %d) and [%d, %d)"%(starts[spos-1], ends[epos], 
                                                                           starts[spos-1], chrom_end_pos,
                                                                           chrom_end_pos, ends[epos])
                        print "Delete %d,%d"%(starts[epos], ends[epos])

                # update epos
                copy_end = ends[epos]
                ends[epos] = chrom_end_pos
                # update epos+1
                starts.insert(epos+1, chrom_end_pos)
                ends.insert(epos+1, copy_end)
                if verbose:
                    print "epos: ", (starts[epos], ends[epos])
                    print "epos+1: ", (starts[epos+1], ends[epos+1])
            elif chrom_end_pos == ends[epos]:
                if verbose:
                    print "already a segment boundary"
    return zip(starts, ends)

## Compute the ploidy vector
def get_ploidy_vector( y, segment_means, segment_bdy, lam ):
    segment_ploidy = np.round(segment_means/lam).astype(int)
    ploidy = -np.ones_like(y)
    for (s,e), p in zip(segment_bdy, segment_ploidy):
        ploidy[s:e] = p
    assert (ploidy < 0).sum( ) == 0, "some bin ploidies not assigned"
    return ploidy

def validate_segment_intervals(segment_bdy, cbdy):
    assert segment_bdy[0][0] == cbdy[0], "segments must start at genome start"
    assert segment_bdy[-1][1] == cbdy[-1], "segments must end where genome ends"
    for i in xrange(len(segment_bdy)-1):
        assert segment_bdy[i][1] == segment_bdy[i+1][0], "segments must be abutting:"\
        " (%d,%d) and (%d,%d)"%(segment_bdy[i][0],
        segment_bdy[i][1],
        segment_bdy[i+1][0],
        segment_bdy[i+1][1])
    total = 0
    for s, e in segment_bdy:
        length = e-s
        assert length > 0, "segments must have nonzero length"
        total += length
    assert total == cbdy[-1], "segments must cover the whole genome"

def objective_function(lam, segment_means, segment_lengths):
    """ Objective function to minimize to determine integer ploidies of segments
    """
    return (segment_lengths*np.power(np.sin(segment_means*np.pi/lam), 2)).sum( )

def log_prior_ploidy(p, mu = 2.0, sig = 0.5):
    """ Prior distribution for ploidy. Assumed to be gaussian."""
    log_val  = scipy.stats.norm.logpdf(p, mu, sig)
    log_norm = np.log1p(np.exp(scipy.stats.norm.logcdf(0.0, mu, sig)))
    return log_val - log_norm

def find_scaling_candidates_v4(y, segment_bdy, segment_means, segment_lengths,
    is_sex_chrom_bin, window, scale_guess, max_ploidy_long=10,
    zero_ploidy_count=50, max_segment_push_to_zero=50, verbose=True,
    log_func=sys.stdout.write):
    """ Find a set of possible scaling factors by minimizing an
        objective function. In addition to the scaling factors, we also
        store other information corresponding to the scale. Same function
        as find_scaling_candidates, but output is a list of dicts."""

    data_cats = ["fopt", "lam", "aploidy", "dom_pdy",
                               "dom_frac", "levels", "entropy"]

    ## create initial guesses by setting the average ploidy based on
    ## the read count
    mean_read_count = y.mean()*window
    possible_ploidies = np.linspace(.5, 12, 20)
    guesses = mean_read_count/possible_ploidies

    ## focus on good segments for scaling (based on dpcv^2)
    dpcv2 = -np.ones_like(segment_means)
    for i, (s,e) in enumerate(segment_bdy):
        counts = y[s:e]
        mean = counts.mean()
        dpcv2[i] = (counts.std()**2 - mean)/np.clip(mean, 1e-2, None)
    median_length = np.median(segment_lengths)
    dpcv2_25_long, dpcv2_75_long = np.percentile(
        dpcv2[segment_lengths >= median_length], [25, 75])
    good_segments = dpcv2 < dpcv2_75_long + 1.5*(dpcv2_75_long-dpcv2_25_long)
    if verbose:
        log_func("Using %.2f pct of segments to scale\n"%(
            good_segments.sum()*100./len(good_segments)))
    segment_means_filtered = segment_means[good_segments]
    segment_lengths_filtered = segment_lengths[good_segments]

    ## Run a minimization algorithm to determine the candidate scaling
    ## factors.
    data_dict = defaultdict(list)
    for g in guesses:
        lams, fopt, _, _, warnflag = scipy.optimize.fmin(lambda lam :
              objective_function(lam, segment_means_filtered,
                                 segment_lengths_filtered),
              g, full_output=True, disp=False)
        lam = lams[0]
        log_func("Guess %.2f, lam %.2f"%(g, lam))

        ## if the minimization failed then stop
        if (warnflag == 1 or warnflag == 2 or (len(guesses) > 1 and 
            (lam > guesses.max() or lam < guesses.min()))):
            log_func(" SKIP\n")
            continue
        log_func(" KEEP\n")
        segment_ploidies = np.round(segment_means/lam).astype(int)

        #pushed_to_zero = np.where((segment_ploidies == 0) &
        #                          (segment_means > zero_ploidy_count))[0]
        #if (len(pushed_to_zero) and
        #    segment_lengths[pushed_to_zero].max( ) > max_segment_push_to_zero):
        #    continue
        
        pcounter = Counter()
        for p, l in zip(segment_ploidies, segment_lengths):
            pcounter[p] += l
        nploidy_levels = len(pcounter.keys( ))

        ## compute the information content of the CNV profile
        entropy_counter = Counter( )
        entropy_norm = 0
        for (s,e), p in zip(segment_bdy, segment_ploidies):
            l = e-s
            ## skip sex chromosomes so both male and female diploids
            ## will have low entropy
            if is_sex_chrom_bin[s:e].sum() > 0.5*l:
                continue
            entropy_counter[p] += l
            entropy_norm += l
        if entropy_norm == 0:
            # if no segment is long, then set to arbitrary high value
            entropy = np.log(len(y)) 
        else:
            entropy = 0.0
            for p, c in entropy_counter.iteritems( ):
                occupancy = float(c)/entropy_norm
                entropy += -occupancy * np.log(occupancy)

        ## find the dominant ploidy and how much of the genome it covers
        for dom_ploidy, dom_count in pcounter.most_common(2):
            if int(dom_ploidy) != 0:
                break

        avg_ploidy = (segment_ploidies.astype(float)*
                      segment_lengths).sum( )/segment_lengths.sum( )
        key = int(np.round(lam))
        dom_frac = dom_count*1.0/segment_lengths.sum()
        data_dict[key].append((fopt, lam, avg_ploidy, dom_ploidy,
                               dom_frac, nploidy_levels, entropy))

    ## if we were not able to find a single scale factor, then
    ## assign the initial segment ploidy 2
    if len(data_dict.keys()) == 0:
        if verbose:
            log_func("Optimization failed, setting median to ploidy 2\n")
        y_agg = aggregate_counts(y, window)
        y_agg_med = np.median(y_agg[y_agg > 0])
        lam = y_agg_med/2.0
        return [dict(zip(data_cats, [None, lam, None,
            None, None, None, None]))]

    data = []
    for k, v in data_dict.iteritems( ):
        data.append(dict(zip(data_cats, np.array(v).mean(axis=0))))

    data.sort( key = lambda x: x["fopt"] )

    if scale_guess is not None and len(scale_guess) == 1:
        sf = scale_guess[0]*window
        best_match = np.argmin([abs(d["lam"] - sf) for d in data])
        if verbose:
            log_func("Override scaling solution. Pick closest to %.2f\n"%sf)
        return [data[best_match]]

    if verbose:
        log_func("%d minima found\n"%len(data))
    return data

def find_best_scale_v15(y, xi, segment_bdy, segment_means, segment_lengths,
    ref, cbdy, window, scale_guess, max_ploidy_long=10, zero_ploidy_count=200, 
    max_segment_push_to_zero=200, prior_params = {"prior_mean": 2.0,
    "prior_std": 0.5}, min_ploidy = None, max_ploidy = None, verbose=True,
    log_func=sys.stdout.write):
    """ Find the best scaling factor that turns a segmented profile
        into an integer copy number profile. min_ploidy and max_ploidy override
        the algorithm. See comments for explanation."""

    def print_scaling_solutions(scaling_data, order, best_index=None):
        log_func("Scaling solutions:\n")
        header = ["lam", "aploidy", "dom_pdy",
                  "dom_frac", "lprior", "entropy", "nfit", "dpcv"]
        log_func("|".join([h.rjust(8) for h in header])+"\n")
        log_func("-"*(len("|".join([h.rjust(8) for h in header])))+"\n")
        for x in order:
            if scaling_data[x]["fopt"] is None:
                log_func("NO DATA\n")
                continue
            row = []
            for h in header:
                value = ("%.2f"%scaling_data[x][h] if h in scaling_data[x]
                    else "NAN")
                row.append(value.rjust(8))
            log_func(" ".join(row)+" ")
            if best_index is not None and best_index == x:
                log_func("X")
            log_func("\n")
    
    is_sex_chrom_bin = np.zeros_like(y, dtype=bool)
    
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)
    for i, chrom in enumerate(chroms):
        if chrom in ref.sex_chromosomes:
            is_sex_chrom_bin[cbdy[i]:cbdy[i+1]] = True


    ## assertions
    assert prior_params["prior_mean"] > 0, "prior mean must be positive"
    assert prior_params["prior_std"] > 0, "prior std. dev. must be positive"
    
    ## find all the candidate scaling factors
    scaling_data0 = find_scaling_candidates_v4(y, segment_bdy, segment_means,
        segment_lengths, is_sex_chrom_bin, window, scale_guess,
        max_ploidy_long=max_ploidy_long, zero_ploidy_count=zero_ploidy_count,
        max_segment_push_to_zero=max_segment_push_to_zero, verbose=verbose,
        log_func=log_func)
    log_func("INITIAL CANDIDATES\n")
    print_scaling_solutions(scaling_data0, range(len(scaling_data0)), None)

    ## if we only find one solution, there's nothing to do
    if len(scaling_data0) == 1:
        log_func("Only 1 solution, nothing to do.\n")
        if verbose:
            print_scaling_solutions(scaling_data0, [0], 0)
        return scaling_data0[0]["lam"], -3, pd.DataFrame(scaling_data0)
    nbins = len(y)
    ## compute DPCV after compensating for ploidy & GC
    for si, data in enumerate(scaling_data0):
        lam = data["lam"]
        
        segment_ploidy = np.round(segment_means/lam).astype(int)

        blocks = defaultdict(list)
        for i, p in enumerate(segment_ploidy):
            blocks[p].append(segment_bdy[i])
        blen = []
        bdpcv = []
        for p, segs in blocks.iteritems():
            seg_len = 0
            for s, e in segs:
                seg_len += e-s
            if seg_len < 0.01*nbins:
                continue
            counts = np.zeros(seg_len)
            pos = 0
            for s, e in segs:
                counts[pos:pos+(e-s)] = y[s:e]/xi[s:e]/max(p, 0.1)
                pos += e-s
            counts_agg = aggregate_counts(counts, window)
            dpcv = np.sqrt(np.abs(counts_agg.std()**2 -
                counts_agg.mean()))/counts_agg.mean()
            blen.append(seg_len)
            bdpcv.append(dpcv)
        blen = np.array(blen)
        bdpcv = np.array(bdpcv)
        data["dpcv"] = (blen*bdpcv).sum()/blen.sum()
        data["lprior"] = log_prior_ploidy(data["aploidy"],
            mu=prior_params["prior_mean"], sig=prior_params["prior_std"])
        data["nfit"] = data["fopt"]/segment_lengths.sum( )
    
    ## Throw away candidates with unusually high DPCV. This removes the ploidy
    ## 1 solutions that we see when there is noise. Also remove high ploidy 
    ## solutions.
    MIN_PLOIDY = 0.0
    MAX_PLOIDY = 8.0
    ploidies = np.array([data["aploidy"] for data in scaling_data0])
    pfilter = ploidies < MAX_PLOIDY

    ## Filter based on the dpcv of likely good solutions. This is mainly to 
    ## eliminate the ploidy 1 case. Take the min of top 3 solutions, since 
    ## very high ploidy solutions will have low dpcv as well, but they are
    ## garbage.
    dpcvs = np.array([data["dpcv"] for data in scaling_data0])
    target_dpcv = np.min(dpcvs[:3])
    entropy = np.array([data["entropy"] for data in scaling_data0])
    target_entropy = np.max(entropy[:3])
    outliers = np.logical_or(dpcvs > target_dpcv + 0.06,
                            entropy < target_entropy - 1.0)
     
    select = ~outliers & pfilter
    if select.sum() == 0:
        if pfilter.sum() > 0:
            select = pfilter
        else:
            select = np.ones(len(scaling_data0), dtype=bool)

    scaling_data = []
    for i, data in enumerate(scaling_data0):
        if select[i]:
            scaling_data.append(data)
           
    ## if min_ploidy and max_ploidy are specified, override the decision making
    ## - find all solutions with ploidy (min_ploidy, max_ploidy)
    ## - if empty find the solutions that are closest to the end points
    ## - amongst the restricted solutions as determined above, pick the solution
    ##   with the least objective function value
    if min_ploidy is not None or max_ploidy is not None:
        pmin = min_ploidy if min_ploidy is not None else MIN_PLOIDY
        pmax = max_ploidy if max_ploidy is not None else MAX_PLOIDY
        assert pmin <= pmax, "min_ploidy <= max_ploidy must be true"

        log_func("Applying min ploidy & max ploidy: %d-%d\n"%(pmin, pmax))

        filtered_data = [data for data in scaling_data
                         if (data["aploidy"] > pmin and
                             data["aploidy"] < pmax)]
        if len(filtered_data) == 0:
            ploidies = np.array([data["aploidy"] for data in scaling_data0])
            pmin_nbr = np.argmin(np.abs(ploidies - pmin))
            pmax_nbr = np.argmin(np.abs(ploidies - pmax))
            if (scaling_data0[pmin_nbr]["fopt"] <
                    scaling_data0[pmax_nbr]["fopt"]):
                best_index = pmin_nbr
            else:
                best_index = pmax_nbr

            if verbose:
                print_scaling_solutions(scaling_data0, [best_index], best_index)
            return (scaling_data0[best_index]["lam"], -1, 
                pd.DataFrame(scaling_data0))
        else:
            best_index = np.argmin([data["fopt"] for data in filtered_data])
            if verbose:
                print_scaling_solutions(filtered_data,
                    range(len(filtered_data)), best_index)
            return (filtered_data[best_index]["lam"], -1, 
                pd.DataFrame(filtered_data))

    ## determine if the solutions are "degenerate", i.e., almost euploid.
    ## if degenerate, then just pick the solution that optimizes the prior
    ## otherwise, pick the best solution or the second best solution based
    ## on whether the second best solution (in terms of the objective
    ## function value) actually is a better fit in terms of fewer segments
    ## whose ploidies don't nicely match the data

    order = np.argsort([x["nfit"] for x in scaling_data])

    topn = order[:3]
    median_dom_frac = np.median([scaling_data[x]["dom_frac"] for x in topn])
    is_degenerate = (median_dom_frac > 0.90)
    if is_degenerate:
        best_prior = np.argmax([x["lprior"] for x in scaling_data0])
        gap = -2
        ## pick the solution with the best prior
        if verbose:
            log_func("DEGENERATE\n")
            print_scaling_solutions([scaling_data0[best_prior]], [0], 
                best_prior)
        return scaling_data0[best_prior]["lam"], gap, pd.DataFrame(scaling_data0)
    
    ## Pick the best fit solution
    best_index = order[0]
    if len(order) > 1:
        gap = scaling_data[order[1]]["nfit"] - scaling_data[best_index]["nfit"]
    else:
        gap = -3

    if verbose:
        print_scaling_solutions(scaling_data, order, best_index)
    
    return scaling_data[best_index]["lam"], gap, pd.DataFrame(scaling_data0)

def call_cnvs(y, xi, ref, cbdy, scale_guess=None, log_func=sys.stdout.write,
    **kwargs):
    ## Compute log likelihood score
    t0 = time.time( )
    ll = llscore6(y, xi, read_threshold=kwargs["ll_read_threshold"], 
        MIN_LOOK = kwargs["ll_min_look"], log_func=log_func)

    ## Heuristics to define breakpoints
    ll_threshold    = kwargs["ll_threshold"]

    ## Define breakpoints
    bp_cands2 = get_breakpoint_positions(y, ll, xi,
        ll_threshold=ll_threshold)

    assert bp_cands2[0] == 0, "genome start must be breakpoint"
    assert bp_cands2[-1] == y.shape[0], "genome end must be breakpoint"

    ## define segments using breakpoints
    segment_bdy0 = []
    for j in xrange(len(bp_cands2)-1):
        segment_bdy0.append((bp_cands2[j], bp_cands2[j+1]))
    log_func("Initial segment count %d\n" % len(segment_bdy0))
    
    ## add chromosome boundaries as mandatory breakpoints
    segment_bdy1 = break_segments_at_points(segment_bdy0, cbdy, verbose=False)

    ## clean up some breakpoints
    segment_bdy2 = delete_unneeded_breakpoints2(y, xi, segment_bdy1, ll, 
        delta_ll_threshold = kwargs["del_bp_ll_threshold"], log_func=log_func)
    log_func("Segments after deletion %d\n" % len(segment_bdy2))

    ## add missed breakpoints
    segment_bdy3 = add_missed_breakpoints2(y, xi, segment_bdy2,
        min_segment_size=kwargs["add_bp_min_segment_size"], 
        min_ll = kwargs["add_bp_min_ll"], log_func=log_func)
    log_func("Adding back segments, total %d\n" % len(segment_bdy3))

    ## clean up some breakpoints
    segment_bdy4 = delete_unneeded_breakpoints2(y, xi, segment_bdy3, ll,
        delta_ll_threshold = 5, log_func=log_func)
    log_func("Segments after second round of deletion %d\n"%len(segment_bdy4))

    ## wiggle breakpoints a bit to see if there is a better positioning
    log_func("Wiggle\n")
    segment_bdy = wiggle_breakpoints(y, xi, segment_bdy4, wiggle_width=2,
        num_iterations=1, verbose=False, log_func=log_func)

    ## make sure segments completely cover genome and are non overlapping
    log_func("Validating %d segments\n"%len(segment_bdy))
    validate_segment_intervals(segment_bdy, cbdy)

    ## transform segment_bdy -> segment_index for h5 storage
    segment_index = np.zeros_like(y, dtype=bool)
    for i, (s, e) in enumerate(segment_bdy):
        segment_index[s:e] = i % 2

    ## Scale the segments to integer ploidies and return additional
    ## information about alternative scaling options

    ## aggregate bins within a segment to resolution given by window
    ## and compute segment mean read counts and lengths
    window = get_segment_window_size(y, kwargs["ll_read_threshold"])
    log_func("window = %d\n" % window)

    segment_means, segment_lengths = compute_segment_data(y, xi, segment_bdy,
        window)
    
    ## Find the scaling factor to produce integer ploidies

    ## Scaling heuristics

    # max ploidy to assign to initially chosen long segment
    max_ploidy_long=kwargs["scale_max_ploidy_long"]

    # max value of segment mean to consider "zero ploidy"
    zero_ploidy_count=kwargs["scale_zero_ploidy_count"]

    # longest segment with segment mean > zero_ploidy_count
    # that we will push to zero ploidy
    max_segment_push_to_zero=kwargs["scale_max_segment_push_to_zero"]*window

    # prior params
    prior_params = {"prior_mean"   : kwargs["scale_prior_mean"],
                    "prior_std"    : kwargs["scale_prior_std"]}

    min_ploidy = kwargs.get("min_ploidy", None)
    max_ploidy = kwargs.get("max_ploidy", None)
    
    lam_best, gap, sdf = find_best_scale_v15( y, xi, segment_bdy, segment_means,
        segment_lengths, ref, cbdy, window, scale_guess, 
        max_ploidy_long=max_ploidy_long, zero_ploidy_count=zero_ploidy_count, 
        prior_params=prior_params,
        max_segment_push_to_zero=max_segment_push_to_zero,
        min_ploidy = min_ploidy, max_ploidy = max_ploidy, verbose=True,
        log_func=log_func)

    log_func("Scaling factor: %d\n"%lam_best)

    ## avoid divide by 0
    lam_best = max(lam_best, 1e-6)

    ## Compute the ploidy vector
    ploidy = get_ploidy_vector(y, segment_means, segment_bdy, lam_best)
    log_func("Est. ploidy %.2f\n" % ploidy.mean())
    log_func("Total %.2f s elapsed\n"%(time.time( ) - t0))
    
    return ploidy, segment_index, gap, sdf, lam_best/window

def get_segment_window_size(y, ll_read_threshold):
    """ Compute the window size to aggregate segments """
    if len(y[y > 0]) == 0:
        return max_look(y)
    window = int(np.round(ll_read_threshold*1.0/np.median(y[y > 0])))
    window = np.clip(window, 1, None)
    return window

def compute_segment_data(y, xi, segment_bdy, window):
    segment_means = np.zeros(len(segment_bdy))
    segment_lengths = np.zeros(len(segment_bdy))
    for i, (s, e) in enumerate(segment_bdy):
        segment_means[i] = (y[s:e].sum()*window/xi[s:e].sum())
        segment_lengths[i] = (e-s)
    return segment_means, segment_lengths

def compute_corr_dist_to_all(y1, y2):
    dy1 = y1 - y1.mean(axis=-1)
    dy2 = y2 - y2.mean(axis=-1)[:, np.newaxis]
    cov = np.mean(dy1*dy2, axis=-1)
    sig1 = np.sqrt(np.mean(np.power(dy1, 2), axis=-1))
    sig2 = np.sqrt(np.mean(np.power(dy2, 2), axis=-1))
    return cov/(sig1*sig2)

def get_segment_bdy_from_index(sindex):
    segment_bdy2 = []
    i = 0
    while i < len(sindex):
        j = i+1
        while j < len(sindex):
            if sindex[i] != sindex[j]:
                break
            j += 1
        segment_bdy2.append((i, j))
        i = j
    return segment_bdy2

