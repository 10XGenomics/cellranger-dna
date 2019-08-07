#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# Mark PCR duplicates in a BAM file
#
import json
import tenkit.dict_utils
import tenkit.bam as tk_bam
import tenkit.lane as tk_lane
import tenkit.coverage
import numpy as np
from crdna.constants import SELF_FIVE_PRIME_POS_TAG

from crdna.duplicates import estimate_flowcell_geometry, compute_null_distribution, DupSummary, broadcast, default_diffusion_threshold
from tenkit.stats import robust_divide

__MRO__ = """
stage DIFFUSION_DUP_ESTIMATOR(
    in  bam     input,
    out json    summary,
    src py      "stages/reads/diffusion_dup_estimator",
) split using (
    in  int     seed,
    in  map     lane_map,
    in  map     flowcell_geometry,
    in  string  chunk_start,
    in  string  chunk_end,
)
"""

MIN_SAMPLES = 500000

def chunk_bound_func(read):
    # Since the reads are sorted by SELF_FIVE_PRIME_POS_TAG tag, use that
    # for the chunk boundary
    if not read.is_unmapped:
        return (read.reference_id, read.get_tag(SELF_FIVE_PRIME_POS_TAG))
    else:
        return None



def split(args):

    # Chunk bam to get 1GB per chunk
    bam_in = tk_bam.create_bam_infile(args.input)
    lane_coord_sys = tk_lane.LaneCoordinateSystem()

    bam_in.reset()
    lane_coord_sys.estimate_tile_extents(bam_in)
    flowcell_geometry = estimate_flowcell_geometry(bam_in, lane_coord_sys)

    print "Flowcell Geometry: ", flowcell_geometry

    if flowcell_geometry is None:
        return {'chunks': [{'seed': None, 'lane_map': None, 'flowcell_geometry': None, 'chunk_start': None, 'chunk_end': None}] }

    chunk_defs = tk_bam.chunk_bam_records(bam_in, chunk_bound_func, chunk_size_gb=0.75)

    for i, chunk in enumerate(chunk_defs):
        chunk['seed'] = i
        chunk['__mem_gb'] = 3

    for chunk in chunk_defs:
        chunk['lane_map'] = lane_coord_sys.to_dict()
        chunk['flowcell_geometry'] = flowcell_geometry

    return {'chunks': chunk_defs, 'join': {'__mem_gb': 1}}

def main(args, outs):

    if args.flowcell_geometry is None:
        return

    lane_coord_sys = tk_lane.LaneCoordinateSystem.from_dict(args.lane_map)

    args.coerce_strings()
    outs.coerce_strings()

    bam_in = tk_bam.create_bam_infile(args.input)
    
    null_distribution = compute_null_distribution(args.flowcell_geometry, seed=args.seed)

    estimator = DupSummary.diffusion_estimator(lane_coord_sys, args.flowcell_geometry)

    consumers = [estimator.read_consumer()]

    source = tk_bam.read_bam_chunk(bam_in, (args.chunk_start, args.chunk_end))
    broadcast(source, consumers)

    # Package up the summaries:
    dup_results = {'null_distribution': null_distribution, "observed_distribution": estimator.result}

    if outs.summary:
        with open(outs.summary, 'w') as f:
            json.dump(dup_results, f, indent=4)

def distribution_from_dict(d):
    
    x = [float(x) for x in d.keys()]
    factor = sum(d.values())
    y = [float(y)/factor for y in d.values()]
    
    X = np.asarray([_x for _x, _y in sorted(zip(x, y))])
    Y = np.asarray([_y for _x, _y in sorted(zip(x, y))])
    return X, Y

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    if chunk_defs[0].flowcell_geometry is None:
        combined_dups = {}
        combined_dups['diffusion'] = {}
        combined_dups['diffusion']['threshold'] =  default_diffusion_threshold(None)

        with open(outs.summary, 'w') as f:
            json.dump(combined_dups, f, indent=4)
        return

    # combine the summary
    dup_summaries = [json.load(open(out.summary)) for out in chunk_outs]
    combined_dups = reduce(lambda x,y: tenkit.dict_utils.add_dicts(x,y,2), dup_summaries)

    null_stat = {}
    for k, v in combined_dups['null_distribution'].items():
        null_stat[float(k)] = v

    cutoff = 0.1 * max(null_stat.keys())
    null_samples = sum(null_stat.values())
    null_success = sum([v for k, v in null_stat.items() if k > cutoff])

    assert null_samples > 0
    assert null_success > 0

    p_null = robust_divide(null_success, null_samples)
    print "P_null = ", p_null


    obs_stat = {}
    for k, v in combined_dups['observed_distribution'].items():
        obs_stat[float(k)] = v
    obs_samples = sum(obs_stat.values())
    obs_success = sum([v for k, v in obs_stat.items() if k > cutoff])

    if obs_samples == 0:
        p_obs = 0.0
    else:
        p_obs = robust_divide(obs_success, obs_samples)

    print "P_obs = ", p_obs

    pcr_dup_fraction = min(1.0, robust_divide(p_obs, p_null))
    diffusion_dup_fraction = 1.0 - pcr_dup_fraction
    print "PCR dup fraction = ", pcr_dup_fraction

    # Compute threshold
    
    # Handle edge case
    diffusion_threshold = default_diffusion_threshold(chunk_defs[0].flowcell_geometry)
    print "Default threshold = ", diffusion_threshold
    if obs_samples > MIN_SAMPLES:

        x_null, y_null = distribution_from_dict(combined_dups['null_distribution'])
        x_obs, y_obs = distribution_from_dict(combined_dups['observed_distribution'])
        y_diff = y_obs - pcr_dup_fraction * y_null

        n_bins = len(x_null)

        # Clamp at zero
        for i in xrange(n_bins):
            y_diff[i] = max(y_diff[i], 0.)
            if i >= 1 and y_diff[i-1]==0.:
                y_diff[i] = 0.

        for i in range(n_bins):
            if y_diff[i] < pcr_dup_fraction*y_null[i]:
                crossover = x_null[i]
                break
        print "Crossover = ", crossover

        factor = sum(y_diff)
        cumulative = 0.0
        for i in xrange(n_bins):
            cumulative = cumulative + y_diff[i]/factor
            if cumulative > 0.99:
                cutoff_99 = x_null[i]
                break
        print "Cutoff(99%) = ", cutoff_99
        
        if crossover == 0.0:
            print "WARNING: There was no crossover. Switching to default threshold"
        else:
            diffusion_threshold = cutoff_99
    else:
        print "WARNING: Not enough observed samples, using the default threshold"

    print "Diffusion threshold = ", diffusion_threshold
    combined_dups['dup_fraction'] = {}
    combined_dups['dup_fraction']['pcr'] = pcr_dup_fraction
    combined_dups['dup_fraction']['diffusion'] = diffusion_dup_fraction
    combined_dups['diffusion'] = {}
    combined_dups['diffusion']['threshold'] = diffusion_threshold

    with open(outs.summary, 'w') as f:
        json.dump(combined_dups, f, indent=4)
