#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

import numpy as np
import pandas as pd
import os
import martian
from collections import defaultdict
import crdna.constants

__MRO__ = """
stage AGGREGATE_TRACKS(
    in  bed      targets,
    in  bed      confident_windows,
    in  h5       map_track,
    in  h5       genome_tracks,
    out h5       tracks,
    src py       "stages/genome_annotator/aggregate_tracks",
)
"""

def split(args):
    raise Exception("Split is unimplemented")

def main(args, outs):
    in_h5_list = [args.map_track, args.genome_tracks]
    out = pd.HDFStore(outs.tracks, "w")
    
    ## first add in mappability and GC tracks
    msizes = {}
    for in_h5 in in_h5_list:
        if in_h5 is None:
            continue
        if not os.path.exists(in_h5):
            martian.exit("Could not find " + in_h5)
            continue
        indata = pd.HDFStore(in_h5, "r")
        for key in indata.keys( ):
            X = indata[key]
            #
            # column-wise concatenation for the constants section
            #
            if key in out:
                out[key] = pd.concat( [out[key],X], axis=0 )
            else:
                out[key] = X 
            chrom = key.split("/")[-1]
            msizes[chrom] = X.shape[0]
        indata.close( )

    ## convert confident windows into numpy and store
    if args.confident_windows is None or not os.path.exists(args.confident_windows):
        for chrom, length in msizes.iteritems():
            out["/CONF/"+chrom] = pd.Series(np.ones(length, dtype=float))
    else:
        conf = defaultdict(list)
        for line in open(args.confident_windows):
            fields = line.strip().split()
            chrom, perc = fields[0], float(fields[3])
            conf[chrom].append(perc)
        for chrom, length in msizes.iteritems():
            X = np.array(conf[chrom])
            cbins = (X > crdna.constants.CONFIDENT_BIN_THRESHOLD).sum( )
            if X.shape[0] == 0 or cbins == 0:
                out["/CONF/"+chrom] = pd.Series(np.ones(length, dtype=float))
            else:
                assert X.shape[0] == length
                out["/CONF/"+chrom] = pd.Series(X)
    out.close( )

def join(args, outs, chunk_defs, chunk_outs):
    raise Exception("Join is unimplemented")

