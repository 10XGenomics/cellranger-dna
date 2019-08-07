# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
__doc__="""These are constants that are specific to the cellranger-DNA pipeline.  Other
general constants can be found in tenkit.constants."""

import os

## Heuristics used by pipeline

## Mapping quality thresholds

# what we allow into read profiles
PROFILE_MAPQ_THRESHOLD     = 30

# what is used in cell detection
CELL_DETECT_MAPQ_THRESHOLD = 30

# what is the longest insert to consider the read still part of a pair
PROFILE_INSERT_THRESHOLD   = 20000

# definition of a confident bin
CONFIDENT_BIN_THRESHOLD    = 0.70

# what defines a good bin @ 20kb
MAPPABILITY_THRESHOLD      = 0.90

# fixed amplicon length used for amp_rate calculation
DEFAULT_AMPLICON_LENGTH = 250

# VESNA settings
MAX_PLOIDY = 12

# N-base cutoff: bins with higher N-base content are filtered out
N_BASE_THRESHOLD = 0.30

## CNV calling thresholds
MIN_CLUSTER_SIZE_TO_CALL = 3

## Heuristics for GC
GC_RES                = 1     ## in multiples of 20 kb
MIN_GC                = 0.30
MAX_GC                = 0.65
NUM_GC_BINS           = 20
MIN_POINTS_PER_BIN    = 30    ## only take bins that have at least this many pts
GC_ORIGIN             = 0.45  ## origin for GC axis

## CNV calling heuristics
BREAKPOINT_READ_THRESHOLD  = 200.0
BREAKPOINT_CALLER_HEURISTICS = {
    "ll_read_threshold": BREAKPOINT_READ_THRESHOLD,
    "ll_min_look": 5,
    "ll_threshold": 5,
    "delta_threshold": 0.0,
    "del_bp_ll_threshold": 5.0,
    "add_bp_min_segment_size": 5,
    "add_bp_min_ll": 5.0,
    "scale_max_ploidy_long": 10,
    "scale_zero_ploidy_count": BREAKPOINT_READ_THRESHOLD/4.0,
    "scale_max_segment_push_to_zero": 5,
    "scale_prior_mean" : 2.0,
    "scale_prior_std" : 1.0
}

# these override the default values from tenkit
RAW_BARCODE_TAG = 'CR'
PROCESSED_BARCODE_TAG='CB'
RAW_BARCODE_QUAL_TAG = 'CY'

# Additional tags for the 5' global position of read and its mate. Used for sorting and marking duplicates
SELF_FIVE_PRIME_POS_TAG = 'GP'
MATE_FIVE_PRIME_POS_TAG = 'MP'
DUPLICATE_COUNT_TAG = 'DC'

# used by REPORT_BASIC
MAX_READ_LENGTH = 300

# maximum supported number of cells
MAX_CELLS = 20000

# Alarm JSON file defines alarms passed to user at completion
CRDNA_ALARMS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "alarms", "crdna.json")

# GC normalization curve limits
MIN_CURVE = 0.1
MAX_CURVE = 10.0
