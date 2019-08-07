#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
# check alarms produced during pipeline run
#
__MRO__ = """
stage MAKE_ALARMS(
   in  json   summary,
   in  string sample_id,
   in  string reference_path,
   out json   alarms,
   out txt    alarms_summary,
   src py     "stages/postprocessing/check_alarms",
)
"""

import tenkit.alarms
import tenkit.safe_json
from crdna.constants import CRDNA_ALARMS
import json
import martian


def main(args, outs):

    with open(args.summary) as sf:
        metrics = json.load(sf)
    
    with open(CRDNA_ALARMS) as af:
        alarms = json.load(af)

    howbadisit = tenkit.alarms.evaluate_alarms(alarms, metrics)

    with open(outs.alarms, 'w') as of:
        of.write(tenkit.safe_json.safe_jsonify(howbadisit))

    with open(outs.alarms_summary, 'w') as sf:
        sf.write("10X Genomics -- Pipeline Run Details\n")
        sf.write("-" * 40 + "\n")
        sf.write("Sample ID: %s\n" % args.sample_id)
        sf.write("Reference: %s\n" % args.reference_path)

        for oopsie in howbadisit:
            sf.write("ALERT: %s is %s. %s\n" % (oopsie["title"], oopsie["formatted_value"], oopsie["message"]))

    if len(howbadisit) > 0:
        martian.alarm("There were %i sequencing alerts. Look at alarms_summary.txt for details.\n" % len(howbadisit))


