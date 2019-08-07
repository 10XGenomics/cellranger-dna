#!/usr/bin/env python
import os
import json
import sys
import websummary.summarize
import numpy as np
import pandas as pd
import numbers
import martian
import tenkit.reference

def make_sample_info(args):

    p = {
        "sample_def": args.sample_def,
        "reference_path": args.reference_path,
        "sample_id": args.sample_id,
        "sample_desc": args.sample_desc,
        "version": martian.get_pipelines_version()
    }
    return p


def fake_sample_info(outs_path):
    return {}

def get_reference_metadata(reference_path):
    path = tenkit.reference.get_metadata(reference_path)
    print "path: %s" % path
    if os.path.exists(path):
        print "found metadata"
        return json.load(open(path))

    return None


def main(args, outs):
    reference_metadata = get_reference_metadata(args.reference_path)
    sample_info = make_sample_info(args)
    body_template = os.path.join(os.path.dirname(__file__), "summary.html")
    with open(body_template, 'r') as infile:
        template = infile.read()

    data = compute_data_blob(args.summary, args.alarms, args.barnyard,
                             args.per_cell_summary_metrics, sample_info, reference_metadata)

    with open(outs.websummary, 'w') as outfile:
        websummary.summarize.generate_html_summary(data, template, None, outfile)


def operate_from_outs(path, output_file="/dev/stdout"):
    body_template = os.path.join(os.path.dirname(__file__), "summary.html")
    with open(body_template, 'r') as infile:
        template = infile.read()
    sample_info = fake_sample_info(path)

    alarms_path = os.path.join(path, "alarms.json")

    data = compute_data_blob(os.path.join(path, "summary.json"), alarms_path,  os.path.join(
        path, "barnyard.csv"), os.path.join(path, "per_cell_summary_metrics.csv"), sample_info, {})

    with open(output_file, 'w') as outfile:
        websummary.summarize.generate_html_summary(data, template, None, outfile)


def compute_data_blob(app_summary_path, app_alarms_path, barnyard_file, per_cell_summary_metrics_file, sample_info, reference_metadata):
    app_summary_info = json.load(open(app_summary_path, 'r'))

    app_alarms = json.load(open(app_alarms_path, 'r'))

    mapd_histogram = make_reads_histogram(per_cell_summary_metrics_file)
    calls_histogram = make_ploidy_histogram(per_cell_summary_metrics_file)

    bc_rank_data = make_bc_rank_data(barnyard_file)
    d = make_data_blob(app_summary_info, app_alarms, bc_rank_data,
                       mapd_histogram, calls_histogram, sample_info, reference_metadata)

    return d


def get_fastq_paths(sample_def):
    if (sample_def == None):
        return ""
    else:
        return [x["read_path"] for x in sample_def]


def getf(f, w, fmt=None):
    if (w in f):
        d = f[w]
        if (isinstance(d, numbers.Number)):
            if (fmt=='f'):
                return "%.1f%%" % (d * 100.0)
            if (fmt==","):
                return "{:,}".format(int(d))
            if (fmt == None):
                return "%.2f" % d
            else:
                return fmt.format(d)
        else:
            return "NOT A NUMBER:" + str(d)
    else:
        return "Missing data for: " + str(w)


def get(f, w):
    if (w in f):
        return f[w]
    else:
        return "Missing data for: " + str(w)


def frac(f, top, bottom):
        if (top in f and bottom in f):
            if (bottom == 0): 
                return "NaN"
            else:
                val = float(f[top])/float(f[bottom])
                return "%.1f%%" % (val*100.0)
        return "No data"


def zerobound(x, y):
    x.insert(0, x[0])
    y.insert(0, 1)
    
    x.append(x[-1])
    y.append(1)


def make_ploidy_histogram(barnyard):
    bc_data = pd.read_csv(barnyard)
    bc_calls = bc_data.mean_ploidy

    h = np.histogram(bc_calls, 25)
    y = h[0].tolist()
    x = h[1].tolist()
    return (x, y)


def make_reads_histogram(barnyard):
    bc_data = pd.read_csv(barnyard)
    total_num_reads = bc_data.total_num_reads

    h = np.histogram(total_num_reads, 25)

    y = h[0].tolist()
    x = h[1].tolist()

    return (x, y)

def halve(a):
    b=[]
    for idx in range(0, len(a)-1, 2):
        mid=(a[idx] + a[idx+1])/2
        b.append(mid)
    if (len(a) % 2 == 1):
        b.append(a[-1])
    return b


def make_bc_rank_data(barnyard):
    bc_data = pd.read_csv(barnyard)
    notcell_mask = bc_data.cell_id == "None"
    cell_mask = bc_data.cell_id != "None"

    cell_reads = bc_data.mapped[cell_mask]
    notcell_reads = bc_data.mapped[notcell_mask]

    cell_reads_depth = list(cell_reads.sort_values(ascending=False))
    notcell_reads_depth = list(notcell_reads.sort_values(ascending=False))
    cell_count = range(1, len(cell_reads_depth) + 1)
    notcell_count = range(len(cell_reads_depth) + 1,
                          len(cell_reads_depth) + len(notcell_reads_depth) + 1)

    while (len(cell_count) > 10000):
            cell_count = halve(cell_count)
            cell_reads_depth = halve(cell_reads_depth)

    while (len(notcell_count) > 10000):
        notcell_count = halve(notcell_count)
        notcell_reads_depth = halve(notcell_reads_depth)

    # connect the distributions
    notcell_count.insert(0, cell_count[-1])
    notcell_reads_depth.insert(0, cell_reads_depth[-1])

    return ((cell_count, cell_reads_depth), (notcell_count, notcell_reads_depth))



def format_quartile(app_summary_info, prefix):
    return "%s, %s, %s" % (getf(app_summary_info, prefix+"_p25"),
                           getf(app_summary_info, prefix+"_p50"),
                           getf(app_summary_info, prefix+"_p75"))


def get_metric_threshold(metric, alarm_info):
    for alarm in alarm_info:
        if alarm["id"] == metric or alarm["parent"] == metric:
            return alarm["level"].lower()
    return "pass"


def make_data_blob(app_summary_info, alarm_info, bc_rank_data, reads_per_cell_histogram_info, calls_histogram_info, sample_info, reference_metadata):

    table_invocation = { 
        "rows": [
                ["Sample ID", get(sample_info, "sample_id")],
                ["Description", get(sample_info, "sample_desc")],
                ["FASTQ path", get_fastq_paths(
                    sample_info.get("sample_def"))],
                ["Reference path", get(sample_info, "reference_path")],
                ["Cell Ranger DNA version", get(sample_info, "version")],
         ]
    }

    if reference_metadata is not None:
        rows = table_invocation["rows"]
        rows.append(["Organism", get(reference_metadata, "organism")])
        rows.append(["Assembly", get(reference_metadata, "assembly")])
        rows.append(["Annotation", get(reference_metadata, "annotation")])

    data = {
        "sample": {
            "pipeline": "Cell Ranger DNA",
            "id": get(sample_info,"sample_id"),
            "description": get(sample_info, "sample_desc")
        },
        "alarms": {"alarms": alarm_info},
        "table-invocation": table_invocation,
        "metric-cells": {"name": 'Estimated Number of Cells',
                    "metric": getf(app_summary_info, "num_cells", ","),
                    "threshold": get_metric_threshold("not_enough_cells", alarm_info)
                    },
        "metric-reads": {
            "name": 'Median effective reads per MB',
            "metric": getf(app_summary_info, "median_effective_reads_per_1Mbp", ","),
            "threshold": get_metric_threshold("reads_per_1mbp", alarm_info)
        },
        "metric-ploidy": {
            "name": "Median ploidy",
            "metric": getf(app_summary_info, "median_ploidy"),
            "threshold": "pass"
        },
        "table-seq": {
            "rows": [

                    ["Fraction of Q30 R1 Bases", 
                        frac(app_summary_info, "total_num_bases_R1_Q30", "total_num_bases_R1")],
                    ["Fraction of Q30 R2 Bases", 
                        frac(app_summary_info, "total_num_bases_R2_Q30", "total_num_bases_R2")],
                    ["Total reads", getf(app_summary_info, "total_num_reads", ",")],
                    ["Total mapped de-duplicated reads in cells", getf(
                        app_summary_info, "total_num_mapped_dedup_reads_in_cells", ",")],
                    ["Fraction of mapped de-duplicated reads in cells", frac(app_summary_info, "total_num_mapped_dedup_reads_in_cells", "total_num_reads")],
                    ["Fraction of reads with valid barcodes", getf(app_summary_info, "correct_bc_rate", "f")],
                    ["Fraction of reads not in cells", getf(app_summary_info, "frac_non_cell_barcode", "f")]
            ]
        },
        "table-cells": {
            "rows": [
                    ["Cells detected", get(app_summary_info,"num_cells")],
                    ["Median effective reads per MB", getf(app_summary_info, "median_effective_reads_per_1Mbp", ",")],
                    ["Median unmapped fraction per cell", getf(app_summary_info, "median_unmapped_frac", "f")],
                    ["Reads from cells", getf(
                        app_summary_info, "total_num_reads_in_cells", ",")],
                    ["Mean mapped de-duplicated reads per cell", getf(
                        app_summary_info, "mean_mapped_dedup_reads_per_cell", ",")],
                    ["Median duplicate fraction per cell", getf(
                        app_summary_info, "median_frac_mapped_duplicates_per_cell", "f")],
                    ["Median average ploidy", getf(app_summary_info, "median_ploidy")],
                    ["MAPD quartiles", format_quartile(app_summary_info, "normalized_mapd")],
                    ["DIMAPD quartiles", format_quartile(app_summary_info, "normalized_dimapd")],
                    ["Average ploidy quartiles", format_quartile(app_summary_info, "mean_ploidy")],
                    ["Fraction of noisy cells", getf(app_summary_info, "frac_noisy_cells", "f")]
            ]
        },

        "plot-readsPerBarcode": {
            "data": [
                {'name': "Noise",
                 'x': bc_rank_data[1][0],
                 'y': bc_rank_data[1][1],
                 'fill': 'tozeroy',
                 'line': {
                     "color": '#bdbdbd',
                     'width': 3
                 },
                 'type': "scatter",
                 },
                {'name': "Cells",
                 'x': bc_rank_data[0][0],
                 'y': bc_rank_data[0][1],
                 'fill': 'tozeroy',
                 'line': {
                     "color": '#bdedbd',
                     'width': 3
                 },
                 'type': "scatter",
                 }
            ],
            "layout": {
                'title': "Barcode Rank",
                'yaxis': {
                    'type': "log",
                    'title': "Mapped Reads",
                },
                'xaxis': {
                    "type": 'log',
                    "title": 'Barcodes',
                },
                'hovermode': 'closest'
            },
            "config": {
                'staticPlot': False
            }
        },
        "plot-readsPerCell": {
            "data": [
                {'name': "Reads-per-cell",
                 'x': reads_per_cell_histogram_info[0],
                 'y': reads_per_cell_histogram_info[1],
                 'line': {
                     "color": '#bdbdbd',
                     'width': 3
                 },
                 'hoverinfo': 'y',
                 'type': "bar",
                 }
            ],
            "layout": {
                'title': "Reads-per-cell Histogram",
                'yaxis': {
                    'type': "log",
                    'title': "Cells",
                },
                'xaxis': {
                    "type": 'linear',
                    "title": 'Reads',
                },
            },
            "config": {
                'staticPlot': False
            }
        },
        "plot-cellPloidy": {
            "data": [
                {'name': "Cell Ploidy",
                 'x': calls_histogram_info[0],
                 'y': calls_histogram_info[1],
                 'line': {
                     "color": '#bdbdbd',
                     'width': 3
                 },
                 'hoverinfo': 'y',
                 'type': "bar",
                 }
            ],
            "layout": {
                'title': "Cell Ploidy Histogram",
                'yaxis': {
                    'type': "linear",
                    'title': "Cells",
                },
                'xaxis': {
                    "type": 'linear',
                    "title": 'Average Ploidy',
                }
            },
            "config": {
                'staticPlot': False
            }
        }
    }

    return data


if __name__ == "__main__":
    operate_from_outs(sys.argv[1])
