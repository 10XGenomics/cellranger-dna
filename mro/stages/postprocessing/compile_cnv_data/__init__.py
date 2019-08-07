#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#
import h5py
import pandas as pd
import shutil
import numpy as np

from longranger.cnv import contig_manager, coverage_matrix
from crdna.clustering import get_descendant_matrix

CHUNK_SIZE=(256, 4096)
COMPRESSION = 2


__MRO__ = """
stage COMPILE_CNV_DATA(
   in  string reference_path,
   in  h5     tracks,
   in  h5     raw_profiles,
   in  h5     raw_node_profiles,
   in  h5     norm_node_profiles,
   in  h5     node_cnv_tracks,
   in  h5     tree_data,
   out h5     cnv_data,
   src py     "stages/postprocessing/compile_cnv_data",
) split using ()
"""

def split(args):
    mat_size_gb = coverage_matrix.get_genome_matrix_size_gb(
       args.raw_profiles)
    chunk_mem_gb = int(np.ceil(mat_size_gb*2))
    return {'chunks' : [{'__mem_gb': chunk_mem_gb}]}

def join(args, outs, chunk_defs, chunk_outs):
    shutil.copy(chunk_outs[0].cnv_data, outs.cnv_data)

def main(args, outs):
    compile_cnv_data(outs.cnv_data,
                     args.reference_path,
                     args.tracks,
                     args.raw_profiles,
                     args.raw_node_profiles,
                     args.norm_node_profiles,
                     args.node_cnv_tracks,
                     args.tree_data)


def compile_cnv_data(output_file,
                     reference_path, tracks,
                     cell_profiles, raw_node_profiles,
                     norm_node_profiles, node_cnv_profiles, tree_data):

    ref = contig_manager.contig_manager(reference_path)
    chroms = ref.primary_contigs(allow_sex_chromosomes=True)

    with h5py.File(output_file, "w") as outf:

        # Copy constant values
        grp = outf.create_group("constants")

        # Load cell barcodes
        with pd.HDFStore(cell_profiles, "r") as profile_store:
            constants = profile_store["constants"]

            num_cells = constants["ncells"]
            bin_size = constants["window_size"]

            outf["cell_barcodes"] = [bc.encode('utf8') for bc in profile_store["barcodes"]]
            assert(num_cells == len(profile_store["barcodes"]))

            grp["bin_size"] = int(bin_size)
            grp["num_cells"] = int(num_cells)

            num_nodes = 2*num_cells - 1
            grp["num_nodes"] = int(num_nodes)
            grp["num_chroms"] = len(chroms)
            grp["chroms"] = [c.encode('utf8') for c in chroms]

            # Compute number of bins in each chromosome
            chrom_bins = []
            for chrom in chroms:
                chrom_length = ref.contig_lengths[chrom]
                nbins = chrom_length/bin_size + int(chrom_length % bin_size!=0)
                chrom_bins.append(nbins)

            grp["num_bins_per_chrom"] = chrom_bins


        # Genome tracks - copy from tracks.h5
        with pd.HDFStore(tracks, "r") as track_store:

            # Put genome tracks in "genome_tracks" group, with harmonized naming scheme
            tracks_grp = outf.create_group("genome_tracks")
            names = [("GC", "gc_fraction"), ("map", "mappability"), ("N", "n_fraction")]

            for (old_name, new_name) in names:

                grp = tracks_grp.create_group(new_name)
                for c in chroms:
                    dataset_name = "/%s/%s" % (old_name, c)

                    data = track_store[dataset_name]
                    mat = data.as_matrix()

                    grp.create_dataset(c,
                        data = mat,
                        shape = mat.shape,
                        maxshape = (None),
                        dtype = mat.dtype,
                        compression = COMPRESSION,
                        shuffle = True,
                        chunks = (min(mat.shape[0], CHUNK_SIZE[0]),))


        # mask from raw_node_profiles -- it's a genome track
        with pd.HDFStore(raw_node_profiles, "r") as raw_store:
            grp = tracks_grp.create_group("is_mappable")
            for c in chroms:
                dataset_name = "/masks/%s" % c

                data = raw_store[dataset_name]
                mat = data.as_matrix()

                grp.create_dataset(c,
                    data = mat,
                    shape = mat.shape,
                    maxshape = (None),
                    dtype = mat.dtype,
                    compression = COMPRESSION,
                    shuffle = True,
                    chunks = (min(mat.shape[0], CHUNK_SIZE[0]),))

        # CNV call data! segment from a genome wide matrix into separate chroms 
        with pd.HDFStore(node_cnv_profiles, "r") as cnv_store:
            cnvs = cnv_store["cnv_tracks"].as_matrix()
            assert(cnvs.shape[0] == num_nodes)

            grp = outf.create_group("cnvs")

            offset = 0
            for (c, bins) in zip(chroms, chrom_bins):

                chrom_cnvs = cnvs[:, offset:(offset+bins)]

                grp.create_dataset(c,
                    data = chrom_cnvs,
                    shape = chrom_cnvs.shape,
                    maxshape = (None, None),
                    dtype = "int32",
                    compression = COMPRESSION,
                    shuffle = True,
                    chunks = CHUNK_SIZE)

                offset += bins

            assert(offset == cnvs.shape[1])


        # Normalized profiles
        with pd.HDFStore(norm_node_profiles, "r") as norm_store:
            grp = outf.create_group("normalized_counts")

            for c in chroms:
                dataset_name = "/contigs/%s" % c

                data = norm_store[dataset_name]
                mat = data.as_matrix()
                assert(mat.shape[0] == num_nodes)

                grp.create_dataset(c,
                    data = mat,
                    shape = mat.shape,
                    maxshape = (None,None),
                    dtype = "float32",
                    compression = COMPRESSION,
                    shuffle = True,
                    chunks = CHUNK_SIZE)


        # Raw profiles
        with pd.HDFStore(raw_node_profiles, "r") as raw_store:
            grp = outf.create_group("raw_counts")

            for c in chroms:
                dataset_name = "/contigs/%s" % c

                data = raw_store[dataset_name]
                mat = data.as_matrix()
                assert(mat.shape[0] == num_nodes)

                grp.create_dataset(c,
                    data = mat,
                    shape = mat.shape,
                    maxshape = (None,None),
                    dtype = "float32",
                    compression = COMPRESSION,
                    shuffle = True,
                    chunks = CHUNK_SIZE)

        # Tree data
        with pd.HDFStore(tree_data, "r") as tree_store:
            tree_grp = outf.create_group("tree")

            # Store the scipy.cluster.hierarchy tree matrix
            dataset_name = "Z"
            Z = tree_store[dataset_name].values
            tree_grp.create_dataset(dataset_name,
                data = Z,
                shape = Z.shape,
                maxshape = (None, None),
                dtype="float64",
                compression = COMPRESSION,
                shuffle = True,
                chunks = CHUNK_SIZE)

            # Store the node -> cells map as a boolean array
            desc = get_descendant_matrix(Z).astype("int32")[num_cells:, :]
            assert(desc.shape[0] == num_nodes-num_cells)

            dataset_name = "is_cell_in_group"
            tree_grp.create_dataset(dataset_name,
                data = desc,
                shape = desc.shape,
                maxshape = (None, None),
                dtype="int32",
                compression = COMPRESSION,
                shuffle = True,
                chunks = CHUNK_SIZE)

            # Store the heterogeneity matrix
            het = tree_store["het"].values
            assert(het.shape[0] == num_nodes-num_cells)

            grp = tree_grp.create_group("heterogeneity")
            offset = 0
            for (c, bins) in zip(chroms, chrom_bins):

                chrom_het = het[:, offset:(offset+bins)]

                grp.create_dataset(c,
                    data = chrom_het,
                    shape = chrom_het.shape,
                    maxshape = (None, None),
                    dtype = "float32",
                    compression = COMPRESSION,
                    shuffle = True,
                    chunks = CHUNK_SIZE)

                offset += bins

            assert(offset == het.shape[1])

