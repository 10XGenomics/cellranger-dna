import pandas as pd
import contig_manager

#..............................................................................
def list_all_contigs(store):
    #store = pandas.HDFStore("matrix_path", "r");
    contigs = store.get_node("/contigs")._v_children.keys()
    return contigs
# list_all_contigs

#..............................................................................
def contig_to_coverage_path(contig):
    return "/contigs/" + contig
# contig_to_coverage_path

#..............................................................................
def contig_to_mask_path(contig):
    return "/masks/" + contig
# contig_to_mask_path

#..............................................................................
def contig_to_coverage(store, contig):
    return store[contig_to_coverage_path(contig)]
# contig_to_coverage

#..............................................................................
def contig_to_mask(store, contig):
    return store[contig_to_mask_path(contig)]
# contig_to_mask

#..............................................................................
def list_primary_contigs(file_name, reference_path, allow_sex_chromosomes=True):
    store = pd.HDFStore(file_name, "r")
    profile_contigs = set(list_all_contigs(store))
    store.close()
    ref = contig_manager.contig_manager(reference_path)
    primary_contigs = ref.primary_contigs(
        species=None, allow_sex_chromosomes=allow_sex_chromosomes)
    filtered = [x for x in primary_contigs if x in profile_contigs]
    return filtered
# list_primary_contigs

#..............................................................................
def load_matrix(file_name, reference_path, start_cell=None, end_cell=None):
    assert (start_cell is None) == (end_cell is None), "specify both "\
        "start_cell and end_cell"
    if start_cell is not None:
        assert start_cell < end_cell, "start_cell < end_cell"
    raw_profiles = []
    mask = []
    reference = contig_manager.contig_manager(reference_path)
    primary_contigs = reference.primary_contigs(allow_sex_chromosomes=True)
    store = pd.HDFStore(file_name, "r");
    for chrom in primary_contigs:
        chrom_profile = contig_to_coverage(store, chrom).values
        chrom_mask = contig_to_mask(store, chrom).values
        raw_profiles.append(chrom_profile[start_cell:end_cell, :].copy())
        mask.append(chrom_mask)
    store.close()
    return (raw_profiles, mask)
# load_matrix

#..............................................................................
def store_matrix(file_name, chroms, profiles, tracks, window_size, masks=None,
    dtype="float"):
    """Creates a profiles.h5 file.
    :param file_name: path to the new HDF5 file
    :param chroms: list of contig names
    :param profiles: list of per-contig profiles
    :param tracks: tracks HDFStore, can be None if masks is set
    :param window_size: window size for the profiles
    :param masks: If masks is set, use this for creating the mask track, not tracks
    """
    store = pd.HDFStore(file_name, "w")
    for chrom_index, chrom in enumerate(chroms):
        store["/contigs/" + chrom] = pd.DataFrame(profiles[chrom_index].astype(dtype))
        if not masks:
            ## number of Ns in each 20 kb bin
            N_content = tracks["/N/" + chrom]
            store["/masks/" + chrom] = pd.Series(N_content < 1.0 / window_size)
        else:
            store["/masks/" + chrom] = pd.Series(masks[chrom_index])
    # for chrom
    ## store the window size in the h5
    store["constants"] = pd.Series({"window_size": window_size})
    store.close()
# store_matrix

#..............................................................................
def get_bin_size(file_name):
    """ 
        Find the bin size associated with a particular h5 file.
    """
    store = pd.HDFStore(file_name, "r")
    assert "constants" in store, "No 'constants' found in %s"%file_name
    bin_size = store["constants"]["window_size"]
    store.close( )
    return bin_size
# get_bin_size

#..............................................................................
def get_num_cells(file_name, reference_path):
    """
        Get the number of cells/number of rows in the profiles h5 data
    """
    ref = contig_manager.contig_manager(reference_path)
    chrom = ref.primary_contigs( )[0]
    store = pd.HDFStore(file_name, "r")
    ncells = store["/contigs/"+chrom].shape[0]
    store.close( )
    return ncells
# get_num_cells

def get_genome_matrix_size_gb(profiles):
    """ Returns the size in GB of the full # cells x # bins matrix """
    store = pd.HDFStore(profiles, "r")
    constants = store["constants"]
    store.close()
    nbins = constants["genomebins"]
    ncells = constants["ncells"]
    mat_size = 4*ncells*nbins/1e9
    return mat_size


