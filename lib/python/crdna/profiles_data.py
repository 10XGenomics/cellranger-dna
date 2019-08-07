#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
__doc__="""Provides classes for reading and manipulating data from profiles.h5 file."""
import os
import numpy as np
import pandas as pd

import crdna.constants
import longranger.cnv.contig_manager as contig_manager
from longranger.cnv.coverage_matrix import load_matrix, list_all_contigs
from longranger.cnv.copy_number_tools import aggregate_counts, aggregate_matrix

class BaseProfilesData(object):
    """Base class for reading a profiles.h5 file into memory and facilitating
    operations on the underlying count matrices.
    """
    def __init__(self, filename, reference_path, load_conf_filter=True, reuse_mask_from=None ):
        """Constructs a ProfilesData object from a profiles.h5 file.

        :param filename: HDF5 filename
        :param reference_path: reference path
        :param load_conf_filter: (optional) whether to auto-discover a confident bins track (tracks.h5)
        :param reuse_mask_from: (optional) use mask in the supplied ProfilesData object instead of the mask
        found in the HDF5 file.
        """
        self._filename = filename
        self._reference_path = reference_path
        self._rebinned = False
        if not os.path.exists(filename):
            raise IOError("Can't find file %s"%filename)
        self._load_hdf5_file(filename,reference_path,reuse_mask_from)
        if load_conf_filter:
            self.load_conf_filter()
        self.reset_mask()

    def load_conf_filter(self,tracks=None):
        if tracks is None:
            tracks = _find_tracks_file(self._filename)
        if tracks is not None:
            self._conf_filter = ConfFilter(tracks,self._reference_path)
        else:
            self._conf_filter = None
        return tracks is not None

    def _load_hdf5_file(self, filename, reference_path, reuse_mask_from=None):
        X, Xm = load_matrix(filename,reference_path)
        if len(X)==0:
            raise ValueError( "Loading profiles from %s returned zero contigs matching the reference at %s"%(filename,reference_path) )
        # it's low cost to hang on to a contig_manager instance
        self._contig_manager = contig_manager.contig_manager(reference_path)
        # this list has to match the list inside load_matrix
        primary_contigs = self._contig_manager.primary_contigs(allow_sex_chromosomes=True)
        store = pd.HDFStore(filename, "r");
        self._window_size = store['constants']['window_size']
        if 'barcodes' in store:
            self._barcodes = np.array(store['barcodes'])
            self._num_cells = len(self._barcodes)
        else:
            self._barcodes = None
            self._num_cells = len(X[0])
        profile_contigs = set(list_all_contigs(store))
        i = 0
        self._contig_list = []
        self._contig_coverage = {}
        self._contig_mask = {}
        self._contig_idx = {}
        for chrom in primary_contigs:
            if chrom in profile_contigs:
                self._contig_list.append(chrom)
                self._contig_coverage[chrom] = X[i]
                if reuse_mask_from:
                    self._contig_mask[chrom] = reuse_mask_from._contig_mask[chrom]
                else:
                    self._contig_mask[chrom] = Xm[i]
                nbins = X[i].shape[1]
                self._contig_idx[chrom] = np.arange(1,nbins*self._window_size,self._window_size)
                i = i + 1
        store.close()

    def reset_mask(self):
        """Resets the mask on this data (as if no mask was applied).
        IMPORTANT: call this before aggregating data."""
        if self._rebinned:
            raise Exception('Resetting mask after re-binning is not supported.  Initialize a new ProfilesData object')
        self._masks = { c : np.ones(len(self._contig_mask[c]),dtype='bool') for c in self._contig_list }

    def apply_mask(self,use_default_mask=True,use_conf_filter=True):
        """Calculates the proper masking to apply to these profiles, allowing
        for optional masks.  Masking is calculated at the original window size
        of the profiles data.

        IMPORTANT: call this before aggregating.

        :param use_default_mask: apply the mask found in profiles.h5
        :param use_conf_filter: apply a confident regions filter
        :returns: tuple of total number of bins, total number of unmasked bins
        """
        if self._rebinned:
            raise Exception('Recalculating mask after re-binning is not supported.  Initialize a new ProfilesData object')
        self._masks = {}
        tot_unmasked, tot_bins = 0, 0
        for chrom in self._contig_list:
            nbins = len(self._contig_mask[chrom])
            mask = np.ones(nbins,dtype='bool')
            if use_default_mask:
                f = self._contig_mask[chrom]
                mask = mask & f
            if use_conf_filter and self._conf_filter:
                f = self._conf_filter.get_filter_for_contig(chrom)
                mask = mask & f
            self._masks[chrom] = mask
            tot_unmasked += mask.sum()
            tot_bins += nbins
        return tot_bins, tot_unmasked

    def get_coverage_for_contig(self,contig):
        """Return the count matrix (ncells x nbins) for a single contig with position
        information.

        Applies any masks which have been set on this object.

        :returns: tuple of contig names (nbins nparray), positions (nbins nparray), 
                  counts (ncells x nbins nparray)
        """ 
        if not contig in self._contig_list:
            return None
        mask = self._masks[contig]
        nbins = mask.sum()
        G = self._contig_coverage[contig][:,mask]
        contigs = np.repeat([contig],nbins)
        pos = self._contig_idx[contig][mask]
        return contigs, pos, G

    def get_stitched_coverage(self,allow_sex_chromosomes=True):
        """Returns stitched together view of coverage across all contigs.

        :param allow_sex_chromosomes: if True, include sex chromosomes
        :returns: tuple of contig names (nbins nparray), positions (nbins nparray), 
                  counts (ncells x nbins nparray)
        """
        contig_list = [ c for c in self._contig_list \
            if self._contig_manager.is_primary_contig(c,allow_sex_chromosomes=allow_sex_chromosomes) ]
        total_bins = sum( self._masks[c].sum() for c in contig_list )
        G = np.zeros( (self._num_cells,total_bins) )
        contigs = np.zeros(total_bins,dtype='object')
        pos = np.zeros(total_bins,dtype='int32')
        i = 0
        for chrom in contig_list:
            _contigs, _pos, _G = self.get_coverage_for_contig(chrom)
            nbins = len(_contigs)
            G[:,i:i+nbins] = _G
            contigs[i:i+nbins] = _contigs
            pos[i:i+nbins] = _pos
            i += nbins
        return contigs, pos, G

    def get_num_cells(self):
        return self._num_cells

    def get_barcodes(self):
        return self._barcodes

    def get_window_size(self):
        return self._window_size

    def get_contigs(self):
        return self._contig_list

class ProfilesData2(BaseProfilesData):
    """Class for reading a profiles.h5 file into memory and facilitating
    operations on the underlying count matrices.

    This class handles masking+aggregation by masking at the finest level and
    then aggregating.  The result is variable bin sizes.

    Example:
        data = ProfilesData2('profiles.h5',reference_path)
        print data.get_contigs()
        data.apply_mask(use_default_mask=True,use_conf_filter=False)
        data.aggregate(factor=5)
        contigs, pos, coverage = data.get_stitched_coverage(allow_sex_chromosomes=False)
    """
    def __init__(self, filename, reference_path, load_conf_filter=True, reuse_mask_from=None ):
        super(ProfilesData2,self).__init__(filename,reference_path,load_conf_filter,reuse_mask_from)

    def aggregate(self, factor=5, mask_factor=None):
        """This is a destructive operation.  Reduce all of the data in this object by
        rebinning by an integer factor.
        
        :param factor: reduce the data by this factor
        :param mask_factor: not used
        """
        self._contig_coverage = { \
              chrom : aggregate_matrix(self._contig_coverage[chrom][:,self._masks[chrom]], factor) \
              for chrom in self._contig_list }
        for chrom in self._contig_list:
            self._contig_idx[chrom] = self._contig_idx[chrom][self._masks[chrom]][::factor]
        self._masks = { \
              c : np.ones(self._contig_coverage[c].shape[1],dtype='bool') \
              for c in self._contig_list }
        self._window_size = self._window_size * factor
        self._rebinned = True

class ProfilesData(BaseProfilesData):
    """Class for reading a profiles.h5 file into memory and facilitating
    operations on the underlying count matrices.

    This class handles masking+aggregation by masking at the super-bin level.
    For larger aggregation factors, this behavior is sub-optimal, and you should
    probably use ProfilesData2.

    Example:
        data = ProfilesData('profiles.h5',reference_path)
        print data.get_contigs()
        data.apply_mask(use_default_mask=True,use_conf_filter=False)
        data.aggregate(factor=5)
        contigs, pos, coverage = data.get_stitched_coverage(allow_sex_chromosomes=False)
    """
    def __init__(self, filename, reference_path, load_conf_filter=True, reuse_mask_from=None ):
        super(ProfilesData,self).__init__(filename,reference_path,load_conf_filter,reuse_mask_from)

    def aggregate(self, factor=5, mask_factor=0.6):
        """This is a destructive operation.  Reduce all of the data in this object by
        rebinning by an integer factor.
        
        :param factor: reduce the data by this factor
        :param mask_factor: when aggregating the mask, this percentage of subbins need to be
                            True to result in a True superbin (default=0.6)
        """
        self._contig_coverage = { \
              chrom : aggregate_matrix(self._contig_coverage[chrom], factor) \
              for chrom in self._contig_list }
        mask_count = int(mask_factor*factor)
        self._masks = { \
              chrom : aggregate_counts(self._masks[chrom].astype(int), factor)>=mask_count \
              for chrom in self._contig_list }
        self._window_size = self._window_size * factor
        for chrom in self._contig_list:
            nbins = self._contig_coverage[chrom].shape[1]
            self._contig_idx[chrom] = np.arange(1,nbins*self._window_size,self._window_size)
        self._rebinned = True

class ConfFilter(object):
    """Handles reading the CONF track from the tracks.h5 file."""
    #...........................................................................
    def __init__(self, filename, reference_path):
        self._filename = filename
        if not os.path.exists(filename):
            raise IOError("Can't find file %s"%filename)
        # if filename
        self._load_hdf5_file(filename,reference_path)
    # __init__

    #...........................................................................
    def _load_hdf5_file(self, filename, reference_path):
        self._contig_manager = contig_manager.contig_manager(reference_path)
        self._contig_list = self._contig_manager.primary_contigs(allow_sex_chromosomes=True)
        store = pd.HDFStore(filename, "r")
        self._window_size = store['constants']['window_size']
        self._conf_filter = {}
        for chrom in self._contig_list:
            cmask = (store["/CONF/"+chrom].values > crdna.constants.CONFIDENT_BIN_THRESHOLD) & \
                    (store["/N/"+chrom].values < 1.0/self._window_size)
            self._conf_filter[chrom] = cmask
        # for chrom
        store.close()
    # _load_hdf5_file

    #...........................................................................
    def get_filter_for_contig(self, contig):
        return self._conf_filter[contig]
    # get_filter_for_contig

    #...........................................................................
    def aggregate(self, factor=5):
        """This is a destructive operation.  Reduce all of the data in this object by
        rebinning by an integer factor."""
        #
        # the logic here is very strict: if any of the sub-bins are not confident then
        # the new superbin is not confident
        #
        self._conf_filter = { \
              chrom : aggregate_counts(self._conf_filter[chrom].astype(int), factor)==factor \
              for chrom in self._contig_list }
        self._window_size = factor * self._window_size
    # aggregate

    #...........................................................................
    def get_stitched_conf_filter(self, allow_sex_chromosomes=True):
        contig_list = [ c for c in self._contig_list \
            if self._contig_manager.is_primary_contig(c,allow_sex_chromosomes=allow_sex_chromosomes) ]
        total_bins = sum( self._conf_filter[c].sum() for c in contig_list )
        conf_filter = np.ones(total_bins, dtype='bool')
        i = 0
        for chrom in contig_list:
            nbins = len(self._conf_filter[chrom])
            conf_filter[ i:i+nbins ] = self._conf_filter[chrom]
            i += nbins
        # for chrom
        return conf_filter
    # get_stitched_conf_filter

# ConfFilter

#...............................................................................
def _find_tracks_file(profiles_h5):
    "Logic for autodiscovering a tracks.h5 file.  Useful for prod dev work."
    tracks = os.path.join(os.path.dirname(profiles_h5),'tracks.h5')
    if os.path.exists(tracks):
        return tracks
    # if tracks
    #
    # fall back to the old location of this file
    #
    tracks = os.path.join(os.path.dirname(profiles_h5),
        '../ANALYZER_SINGLECELL_PD/CNV_CALLER_SINGLECELL/_REPORTER_SINGLECELL/_GENOME_ANNOTATOR/AGGREGATE_TRACKS/fork0/files/',
        'tracks.h5')
    if os.path.exists(tracks):
        return tracks
    # if tracks
    return None
# _find_tracks_file

