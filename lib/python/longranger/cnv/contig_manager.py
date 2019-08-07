#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#

import json
import os.path


def check_aos(a):
    if not isinstance(a, list):
        return False

    for k in a:
        if not isinstance(k, unicode):
            return False
    return True

# This function verifies that a contig definitions file is correctly formed
def verify_contig_defs(contig_defs_path, fasta_index_path):

    # Step 1 does the file exist?
    if not os.path.exists(contig_defs_path):
        return "Contig definitions file '%s' is missing" % (contig_defs_path)

    # Step 2 is it valid JSON?
    try: 
        contigs = json.load(open(contig_defs_path))
    except ValueError,e:
        return "Malformed json in '%s': %s" % (contig_defs_path, e)

    # Step 3: species_prefixes must be an array of strings
    species_prefixes = contigs.get("species_prefixes", [u""])
    if not check_aos(species_prefixes):
        return "Species_prefixes must be an array of strings"

    # Step 4: primary contigs must be an array of strings 
    if not check_aos(contigs.get("primary_contigs")):
        return "Primary_contigs must be an array of strings"

    # Step 5: prefix contigs can not be prefixes of themselves
    for p1 in species_prefixes:
        for p2 in species_prefixes:
            if p1 != p2 and p1.startswith(p2):
                return "Species_prefixes are ambiguous. No prefix may be a prefix of another prefix."

    # step 6: every primary contig must be prefixed by a species prefix
    for c in contigs["primary_contigs"]:
        ok = False
        for p in species_prefixes:
            if c.startswith(p):
                ok = True
        if not ok:
            return "Each primary contig must be prefixed by a species prefix"

    # Step 7: sex_chromosomes must be a map of maps; Each sub-maps keys must be a primary contig; each sub-maps values must be an integer
    if contigs["sex_chromosomes"] is None or not isinstance(contigs["sex_chromosomes"], dict):
        return "Sex chromosomes must be an object of objects of integers."

    for k in contigs["sex_chromosomes"]:
        v = contigs["sex_chromosomes"][k]
        if not isinstance(v, dict):
            return "Sex chromosomes must be an object of objects of integers. {} is not an object.".format(k)
        for sx in v:
            if not sx in contigs["primary_contigs"]:
                return "Sex chromosomes must be primary contigs. %s is not a primary contig." % sx
            if not isinstance(v[sx], int):
                return "Sex chromosome ploidies must be integers."


    # Step 8: every primary contig must exist in the reference
    all_fasta_contigs = []

    for line in open(fasta_index_path):
        fields = line.strip().split()
        all_fasta_contigs.append(fields[0])

    for c in contigs["primary_contigs"]:
        if not c in (all_fasta_contigs):
            return "Primary contig %s is in the %s but not in the reference" % (c, contig_defs_path)

    # Step 9: there must be a primary contig
    if len(contigs["primary_contigs"]) == 0:
        return "At least one contig must be primary."

    return None

# This implements reference "metadata" operations for single species and barnyard references.
# The names for the contigs in the reference are formated as either: "chr<XXX>" or "<species>_chr<XXX>".
# This operates on a file, "contigs.json" found in the reference fasta directory that describes
# the layout.


class contig_manager:
    def __init__(self, path):
        contigs_def_path = os.path.join(path, 'fasta', 'contig-defs.json')
        fasta_index_path = os.path.join(path, "fasta", "genome.fa.fai")
        assert os.path.exists(contigs_def_path)
        assert os.path.exists(fasta_index_path)

        self.contigs = json.load(open(contigs_def_path))
        ## add a species_prefixes key if not present
        if "species_prefixes" not in self.contigs:
            self.contigs["species_prefixes"] = [""]

        self.sex_chromosomes = {}
        for s in self.contigs['sex_chromosomes']:
            for c in self.contigs['sex_chromosomes'][s]:
                self.sex_chromosomes[c] = True

        all_contigs    = []
        contig_lengths = {}
        for line in open(fasta_index_path):
            fields = line.strip().split()
            all_contigs.append(fields[0])
            contig_lengths[fields[0]] = int(fields[1])
        self.all_contigs    = all_contigs
        self.contig_lengths = contig_lengths
    
    def list_all_contigs(self, species=None):
        """List all contigs, optionally filter by species."""
        if (species != None):
            if not species in self.contigs['species_prefixes']:
                raise Exception('Unknown species')
        return [x for x in self.all_contigs if (species is None or self.species_from_contig(x) == species)]
    
    def get_contig_lengths(self):
        """Return dictionary of contig: contig length"""
        return self.contig_lengths

    def primary_contigs(self, species=None, allow_sex_chromosomes=False):
        """List all primary contigs, optionally filter by species."""
        if (species != None):
            if not species in self.contigs['species_prefixes']:
                raise Exception('Unknown species')
        return [x for x in self.contigs['primary_contigs'] if self.is_primary_contig(x, species=species, allow_sex_chromosomes=allow_sex_chromosomes)]

    def expected_ploidy(self, contig, sex):
        """Return the expected ploidy for a given contig and a particular
        sex."""
        if not contig in self.all_contigs:
            raise Exception('Unknown contig')
        else:
            if sex is None:
                return 2
            species = self.species_from_contig(contig)
            key = species + '_' + sex
            if not key in self.contigs['sex_chromosomes']:
                raise Exception('Unknown species or sex')
            if (contig in self.contigs['sex_chromosomes'][key]):
                return self.contigs['sex_chromosomes'][key][contig]
            else:
                return 2

    def is_primary_contig(self, contig_name, species=None, allow_sex_chromosomes=False):
        if not contig_name in self.contigs['primary_contigs']:
            return False

        if (species != None) and self.species_from_contig(contig_name) != species:
            return False

        if (allow_sex_chromosomes == False):
            if (contig_name in self.sex_chromosomes):
                return False
        return True

    def list_species(self):
        """List all speces contained in this reference."""
        return self.contigs['species_prefixes']

    def species_from_contig(self, contig):
        """Given the name of a contig, extract the species."""
        if self.contigs['species_prefixes'] == [""]:
            return ""

        s = contig.split('_')
        if len(s) > 1:
            return s[0]
        else:
            return ''

    def non_nuclear_contigs(self):
        """ Lists non-nuclear contigs (e.g. mitochondria and chloroplasts) """
        if self.contigs.has_key("non_nuclear_contigs"):
            return self.contigs["non_nuclear_contigs"]
        else:
            None
