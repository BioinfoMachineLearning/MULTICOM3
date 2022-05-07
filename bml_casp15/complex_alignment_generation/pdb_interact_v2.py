import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *
import os


class PDB_interact_v2:

    def __init__(self, uniprot2pdb_mapping_file, complexes_cm_file):
        self.uniprot2pdb_mapping_file = uniprot2pdb_mapping_file
        self.complexes_cm_file = complexes_cm_file
        self.uniprot2pdb_map = {}

    def load_data(self):
        with open(self.uniprot2pdb_mapping_file) as f:
            for line in f:
                # UniRef100_P13744        2E9QA
                uniprot_id, pdbcode = line.rstrip('\n').split()
                uniprot_id = uniprot_id.split('_')[1]
                if uniprot_id in self.uniprot2pdb_map:
                    self.uniprot2pdb_map[uniprot_id] += "," + pdbcode
                else:
                    self.uniprot2pdb_map[uniprot_id] = pdbcode

    def get_interactions_from_ids(self, ids, is_homomers):
        pdbcodes = self.uniprot2pdb_map[ids[0]]
        prev_df = pd.DataFrame({'pdbcode': [pdbcode[0:4] for pdbcode in pdbcodes]})
        for i in range(1, len(ids)):
            pdbcodes = self.uniprot2pdb_map[ids[i]]
            curr_df = pd.DataFrame({'pdbcode': [pdbcode[0:4] for pdbcode in pdbcodes]})
            prev_df = prev_df.merge(curr_df, on="pdbcode")

        if is_homomers:
            keep_indices = []
            for i in range(len(prev_df)):
                cmd = f"grep {prev_df.loc[i, 'pdbcode']} {self.complexes_cm_file}"
                contents = os.popen(cmd).read().split('\n')
                is_complex = len(contents[0].split()) > 1
                if is_complex:
                    keep_indices += [i]
            return prev_df.iloc[keep_indices]

        return prev_df

    def get_interactions_from_two_alignments(self, alignment1, alignment2, suffix, is_homomers):
        id_1 = []
        id_2 = []
        for id1 in alignment1.ids:
            for id2 in alignment2.ids:
                if id1 == id2 and not is_homomers:
                    id_1 += [id1]
                    id_2 += [id2]
                elif id1 in self.uniprot2pdb_map and id2 in self.uniprot2pdb_map:
                    if len(self.get_interactions_from_ids([id1, id2], is_homomers)) > 0:
                        id_1 += [id1]
                        id_2 += [id2]
        return pd.DataFrame({f"id_{suffix}": id_1, f"id_{suffix + 1}": id_2})

    def get_interactions(self, alignments, is_homomers=False):
        prev_df = None
        for i in range(0, len(alignments) - 1):
            curr_df = self.get_interactions_from_two_alignments(alignments[i], alignments[i + 1], i + 1, is_homomers)
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, on=f"id_{i + 1}")
        return prev_df

    def get_interactions_large(self, alignments, is_homomers=False):
        def create_df(alignment):
            ids = []
            pdbcodes = []
            for id in alignment.ids:
                if id in self.uniprot2pdb_map:
                    pdbcodes_for_id = [pdbcode[0:4] for pdbcode in self.uniprot2pdb_map[id]]
                    print(pdbcodes_for_id)
                    pdbcodes += pdbcodes_for_id
                    for i in range(len(pdbcodes_for_id)):
                        ids += [id]
            return pd.DataFrame({'id': ids, 'pdbcode': pdbcodes})

        prev_df = create_df(alignments[0])
        print(prev_df)
        for i in range(1, len(alignments)):
            curr_df = create_df(alignments[i])
            prev_df = prev_df.merge(curr_df, on="pdbcode", suffixes=(f"_{i}", f"_{i + 1}"))
        return prev_df
