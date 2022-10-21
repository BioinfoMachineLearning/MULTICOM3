import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *
import os

HHBLITS_AA_TO_ID = {
    'A': 0,
    'B': 2,
    'C': 1,
    'D': 2,
    'E': 3,
    'F': 4,
    'G': 5,
    'H': 6,
    'I': 7,
    'J': 20,
    'K': 8,
    'L': 9,
    'M': 10,
    'N': 11,
    'O': 20,
    'P': 12,
    'Q': 13,
    'R': 14,
    'S': 15,
    'T': 16,
    'U': 1,
    'V': 17,
    'W': 18,
    'X': 20,
    'Y': 19,
    'Z': 3,
    '-': 21,
}

def create_species_dict(msa_df):
    species_lookup = {}
    #print(msa_df)
    for species, species_df in msa_df.groupby('pdbcode'):
        species_lookup[species] = species_df
        #print(species_df)
    return species_lookup

def match_rows_by_sequence_similarity(this_species_msa_dfs):
    all_paired_msa_rows = []
    num_seqs = [len(species_df) for species_df in this_species_msa_dfs if species_df is not None]
    take_num_seqs = np.min(num_seqs)
    sort_by_similarity = (lambda x: x.sort_values('msa_similarity', axis=0, ascending=False))
    for species_df in this_species_msa_dfs:
        if species_df is not None:
            species_df_sorted = sort_by_similarity(species_df)
            msa_rows = species_df_sorted.msa_row.iloc[:take_num_seqs].values
        else:
            msa_rows = [-1] * take_num_seqs  # take the last 'padding' row
        all_paired_msa_rows.append(msa_rows)
    all_paired_msa_rows = list(np.array(all_paired_msa_rows).transpose())
    return all_paired_msa_rows


def reorder_paired_rows(all_paired_msa_rows_dict):
    all_paired_msa_rows = []

    for num_pairings in sorted(all_paired_msa_rows_dict, reverse=True):
        paired_rows = all_paired_msa_rows_dict[num_pairings]
        print(paired_rows)
        paired_rows_product = abs(np.array([np.prod(rows) for rows in paired_rows]))
        print("test")
        print(paired_rows_product)
        paired_rows_sort_index = np.argsort(paired_rows_product)
        all_paired_msa_rows.extend(paired_rows[paired_rows_sort_index])

    return np.array(all_paired_msa_rows)

def parse_header(header):
    # discard anything in header that might come after the
    # first whitespace (by convention this is typically annotation)
    header = header.split()[0]

    # try to extract region from sequence header
    m = re.search("(.+)/(\d+)-(\d+)", header)
    if m:
        id_, start_str, end_str = m.groups()
        region_start, region_end = int(start_str), int(end_str)
        return id_, region_start, region_end
    else:
        # cannot find region, so just give back sequence iD
        return header, None, None

class PDB_interact_v3:

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

    def make_msa_df(self, alignment):
        print("make msa df111")
        ids = []
        pdbcodes = []
        msa_similarity = []
        msa_row = []
        for i, id in enumerate(alignment.ids):
            ids += [id]
            pdbcodes += [id]
            per_seq_similarity = len([1 for j in range(len(alignment.main_seq)) if alignment.main_seq[j] == alignment.seqs[i][j]]) / float(len(alignment.main_seq))
            msa_similarity += [per_seq_similarity]
            msa_row += [i]

            if id not in self.uniprot2pdb_map:
                continue

            unique_pdbcodes = set([pdbcode[0:4] for pdbcode in self.uniprot2pdb_map[id].split(',')])
            print(unique_pdbcodes)
            for pdbcode in unique_pdbcodes:
                ids += [id]
                pdbcodes += [pdbcode]
                per_seq_similarity = len([1 for j in range(len(alignment.main_seq)) if alignment.main_seq[j] == alignment.seqs[i][j]]) / float(len(alignment.main_seq))
                msa_similarity += [per_seq_similarity]
                msa_row += [i]

        print("make msa df222")
        return pd.DataFrame({'id': ids, 'pdbcode': pdbcodes, 'msa_similarity': msa_similarity, 'msa_row': msa_row})

    def get_interactions_v2(self, alignments, is_homomers=False):
        print("111111111111111111111111111")
        num_examples = len(alignments)
        print(len(alignments))

        all_chain_species_dict = []
        common_species = set()
        for chain_alignment in alignments:
            msa_df = self.make_msa_df(chain_alignment)
            species_dict = create_species_dict(msa_df)
            all_chain_species_dict.append(species_dict)
            common_species.update(set(species_dict))

        common_species = sorted(common_species)

        if is_homomers:
            common_species_filtered = set()
            for pdbcode in common_species:
                cmd = f"grep {pdbcode} {self.complexes_cm_file}"
                contents = os.popen(cmd).read().split('\n')
                is_complex = len(contents[0].split()) > 1
                if is_complex:
                    common_species_filtered.add(pdbcode)
            common_species = common_species_filtered

        all_paired_msa_rows = [np.zeros(num_examples, int)]
        all_paired_msa_rows_dict = {k: [] for k in range(num_examples)}
        all_paired_msa_rows_dict[num_examples] = [np.zeros(num_examples, int)]

        print(common_species)

        print("22222222222222222222222222222")
        for species in common_species:
            if not species:
                continue
            this_species_msa_dfs = []
            species_dfs_present = 0
            for species_dict in all_chain_species_dict:
                if species in species_dict:
                    this_species_msa_dfs.append(species_dict[species])
                    species_dfs_present += 1
                else:
                    this_species_msa_dfs.append(None)

            # Skip species that are present in only one chain.
            if species_dfs_present <= 1:
                continue

            paired_msa_rows = match_rows_by_sequence_similarity(this_species_msa_dfs)
            all_paired_msa_rows.extend(paired_msa_rows)
            all_paired_msa_rows_dict[species_dfs_present].extend(paired_msa_rows)

        all_paired_msa_rows_dict = {
            num_examples: np.array(paired_msa_rows) for
            num_examples, paired_msa_rows in all_paired_msa_rows_dict.items()
        }
        
        print(all_paired_msa_rows_dict)

        paired_rows = reorder_paired_rows(all_paired_msa_rows_dict)

        print(paired_rows)

        # make sure each id only apprears once in the concatenated alignment
        seen_ids = {k: [] for k in range(num_examples)}
        paired_rows_filtered = []
        for pair_index in range(len(paired_rows[:, 0])):
            if pair_index == 0:
                paired_rows_filtered += [paired_rows[pair_index, :]]
                continue
            print('test222222222222222222')
            row_indices = list(paired_rows[pair_index, :])
            pass_filter = True
            for j in range(len(alignments)):
                index = row_indices[j]
                if index == -1:
                    continue
                else:
                    header, start, end = parse_header(alignments[j].headers[index])
                    if f"{header}_{start}-{end}" in seen_ids[j]:
                        pass_filter = False
                        break
            
            if pass_filter:
                paired_rows_filtered += [paired_rows[pair_index]]
                for j in range(len(alignments)):
                    index = row_indices[j]
                    if index == -1:
                        continue
                    else:
                        header, start, end = parse_header(alignments[j].headers[index])
                        seen_ids[j] += [f"{header}_{start}-{end}"]

        return np.array(paired_rows_filtered)

