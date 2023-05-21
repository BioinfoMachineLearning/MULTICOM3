import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from multicom3.monomer_alignment_generation.alignment import *

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
    for species, species_df in msa_df.groupby('ox'):
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

def make_msa_df(alignment):
    print("make msa df111")
    ids = []
    ox_species = []
    msa_similarity = []
    msa_row = []
    for i, id in enumerate(alignment.ids):
        ox = -1
        gap_fraction = alignment[id].count('-') / float(len(alignment[id]))
        if gap_fraction <= 0.9:  # Only use the lines with less than 90 % gaps
            header = alignment.headers[i]
            if 'OX=' in header:
                OX = header.split('OX=')[1]
                if len(OX) > 0:
                    ox = int(OX.split(' ')[0])
            
        if ox != -1:  
            ids += [id]
            ox_species += [ox]
            per_seq_similarity = len([1 for j in range(len(alignment.main_seq)) if alignment.main_seq[j] == alignment.seqs[i][j]]) / float(len(alignment.main_seq))
            msa_similarity += [per_seq_similarity]
            msa_row += [i]

    print("make msa df222")
    return pd.DataFrame({'id': ids, 'ox': ox_species, 'msa_similarity': msa_similarity, 'msa_row': msa_row})


class UNICLUST_oxmatch_v3:

    def get_interactions_v2(alignments):
        print("111111111111111111111111111")
        num_examples = len(alignments)
        print(len(alignments))

        all_chain_species_dict = []
        common_species = set()
        for chain_alignment in alignments:
            msa_df = make_msa_df(chain_alignment)
            species_dict = create_species_dict(msa_df)
            all_chain_species_dict.append(species_dict)
            common_species.update(set(species_dict))

        common_species = sorted(common_species)

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

        return paired_rows