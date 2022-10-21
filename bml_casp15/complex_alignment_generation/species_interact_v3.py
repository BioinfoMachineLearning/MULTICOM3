import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *

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

def read_species_annotation_table(annotation):

    SPECIES_ANNOTATION_COLUMNS = ["OS", "Tax"]

    data = annotation

    # initialize the column to extract the species information from
    annotation_column = None
    current_num_annotations = 0

    # Determine whether to extract based on the "OS" field
    # or the "Tax" field. Generally, OS is for Uniprot
    for column in SPECIES_ANNOTATION_COLUMNS:
        # if this column contains more non-null values
        if column not in data:
            continue

        num_annotations = sum(data[column].notnull())
        if num_annotations > current_num_annotations:
            # use that column to extract data
            annotation_column = column
            current_num_annotations = num_annotations

    # if we did not find an annotation column, return an error
    if annotation_column is None:
        return None

    # creates a new column called species with the species annotations
    data.loc[:, "species"] = data.loc[:, annotation_column]

    return data[["id", "name", "species"]]

def extract_header_annotation(alignment):
    columns = [
            ("GN", "gene"),
            ("OS", "organism"),
            ("PE", "existence_evidence"),
            ("SV", "sequence_version"),
            ("n", "num_cluster_members"),
            ("Tax", "taxon"),
            ("RepID", "representative_member")
        ]

    col_to_descr = OrderedDict(columns)
    regex = re.compile("\s({})=".format("|".join(col_to_descr.keys())))

    res = []
    for seq_idx, seq_id in enumerate(alignment.ids):
        full_header = alignment.headers[seq_idx]
        anno = None
        if ("GS" in alignment.annotation and
                    seq_id in alignment.annotation["GS"] and
                    "DE" in alignment.annotation["GS"][seq_id]):
            anno = alignment.annotation["GS"][seq_id]["DE"]
        else:
            split = full_header.split(maxsplit=1)
            if len(split) == 2:
                _, anno = split

        # extract info from line if we got one
        if anno is not None:
            # do split on known field names o keep things
            # simpler than a gigantic full regex to match
            # (some fields are allowed to be missing)
            pairs = re.split(regex, anno)
            pairs = ["id", seq_id, "name"] + pairs

            # create feature-value map
            feat_map = dict(zip(pairs[::2], pairs[1::2]))
            res.append(feat_map)
        else:
            res.append({"id": seq_id})

    df = pd.DataFrame(res)
    return df.reindex(
            ["id", "name"] + list(col_to_descr.keys()),
            axis=1
        )

def create_species_dict(msa_df):
    species_lookup = {}
    #print(msa_df)
    for species, species_df in msa_df.groupby('species'):
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
        
    annotation = extract_header_annotation(alignment)

    annotation_table = read_species_annotation_table(annotation)

    print(annotation_table)
        
    ids = []
    species = []
    msa_similarity = []
    msa_row = []

    if annotation_table is None:
        return pd.DataFrame({'id': ids, 'species': species, 'msa_similarity': msa_similarity, 'msa_row': msa_row})

    for index, row in annotation_table.iterrows():
        ids += [row.id]
        species += [row.species]
        per_seq_similarity = len([1 for j in range(len(alignment.main_seq)) if alignment.main_seq[j] == alignment.seqs[index][j]]) / float(len(alignment.main_seq))
        msa_similarity += [per_seq_similarity]
        msa_row += [index]

    print("make msa df222")
    return pd.DataFrame({'id': ids, 'species': species, 'msa_similarity': msa_similarity, 'msa_row': msa_row})

class Species_interact_v3:

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

