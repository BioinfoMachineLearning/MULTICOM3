import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *


class STRING_interact_v2:

    def __init__(self, string2uniprot_map_file):

        self.mapping_file = string2uniprot_map_file

    def load_data(self, score_threshold=0.5):

        self.interaction_map = {}

        with open(self.mapping_file) as f:
            for line in f:
                uniprotid1, uniprotid2, score = line.rstrip('\n').split()
                if int(score) < score_threshold:
                    continue
                uniprotid1 = uniprotid1.split('_')[1]
                uniprotid2 = uniprotid2.split('_')[1]

                if uniprotid1 in self.interaction_map:
                    self.interaction_map[uniprotid1] += "," + uniprotid2
                else:
                    self.interaction_map[uniprotid1] = uniprotid2

                if uniprotid2 in self.interaction_map:
                    self.interaction_map[uniprotid2] += "," + uniprotid1
                else:
                    self.interaction_map[uniprotid2] = uniprotid1

    def get_interactions_from_two_alignments(self, alignment1, alignment2, suffix):
        id_1 = []
        id_2 = []
        for id1 in alignment1.ids:
            interact_ids = []
            if id1 in self.interaction_map:
                interact_ids = self.interaction_map[id1].split(',')
            for id2 in alignment2.ids:
                if id2 == id1 or id2 in interact_ids:
                    id_1 += [id1]
                    id_2 += [id2]
        return pd.DataFrame({f"id_{suffix}": id_1, f"id_{suffix + 1}": id_2})

    def get_interactions(self, alignments):
        prev_df = None
        for i in range(0, len(alignments) - 1):
            curr_df = self.get_interactions_from_two_alignments(alignments[i], alignments[i+1], suffix=i+1)
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, on=f"id_{i + 1}")
        return prev_df
