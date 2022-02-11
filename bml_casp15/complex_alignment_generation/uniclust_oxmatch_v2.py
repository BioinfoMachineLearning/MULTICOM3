import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *


class UNICLUST_oxmatch_v2:

    def get_interactions_from_two_alignment(alignment1, alignment2, suffix):

        def match_top_ox(ox1, ox2):
            matching_ox = np.intersect1d(ox1, ox2)

            ind1 = []  # Index to select from the individual MSAs
            ind2 = []
            # Go through all matching and select the first (top) hit
            for ox in matching_ox:
                if ox == -1:
                    continue
                ind1.append(min(np.argwhere(ox1 == ox)[:, 0]))
                ind2.append(min(np.argwhere(ox2 == ox)[:, 0]))

            return ind1, ind2

        def get_ox(alignment):
            species = []
            for id in alignment.ids:
                ox = -1
                gap_fraction = alignment[id].count('-') / float(len(alignment[id]))
                if gap_fraction <= 0.9:  # Only use the lines with less than 90 % gaps
                    header = alignment.headers[id][0]
                    if 'OX=' in header:
                        OX = header.split('OX=')[1]
                        if len(OX) > 0:
                            ox = int(OX.split(' ')[0])
                species.append(ox)
            return species

        ox1 = get_ox(alignment1)
        ox2 = get_ox(alignment2)

        idx1, idx2 = match_top_ox(ox1, ox2)

        return pd.DataFrame(
            {f"id_{suffix}": [alignment1.ids[i] for i in idx1], f"id_{suffix+1}": [alignment2.ids[j] for j in idx2]}
        )

    def get_interactions(alignments):
        prev_df = None
        for i in range(0, len(alignments) - 1):
            curr_df = UNICLUST_oxmatch_v2.get_interactions_from_two_alignment(alignments[i], alignments[i+1], suffix=i+1)
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, on=f"id_{i + 1}")
        return prev_df
