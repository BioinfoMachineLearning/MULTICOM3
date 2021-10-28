import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.alignment import alignment


class UNICLUST_oxmatch:

    def get_interactions(alignment1, alignment2):

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
            {"id_1": [alignment1.ids[i] for i in idx1], "id_2": [alignment2.ids[j] for j in idx2]}
        )
