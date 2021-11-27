import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *


class UNIPROT_distance:
    ### ColabFold
    def uni_num(ids):
        ########################################
        pa = {a: 0 for a in ascii_uppercase}
        for a in ["O", "P", "Q"]: pa[a] = 1
        ma = [[{} for k in range(6)], [{} for k in range(6)]]
        for n, t in enumerate(range(10)):
            for i in [0, 1]:
                for j in [0, 4]: ma[i][j][str(t)] = n

        for n, t in enumerate(list(ascii_uppercase) + list(range(10))):
            for i in [0, 1]:
                for j in [1, 2]: ma[i][j][str(t)] = n
            ma[1][3][str(t)] = n
        for n, t in enumerate(ascii_uppercase):
            ma[0][3][str(t)] = n
            for i in [0, 1]: ma[i][5][str(t)] = n
        ########################################
        nums = []
        for uni in ids:
            p = pa[uni[0]]
            tot, num = 1, 0
            if len(uni) == 10:
                for n, u in enumerate(reversed(uni[-4:])):
                    num += ma[p][n][u] * tot
                    tot *= len(ma[p][n].keys())
            for n, u in enumerate(reversed(uni[:6])):
                num += ma[p][n][u] * tot
                tot *= len(ma[p][n].keys())
            nums.append(num)
        return nums

    ### RoseTTAFold
    def uni2idx(ids):
        '''convert uniprot ids into integers according to the structure
        of uniprot accession numbers'''
        # ipdb.set_trace()
        ids2 = [i + 'AAA0' if len(i) == 6 else i for i in ids]
        arr = np.array([list(s) for s in ids2], dtype='|S1').view(np.uint8)

        for i in [1, 5, 9]:
            arr[:, i] -= ord('0')

        arr[arr >= ord('A')] -= ord('A')
        arr[arr >= ord('0')] -= ord('0') - 26

        arr[:, 0][arr[:, 0] > ord('Q') - ord('A')] -= 3

        arr = arr.astype(np.int64)

        coef = np.array([23, 10, 26, 36, 36, 10, 26, 36, 36, 1], dtype=np.int64)
        coef = np.tile(coef[None, :], [len(ids), 1])

        c1 = [i for i, id_ in enumerate(ids) if id_[0] in 'OPQ' and len(id_) == 6]
        c2 = [i for i, id_ in enumerate(ids) if id_[0] not in 'OPQ' and len(id_) == 6]

        coef[c1] = np.array([3, 10, 36, 36, 36, 1, 1, 1, 1, 1])
        coef[c2] = np.array([23, 10, 26, 36, 36, 1, 1, 1, 1, 1])

        for i in range(1, 10):
            coef[:, -i - 1] *= coef[:, -i]

        return np.sum(arr * coef, axis=-1)

    def filter_ids(ids):
        filtered_ids = []
        regex = re.compile("^[A-Za-z0-9]+$")
        for id in ids:
            add = True
            if re.fullmatch(regex, id.split('_')[-1]) is None:
                add = False
            if id.split('_')[-1].startswith("UPI"):
                add = False
            if add:
                filtered_ids += [id]
        return filtered_ids

    def get_interactions(alignment1, alignment2):

        ids1 = UNIPROT_distance.filter_ids(alignment1.ids)
        ids2 = UNIPROT_distance.filter_ids(alignment2.ids)

        if len(ids1) == 0 or len(ids2) == 0:
            return pd.DataFrame({"id_1": [], "id_2": []})

        hash1 = UNIPROT_distance.uni2idx(ids1)
        hash2 = UNIPROT_distance.uni2idx(ids2)

        df_dict = {"id_1": [], "id_2": []}
        for i in range(len(hash1)):
            id1 = hash1[i]
            for j in range(len(hash2)):
                id2 = hash2[j]
                if abs(id1-id2) < 10:
                    df_dict['id_1'] += [ids1[i]]
                    df_dict['id_2'] += [ids2[j]]
                    
        return pd.DataFrame(df_dict)
