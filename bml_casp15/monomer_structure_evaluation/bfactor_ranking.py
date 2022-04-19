import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from bml_casp15.common.protein import parse_pdb_row


def get_avg_factor(infile):
    bfactors = []
    read = set()
    for line in open(infile):
        if line.startswith('ATOM'):
            resn = line[17:20]
            chain = line[21]
            resi = line[22:26]
            icode = line[26]
            bfactor = line[60:66]
            r_uid = (resn, chain, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
                #print(read)
            else:
                continue
            bfactors += [float(bfactor)]
    return np.mean(np.array(bfactors))


class Bfactor_qa:
    def __init__(self, params):
        self.params = params

    def run(self, input_dir):
        models = []
        avg_bfactors = []
        for pdb in os.listdir(input_dir):
            models += [pdb]
            avg_bfactors += [get_avg_factor(input_dir + '/' + pdb)]

        ranking_df = pd.DataFrame({'model': models, 'bfactor': avg_bfactors})
        ranking_df = ranking_df.sort_values(by=['bfactor'], ascending=False, ignore_index=True)
        return ranking_df
