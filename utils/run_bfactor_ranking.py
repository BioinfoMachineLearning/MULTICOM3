import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.quaternary_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa


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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inpdb', type=is_file, required=True)
    args = parser.parse_args()
    print(get_avg_factor(args.inpdb))

