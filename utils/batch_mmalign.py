import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True)
    parser.add_argument('--native_pdb', type=is_file, required=True)
    parser.add_argument('--mmalign_program', type=is_file, required=True)

    args = parser.parse_args()

    models = []
    scores = []

    for model in os.listdir(args.indir):
        ref_tmscore = cal_tmscore(args.mmalign_program, args.indir + '/' + model, args.native_pdb)
        models += [model]
        scores += [ref_tmscore]

    df = pd.DataFrame({'model': models, 'score': scores})

    print(df)

