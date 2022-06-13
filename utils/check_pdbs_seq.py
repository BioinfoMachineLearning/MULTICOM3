import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_file, required=True)
    parser.add_argument('--seq', type=str, required=True)
    args = parser.parse_args()

    pdb2seq = '/home/bml_casp15/BML_CASP15/utils/pdb2seq.pl'

    for pdb in os.listdir(args.indir):
        cmd = f"perl {pdb2seq} {args.indir}/{pdb}"
        print(cmd)
        seq = os.popen(cmd).read().strip()
        print(seq)
        if seq != args.seq:
            raise Exception(f"sequence for {pdb} is wrong")
