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
    parser.add_argument('--inpdb', type=str, required=True)
    parser.add_argument('--chain', type=str, required=True)
    parser.add_argument('--outpdb', type=str, required=True)
    args = parser.parse_args()

    chain_contents = {}
    for line in open(args.inpdb):
        if line.startswith('ATOM'):
            chain_id = line[21]
            if chain_id not in chain_contents:
                chain_contents[chain_id] = [line]
            else:
                chain_contents[chain_id] += [line]

    with open(args.outpdb, 'a') as fw:
        for line in chain_contents[args.chain]:
            if line.startswith('ATOM'):
                fw.write(line)
