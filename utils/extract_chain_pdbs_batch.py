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
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--chain_ids', type=str, required=True)
    args = parser.parse_args()

    chain_ids = args.chain_ids.split(',')

    makedir_if_not_exists(args.outdir)
    
    for inpdb in os.listdir(args.indir):
        chain_contents = {}
        for line in open(args.indir + '/' + inpdb):
            if line.startswith('ATOM'):
                chain_id = line[21]
                if chain_id not in chain_contents:
                    chain_contents[chain_id] = [line]
                else:
                    chain_contents[chain_id] += [line]

        with open(args.outdir + '/' + inpdb, 'w') as fw:
            for map in args.mapping.split(','):
                chain_src, chain_trg = map.split('_')
                for line in chain_contents[chain_src]:
                    fw.write(line[:21] + chain_trg + line[22:])
                fw.write("TER\n")
