import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--cmd', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--times', type=int, default=5)

    args = parser.parse_args()

    makedir_if_not_exists(args.outdir)

    with open(args.outdir + '/run_command.sh', 'w') as fw:
        for i in range(args.times):
            fw.write(args.cmd + '\n')

    print(args.outdir + '/run_command.sh')

