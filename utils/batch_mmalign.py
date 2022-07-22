import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir


def run_command(inparams):
    mmalign_program, pdb1, pdb2 = inparams
    cmd = mmalign_program + ' ' + pdb1 + ' ' +  pdb2 + " | grep TM-score | awk '{print $2}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[1].rstrip('\n'))
    return pdb1, pdb2, tmscore


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True)
    parser.add_argument('--native_pdb', type=is_file, required=True)
    parser.add_argument('--mmalign_program', type=is_file, required=True)

    args = parser.parse_args()

    models = []
    scores = []
    process_list = []
    for model in os.listdir(args.indir):
        process_list.append([args.mmalign_program, args.native_pdb, args.indir + '/' + model])

    pool = Pool(processes=40)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()

    for result in results:
        pdb1, pdb2, tmscore = result
        models += [pdb2]
        scores += [tmscore]

    df = pd.DataFrame({'model': models, 'score': scores})

    print(df)

