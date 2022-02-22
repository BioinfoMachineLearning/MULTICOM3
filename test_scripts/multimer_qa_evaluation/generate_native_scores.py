import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
import torch
import torch.nn as nn
from bml_casp15.quaternary_structure_evaluation.pairwise_dockq import *
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir


def combine_pdb(pdb1, pdb2, combine_pdb):
    with open(combine_pdb, 'w') as out:
        for line in open(pdb1):
            if not line.startswith('ATOM'): continue
            out.write(line[:21] + 'B' + line[22:])
        for line in open(pdb2):
            if not line.startswith('ATOM'): continue
            out.write(line[:21] + 'C' + line[22:])

def cal_scores(dockq_program, model, combine_atom):
    cmd = "python " + dockq_program + " -short " + model + " " + combine_atom + "| grep DockQ | awk " \
                                                                                "'{print $2}' "
    # print(cmd)
    score = os.popen(cmd).read()
    if len(score.rstrip('\n')) == 0:
        return -1
    return float(score)

def generate_native_score(inparams):

    atomdir, dimer, inputdir, dockq_program = inparams

    chain1, chain2 = dimer.split('_')

    chain1_atom = f"{atomdir}/{chain1}.atom"

    chain2_atom = f"{atomdir}/{chain2}.atom"

    combine_atom = f"{inputdir}/{dimer}/{chain1}_{chain2}.atom"

    combine_pdb(chain1_atom, chain2_atom, combine_atom)

    result_csv = f"{inputdir}/{dimer}/native_scores.csv"
    if os.path.exists(result_csv):
        print(f"The native scores for {dimer} has been generated: {result_csv}")
        return

    models = []
    dockq_scores = []
    for pdb in os.listdir(inputdir + '/' + dimer + '/pdb'):
        dockq_score = cal_scores(dockq_program, f"{inputdir}/{dimer}/pdb/{pdb}", combine_atom)
        models += [pdb]
        dockq_scores += [dockq_score]

    df = pd.DataFrame({'model': models, 'dockq_score': dockq_scores})
    df = df.sort_values(by=['dockq_score'], ascending=False)
    df.to_csv(result_csv)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputdir', type=is_dir, required=True)
    parser.add_argument('--atomdir', type=is_dir, required=True)
    args = parser.parse_args()

    dockq_program = '/home/bml_casp15/BML_CASP15/tools/DockQ/DockQ.py'

    process_list = []
    for dimer in os.listdir(args.inputdir):
        process_list.append([args.atomdir, dimer, args.inputdir, dockq_program])

    pool = Pool(processes=10)
    results = pool.map(generate_native_score, process_list)
    pool.close()
    pool.join()