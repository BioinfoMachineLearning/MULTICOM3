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


def cal_scores(tmscore_program, lddt_program, inpdb, nativepdb):
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $3}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[2].rstrip('\n'))
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep GDT-score | awk '{print $3}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    gdtscore = float(tmscore_contents[0].rstrip('\n'))

    # cmd = lddt_program + ' ' + inpdb + ' ' + nativepdb + " | grep Global | awk '{print $4}' "
    # print(cmd)
    # lddt_contents = os.popen(cmd).read().split('\n')
    # lddtscore = float(lddt_contents[0].rstrip('\n'))
    lddtscore=0
    return tmscore, gdtscore, lddtscore

def cal_batch(inparams):
    inputdir, target, native_pdb, tmscore_program, lddt_program, result_csv = inparams
    models = []
    gdt_scores = []
    tmscores = []
    lddtscores = []
    for pdb in os.listdir(args.inputdir + '/' + target + '/pdb'):
        tmscore, gdtscore, lddtscore = cal_scores(tmscore_program, lddt_program, f"{args.inputdir}/{target}/pdb/{pdb}",
                                                  native_pdb)
        models += [pdb]
        gdt_scores += [gdtscore]
        tmscores += [tmscore]
        lddtscores += [lddtscore]

    df = pd.DataFrame({'model': models, 'gdtscore': gdt_scores, 'tmscore': tmscores, 'global_lddt': lddtscores})
    df.to_csv(result_csv)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputdir', type=is_dir, required=True)
    parser.add_argument('--atomdir', type=is_dir, required=True)
    args = parser.parse_args()
    tmscore_program = '/home/bml_casp15/BML_CASP15/tools/TMscore'
    lddt_program = '/home/bml_casp15/BML_CASP15/tools/lddt'

    process_list = []

    for target in os.listdir(args.inputdir):

        native_pdb = f"{args.atomdir}/{target}.atom"
        if not os.path.exists(native_pdb):
            raise Exception(f"Cannot find the native pdb for {target}: {native_pdb}")

        result_csv = f"{args.inputdir}/{target}/native_scores.csv"
        # if os.path.exists(result_csv):
        #     df = pd.read_csv(result_csv)
        #     if 'global_lddt' in df:
        #         print(f"The native scores for {target} has been generated: {result_csv}")
        #         continue
        process_list.append([args.inputdir, target, native_pdb, tmscore_program, lddt_program, result_csv])

    pool = Pool(processes=40)
    results = pool.map(cal_batch, process_list)
    pool.close()
    pool.join()
