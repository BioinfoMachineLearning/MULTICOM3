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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputdir', type=is_file, required=True)
    parser.add_argument('--atom_dir', type=is_dir, required=True)
    args = parser.parse_args()
    tmscore_program = '/home/bml_casp15/BML_CASP15/tools/TMscore'

    for target in os.listdir(args.inputdir):

        summary_csv = args.inputdir + '/' + target + '/summary.csv'

        summary_df = pd.read_csv(summary_csv).T



    for method in os.listdir(args.output_dir):
        corr_summary_all = pd.DataFrame(columns=['dockq', 'tmscore1', 'tmscore2'])
        corr_summary_dimer = pd.DataFrame(columns=['dockq', 'tmscore1', 'tmscore2'])

        process_list = []
        for line in open(args.dimerlist).readlines():
            chain1, chain2 = line.rstrip('\n').split()
            process_list.append([f"{chain1}_{chain2}", args.atom_dir, f"{args.output_dir}/{method}/{chain1}_{chain2}",
                                 args.tmpdir, dockq_program, tmscore_program])
        pool = Pool(processes=20)
        results = pool.map(evaluate_dimer, process_list)
        pool.close()
        pool.join()

        for result in results:
            corr_summary_all = corr_summary_all.append(result['corr_summary_all_pd'])
            corr_summary_dimer = corr_summary_dimer.append(result['corr_summary_dimer_pd'])

        # corr_summary_all.to_csv(f'{method}_all.csv')
        corr_summary_dimer.to_csv(f'{method}_avg.csv')
