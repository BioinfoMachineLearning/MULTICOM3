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


def run_cmd(cmd):
    print(cmd)
    try:
        os.system(cmd)
    except Exception as e:
        print(e)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fastadir', type=is_dir, required=True)
    parser.add_argument('--monomer_model_dir', type=is_dir, required=True)
    parser.add_argument('--modeldir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=is_dir, required=True)
    args = parser.parse_args()

    eva_script = '/home/bml_casp15/BML_CASP15/bin/N8_evaluate_multimer_structures.py'
    option_file = '/home/bml_casp15/BML_CASP15/bin/db_option_test'

    process_list = []
    for dimer in os.listdir(args.modeldir):
        cmd = f"python {eva_script} --option_file {option_file} --fasta_path {args.fastadir}/{dimer}.fasta " \
              f"--input_dir {args.modeldir}/{dimer} " \
              f"--output_dir {args.output_dir}/{dimer} --monomer_model_dir {args.monomer_model_dir} " \
              f"--stoichiometry A1B1 > /dev/null 2>/dev/null "
        process_list += [cmd]

    pool = Pool(processes=5)
    results = pool.map(run_cmd, process_list)
    pool.close()
    pool.join()
