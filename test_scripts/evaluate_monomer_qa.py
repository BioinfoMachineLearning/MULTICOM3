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
    parser.add_argument('--dimerlist', type=is_file, required=True)
    parser.add_argument('--atom_dir', type=is_dir, required=True)
    parser.add_argument('--fastadir', type=is_dir, required=True)
    parser.add_argument('--monomerdir', type=is_dir, required=True)
    parser.add_argument('--dimerdir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=is_dir, required=True)
    args = parser.parse_args()

    eva_script = '/home/bml_casp15/BML_CASP15/bin/N7_evaluate_monomer_structures.py'
    option_file = '/home/bml_casp15/BML_CASP15/bin/db_option'

    for dimer in open(args.dimerlist):
        chain1, chain2 = dimer.rstrip('\n').split()
        chain_in_multimer_dict = {chain1: 'B', chain2: 'C'}
        input_multimer_dir = f"{args.dimerdir}/{chain1}_{chain2}"
        for chain in [chain1, chain2]:
            fastafile = f"{args.fastadir}/{chain}.fasta"
            input_monomer_dir = f"{args.monomerdir}/{chain}"
            chain_in_multimer = chain_in_multimer_dict[chain]
            cmd = f"python {eva_script} --option_file {option_file} --targetname {chain} --fasta_file {fastafile} " \
                  f"--input_monomer_dir {input_monomer_dir} --input_multimer_dir {input_multimer_dir} " \
                  f"--chain_in_multimer {chain_in_multimer} --output_dir {args.output_dir} --use_gpu=false > " \
                  f"/dev/null 2>/dev/null "
            print(cmd)
            try:
                os.system(cmd)
            except Exception as e:
                print(e)