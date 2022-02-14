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


def split_pdb(complex_pdb, outpdbs):
    pre_chain = None
    i = 0
    fw = open(outpdbs[i], 'w')

    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if pre_chain is None:
            pre_chain = chain_name
            fw.write(line)
        elif chain_name == pre_chain:
            fw.write(line)
        else:
            fw.close()
            i = i + 1
            fw = open(outpdbs[i], 'w')
            fw.write(line)
            pre_chain = chain_name
    fw.close()


def complete_result(outputdir):
    complete = True
    for i in range(1, 6):
        model = f'{outputdir}/unrelaxed_model_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def get_range_top_of_contact_map(input_map, range='long', top='l2'):
    length = input_map.shape[0]
    map_copy = np.copy(input_map)
    if range == 'long':
        map_copy = np.triu(map_copy, 23)
    elif range == 'medium':
        map_copy = np.triu(map_copy, 11) - np.triu(map_copy, 23)

    map_copy = map_copy.flatten()
    if top == 'l2':
        map_copy[map_copy.argsort()[::-1][0:int(length / 2)]] = 1
        map_copy[map_copy < 1] = 0
    elif top == 'l5':
        map_copy[map_copy.argsort()[::-1][0:int(length / 5)]] = 1
        map_copy[map_copy < 1] = 0
    elif top == 'l':
        map_copy[map_copy.argsort()[::-1][0:int(length)]] = 1
        map_copy[map_copy < 1] = 0
    map_copy = map_copy.reshape(length, length)
    return map_copy


def evaluate_monomer(inparams):

    monomer, atomdir, outputdir, tmscore_program = inparams

    native_atom = f"{atomdir}/{monomer}.atom"

    if not os.path.exists(native_atom):
        print(f"Cannot find {native_atom} for {monomer}")
        return None, None

    avg_tmscore = 0
    if complete_result(outputdir):
        for i in range(1, 6):
            model = f'{outputdir}/unrelaxed_model_{i}.pdb'
            cmd = tmscore_program + f' {model} {native_atom}' + " | grep TM-score | awk '{print $3}' "
            #print(cmd)
            tmscore_contents = os.popen(cmd).read().split('\n')
            tmscore = float(tmscore_contents[2].rstrip('\n'))
            avg_tmscore += tmscore
        avg_tmscore = avg_tmscore / 5

    return monomer, avg_tmscore


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--atom_dir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=is_dir, required=True)
    parser.add_argument('--monomer_list', type=is_file, required=False)

    args = parser.parse_args()

    tmscore_program = '/home/bml_casp15/BML_CASP15/tools/TMscore'

    monomers_list = os.listdir(args.output_dir)
    if args.monomer_list is not None:
        monomers_list = [line.rstrip() for line in open(args.monomer_list).readlines()]
    print(monomers_list)
    process_list = []
    for monomer in os.listdir(args.output_dir):
        if monomer in monomers_list:
            process_list.append([monomer, args.atom_dir, args.output_dir + '/' + monomer, tmscore_program])

    pool = Pool(processes=15)
    results = pool.map(evaluate_monomer, process_list)
    pool.close()
    pool.join()

    monomers = []
    tmscores = []
    for result in results:
        monomer, tmscore = result
        if monomer is not None and tmscore is not None:
            monomers += [monomer]
            tmscores += [tmscore]

    df = pd.DataFrame({'monomer': monomers, 'tmscore': tmscores})
    df.to_csv(f'monomer_tmscore.csv')
