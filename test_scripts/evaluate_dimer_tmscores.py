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
        model = f'{outputdir}/unrelaxed_model_{i}_multimer.pdb'
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


def evaluate_dimer(inparams):

    dimer, atomdir, outputdir, tmpdir, dockq_program, tmscore_program, monomer_list = inparams

    chain1, chain2 = dimer.split('_')

    chain1_atom = f"{atomdir}/{chain1}.atom"

    chain2_atom = f"{atomdir}/{chain2}.atom"

    combine_atom = f"{outputdir}/{dimer}/{chain1}_{chain2}.atom"

    combine_pdb(chain1_atom, chain2_atom, combine_atom)

    tmpdir = tmpdir + '/' + dimer

    monomer_count = 0
    corr_summary_monomers = pd.DataFrame(columns=['monomer', 'tmscore'])

    for method in ["alphafold"]:
        if complete_result(f'{outputdir}/{dimer}/{method}'):
            monomer_result_1 = {'monomer': chain1, 'tmscore': 0}
            monomer_result_2 = {'monomer': chain2, 'tmscore': 0}

            for i in range(1, 6):
                model = f'{outputdir}/{dimer}/{method}/unrelaxed_model_{i}_multimer.pdb'

                splitdir = f'{outputdir}/{dimer}/{method}/u_model_{i}_s'
                makedir_if_not_exists(splitdir)
                split_pdb(model, [f'{outputdir}/{dimer}/{method}/u_model_{i}_s/B.pdb',
                                  f'{outputdir}/{dimer}/{method}/u_model_{i}_s/C.pdb'])

                cmd = tmscore_program + f' {outputdir}/{dimer}/{method}/u_model_{i}_s/B.pdb ' + chain1_atom + " | grep TM-score | awk '{print $3}' "
                # print(cmd)
                tmscore_contents = os.popen(cmd).read().split('\n')
                tmscore1 = float(tmscore_contents[2].rstrip('\n'))
                monomer_result_1['tmscore'] += tmscore1

                cmd = tmscore_program + f' {outputdir}/{dimer}/{method}/u_model_{i}_s/C.pdb ' + chain2_atom + " | grep TM-score | awk '{print $3}' "
                # print(cmd)
                tmscore_contents = os.popen(cmd).read().split('\n')
                tmscore2 = float(tmscore_contents[2].rstrip('\n'))
                monomer_result_2['tmscore'] += tmscore2

            monomer_result_1['tmscore'] = monomer_result_1['tmscore'] / 5
            monomer_result_2['tmscore'] = monomer_result_2['tmscore'] / 5

            if monomer_list is not None:
                if chain1 in monomer_list:
                    corr_summary_monomers = corr_summary_monomers.append(
                        pd.DataFrame(monomer_result_1, index=[monomer_count]))
                    monomer_count += 1
                if chain2 in monomer_list:
                    corr_summary_monomers = corr_summary_monomers.append(
                        pd.DataFrame(monomer_result_2, index=[monomer_count]))
                    monomer_count += 1
            else:
                corr_summary_monomers = corr_summary_monomers.append(pd.DataFrame(monomer_result_1, index=[monomer_count]))
                monomer_count += 1
                corr_summary_monomers = corr_summary_monomers.append(pd.DataFrame(monomer_result_2, index=[monomer_count]))
                monomer_count += 1

    return corr_summary_monomers


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dimerlist', type=is_file, required=True)
    parser.add_argument('--atom_dir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=is_dir, required=True)
    parser.add_argument('--tmpdir', type=is_dir, required=True)
    parser.add_argument('--monomer_list', type=is_file, required=False)

    args = parser.parse_args()

    dockq_program = '/home/bml_casp15/BML_CASP15/tools/DockQ/DockQ.py'
    tmscore_program = '/home/bml_casp15/BML_CASP15/tools/TMscore'

    pairwise_qa = Pairwise_dockq_qa(dockq_program)

    methods = ['no_templates',
               'original_template_pipeline',
               'custom_template_file',
               'custom_template_with_alphafold_models']

    monomer_list = []
    if args.monomer_list is not None:
        monomers_list = [line.rstrip() for line in open(args.monomer_list).readlines()]
    print(monomers_list)

    for method in methods:
        corr_summary_all = pd.DataFrame(columns=['monomer', 'tmscore'])
        process_list = []
        for line in open(args.dimerlist).readlines():
            chain1, chain2 = line.rstrip('\n').split()
            process_list.append([f"{chain1}_{chain2}", args.atom_dir, args.output_dir + '/' + method,
                                 args.tmpdir, dockq_program, tmscore_program, monomers_list])
        pool = Pool(processes=30)
        results = pool.map(evaluate_dimer, process_list)
        pool.close()
        pool.join()

        for result in results:
            corr_summary_all = corr_summary_all.append(result)

        corr_summary_all.to_csv(f'{method}_tmscore.csv')
