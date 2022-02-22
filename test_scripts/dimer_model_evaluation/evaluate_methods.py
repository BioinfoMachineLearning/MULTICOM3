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


def evaluate_dimer(inparams):
    dimer, atomdir, outputdir, tmpdir, dockq_program, tmscore_program, mmalign_program = inparams

    chain1, chain2 = dimer.split('_')

    chain1_atom = f"{atomdir}/{chain1}.atom"

    chain2_atom = f"{atomdir}/{chain2}.atom"

    combine_atom = f"{outputdir}/{chain1}_{chain2}.atom"

    combine_pdb(chain1_atom, chain2_atom, combine_atom)

    tmpdir = tmpdir + '/' + dimer

    model_count = 0
    dimer_count = 0
    corr_summary_all = pd.DataFrame(columns=['dockq', 'tmscore1', 'tmscore2', 'mmalign'])
    corr_summary_dimer = pd.DataFrame(columns=['dockq', 'tmscore1', 'tmscore2', 'mmalign'])


    dockQ_scores = []

    dimer_result = {'dimer': dimer,
                    'dockq': 0,
                    'tmscore1': 0,
                    'tmscore2': 0,
                    'mmalign': 0}

    if complete_result(f'{outputdir}'):
        for i in range(1, 6):
            model = f'{outputdir}/unrelaxed_model_{i}_multimer.pdb'
            cmd = "python " + dockq_program + " -short " + model + " " + combine_atom + "| grep DockQ | awk " \
                                                                                        "'{print $2}' "
            # print(cmd)
            score = os.popen(cmd).read()
            if len(score.rstrip('\n')) == 0:
                bfail = True
                break
            dockQ_scores += [float(score.rstrip('\n'))]

            splitdir = f'{outputdir}/u_model_{i}_s'
            makedir_if_not_exists(splitdir)
            split_pdb(model, [f'{outputdir}/u_model_{i}_s/B.pdb',
                              f'{outputdir}/u_model_{i}_s/C.pdb'])

            cmd = tmscore_program + f' {outputdir}/u_model_{i}_s/B.pdb ' + chain1_atom + " | grep TM-score | awk '{print $3}' "
            # print(cmd)
            tmscore_contents = os.popen(cmd).read().split('\n')
            tmscore1 = float(tmscore_contents[2].rstrip('\n'))

            cmd = tmscore_program + f' {outputdir}/u_model_{i}_s/C.pdb ' + chain2_atom + " | grep TM-score | awk '{print $3}' "
            # print(cmd)
            tmscore_contents = os.popen(cmd).read().split('\n')
            tmscore2 = float(tmscore_contents[2].rstrip('\n'))

            cmd = mmalign_program + f' {model} ' + combine_atom + " | grep TM-score | awk '{print $2}' "
            # print(cmd)
            mmalign_contents = os.popen(cmd).read().split('\n')
            mmalign_score = float(mmalign_contents[1].rstrip('\n'))

            sum_dict = {'dimer': dimer,
                        'model': f"{method}_{i}",
                        'dockq': float(score.rstrip('\n')),
                        'tmscore1': tmscore1,
                        'tmscore2': tmscore2,
                        'mmalign': mmalign_score}

            corr_summary_all = corr_summary_all.append(pd.DataFrame(sum_dict, index=[model_count]))

            model_count += 1

            dimer_result['dockq'] += float(score.rstrip('\n'))
            dimer_result['tmscore1'] += tmscore1
            dimer_result['tmscore2'] += tmscore2
            dimer_result['mmalign'] += mmalign_score

    dimer_result['dockq'] = dimer_result['dockq'] / 5
    dimer_result['tmscore1'] = dimer_result['tmscore1'] / 5
    dimer_result['tmscore2'] = dimer_result['tmscore2'] / 5
    dimer_result['mmalign'] = dimer_result['mmalign'] / 5
    corr_summary_dimer = corr_summary_dimer.append(pd.DataFrame(dimer_result, index=[dimer_count]))
    dimer_count += dimer_count

    return {'corr_summary_all_pd': corr_summary_all,
            'corr_summary_dimer_pd': corr_summary_dimer}


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dimerlist', type=is_file, required=True)
    parser.add_argument('--atom_dir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=is_dir, required=True)
    parser.add_argument('--tmpdir', type=is_dir, required=True)

    args = parser.parse_args()

    dockq_program = '/home/bml_casp15/BML_CASP15/tools/DockQ/DockQ.py'
    tmscore_program = '/home/bml_casp15/BML_CASP15/tools/TMscore'
    mmalign_program = '/home/bml_casp15/BML_CASP15/tools/MMalign'

    for method in os.listdir(args.output_dir):
        corr_summary_all = pd.DataFrame(columns=['dockq', 'tmscore1', 'tmscore2', 'mmalign'])
        corr_summary_dimer = pd.DataFrame(columns=['dockq', 'tmscore1', 'tmscore2', 'mmalign'])

        process_list = []
        for line in open(args.dimerlist).readlines():
            chain1, chain2 = line.rstrip('\n').split()
            process_list.append([f"{chain1}_{chain2}", args.atom_dir, f"{args.output_dir}/{method}/{chain1}_{chain2}",
                                 args.tmpdir, dockq_program, tmscore_program, mmalign_program])
        pool = Pool(processes=20)
        results = pool.map(evaluate_dimer, process_list)
        pool.close()
        pool.join()

        for result in results:
            corr_summary_all = corr_summary_all.append(result['corr_summary_all_pd'])
            corr_summary_dimer = corr_summary_dimer.append(result['corr_summary_dimer_pd'])

        # corr_summary_all.to_csv(f'{method}_all.csv')
        corr_summary_dimer.to_csv(f'{method}_avg.csv')
