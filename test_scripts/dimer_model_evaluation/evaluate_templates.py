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
    dimer, atomdir, outputdir, tmpdir, dockq_program, tmscore_program, mmalign_program = inparams

    chain1, chain2 = dimer.split('_')

    chain1_atom = f"{atomdir}/{chain1}.atom"

    chain2_atom = f"{atomdir}/{chain2}.atom"

    combine_atom = f"{outputdir}/{dimer}/{chain1}_{chain2}.atom"

    combine_pdb(chain1_atom, chain2_atom, combine_atom)

    tmpdir = tmpdir + '/' + dimer

    model_count = 0
    dimer_count = 0
    corr_summary_all = pd.DataFrame(
        columns=['dockq', 'tmscore1', 'tmscore2', 'mmalign', 'lddt', 'avg_l2_prob', 'avg_l5_prob', 'ptmscore'])
    corr_summary_dimer = pd.DataFrame(
        columns=['dockq', 'tmscore1', 'tmscore2', 'mmalign', 'lddt', 'avg_l2_prob', 'avg_l5_prob', 'ptmscore'])

    for method in ["alphafold"]:
        dockQ_scores = []
        if complete_result(f'{outputdir}/{dimer}/{method}'):
            dimer_result = {'dimer': dimer,
                            'dockq': 0,
                            'tmscore1': 0,
                            'tmscore2': 0,
                            'mmalign': 0,
                            'lddt': 0,
                            'avg_l2_prob': 0,
                            'avg_l5_prob': 0,
                            'ptmscore': 0}

            for i in range(1, 6):
                model = f'{outputdir}/{dimer}/{method}/unrelaxed_model_{i}_multimer.pdb'
                cmd = "python " + dockq_program + " -short " + model + " " + combine_atom + "| grep DockQ | awk " \
                                                                                            "'{print $2}' "
                # print(cmd)
                score = os.popen(cmd).read()
                if len(score.rstrip('\n')) == 0:
                    bfail = True
                    break
                dockQ_scores += [float(score.rstrip('\n'))]

                splitdir = f'{outputdir}/{dimer}/{method}/u_model_{i}_s'
                makedir_if_not_exists(splitdir)
                split_pdb(model, [f'{outputdir}/{dimer}/{method}/u_model_{i}_s/B.pdb',
                                  f'{outputdir}/{dimer}/{method}/u_model_{i}_s/C.pdb'])

                cmd = tmscore_program + f' {outputdir}/{dimer}/{method}/u_model_{i}_s/B.pdb ' + chain1_atom + " | grep TM-score | awk '{print $3}' "
                # print(cmd)
                tmscore_contents = os.popen(cmd).read().split('\n')
                tmscore1 = float(tmscore_contents[2].rstrip('\n'))

                cmd = tmscore_program + f' {outputdir}/{dimer}/{method}/u_model_{i}_s/C.pdb ' + chain2_atom + " | grep TM-score | awk '{print $3}' "
                # print(cmd)
                tmscore_contents = os.popen(cmd).read().split('\n')
                tmscore2 = float(tmscore_contents[2].rstrip('\n'))

                cmd = mmalign_program + f' {outputdir}/{dimer}/{method}/u_model_{i}_s/C.pdb ' + combine_atom + " | grep TM-score | awk '{print $2}' "
                # print(cmd)
                mmalign_contents = os.popen(cmd).read().split('\n')
                mmalign_score = float(mmalign_contents[2].rstrip('\n'))

                model_result = f'{outputdir}/{dimer}/{method}/result_model_{i}_multimer.pkl'
                with open(model_result, 'rb') as f:
                    prediction_result = pickle.load(f)

                    dist_bins = np.append(0, prediction_result["distogram"]["bin_edges"])
                    dist_logits = prediction_result["distogram"]["logits"]
                    dist_mtx = dist_bins[dist_logits.argmax(-1)]
                    contact_mtx = nn.functional.softmax(torch.FloatTensor(dist_logits))[:, :, dist_bins < 8].sum(
                        -1).numpy()
                    contact_mtx[contact_mtx > 1] = 1

                    topl2 = get_range_top_of_contact_map(contact_mtx, 'long', 'l2')
                    number_of_topl2 = np.sum(topl2)
                    sum_prob_of_topl2 = np.sum(contact_mtx * topl2)
                    avg_prob_l2 = sum_prob_of_topl2 / number_of_topl2

                    topl5 = get_range_top_of_contact_map(contact_mtx, 'long', 'l5')
                    number_of_topl5 = np.sum(topl5)
                    sum_prob_of_topl5 = np.sum(contact_mtx * topl5)
                    avg_prob_l5 = sum_prob_of_topl5 / number_of_topl5

                    sum_dict = {'dimer': dimer,
                                'model': f"{method}_{i}",
                                'dockq': float(score.rstrip('\n')),
                                'tmscore1': tmscore1,
                                'tmscore2': tmscore2,
                                'mmalign': mmalign_score,
                                'lddt': np.mean(prediction_result['plddt']),
                                'avg_l2_prob': avg_prob_l2,
                                'avg_l5_prob': avg_prob_l5,
                                'ptmscore': float(prediction_result['ptm'])}

                    corr_summary_all = corr_summary_all.append(pd.DataFrame(sum_dict, index=[model_count]))

                    model_count += 1

                    dimer_result['dockq'] += float(score.rstrip('\n'))
                    dimer_result['tmscore1'] += tmscore1
                    dimer_result['tmscore2'] += tmscore2
                    dimer_result['lddt'] += np.mean(prediction_result['plddt'])
                    dimer_result['avg_l2_prob'] += avg_prob_l2
                    dimer_result['avg_l5_prob'] += avg_prob_l5
                    dimer_result['ptmscore'] += float(prediction_result['ptm'])

            dimer_result['dockq'] = dimer_result['dockq'] / 5
            dimer_result['tmscore1'] = dimer_result['tmscore1'] / 5
            dimer_result['tmscore2'] = dimer_result['tmscore2'] / 5
            dimer_result['mmalign'] = dimer_result['mmalign'] / 5
            dimer_result['lddt'] = dimer_result['lddt'] / 5
            dimer_result['avg_l2_prob'] = dimer_result['avg_l2_prob'] / 5
            dimer_result['avg_l5_prob'] = dimer_result['avg_l5_prob'] / 5
            dimer_result['ptmscore'] = dimer_result['ptmscore'] / 5
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

    pairwise_qa = Pairwise_dockq_qa(dockq_program)

    methods = ['no_templates',
               'original_template_pipeline',
               'custom_template_file',
               'custom_template_with_alphafold_models']

    for method in methods:
        corr_summary_all = pd.DataFrame(columns=['dockq', 'lddt', 'avg_l2_prob', 'avg_l5_prob', 'tmscore1', 'tmscore2', 'mmalign'])
        corr_summary_dimer = pd.DataFrame(
            columns=['dockq', 'lddt', 'avg_l2_prob', 'avg_l5_prob', 'tmscore1', 'tmscore2', 'mmalign'])

        process_list = []
        for line in open(args.dimerlist).readlines():
            chain1, chain2 = line.rstrip('\n').split()
            process_list.append([f"{chain1}_{chain2}", args.atom_dir, args.output_dir + '/' + method,
                                 args.tmpdir, dockq_program, tmscore_program, mmalign_program])
        pool = Pool(processes=30)
        results = pool.map(evaluate_dimer, process_list)
        pool.close()
        pool.join()

        for result in results:
            corr_summary_all = corr_summary_all.append(result['corr_summary_all_pd'])
            corr_summary_dimer = corr_summary_dimer.append(result['corr_summary_dimer_pd'])

        corr_summary_all.to_csv(f'{method}_all.csv')
        corr_summary_dimer.to_csv(f'{method}_avg.csv')
