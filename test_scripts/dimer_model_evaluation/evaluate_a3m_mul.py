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

    dimer, atomdir, outputdir, tmpdir, dockq_program = inparams

    chain1, chain2 = dimer.split('_')

    chain1_atom = f"{atomdir}/{chain1}.atom"

    chain2_atom = f"{atomdir}/{chain2}.atom"

    combine_atom = f"{outputdir}/{dimer}/{chain1}_{chain2}.atom"

    combine_pdb(chain1_atom, chain2_atom, combine_atom)

    ptm_result = {'dimer': dimer}
    plddt_result = {'dimer': dimer}
    topl2_result = {'dimer': dimer}
    topl5_result = {'dimer': dimer}
    pairwise_result = {'dimer': dimer}
    average_result = {'dimer': dimer}
    tmpdir = tmpdir + '/' + dimer

    model_count = 0
    target_model_dockq = pd.DataFrame(columns=['dockq', 'lddt', 'avg_l2_prob', 'avg_l5_prob', 'tmscore'])
    corr_summary = pd.DataFrame(columns=['dockq', 'lddt', 'avg_l2_prob', 'avg_l5_prob', 'tmscore'])
    all_models_path = {}
    #print(methods)
    for method in methods:
        dimer_chainA = chain1[4]
        dimer_chainB = chain2[4]
        dockQ_scores = []
        method_models_path = {}

        if not complete_result(f'{args.output_dir}/{dimer}/{method}'):
            lddt_top1_score = -1
            l2_top1_score = -1
            l5_top1_score = -1
            tmscore_top1_score = -1
            pairwise_top1_score = -1
            avg_top1_score = -1
            bestof5_score = -1
        else:
            bfail = False
            method_model_dockq = pd.DataFrame(columns=['dockq', 'lddt', 'avg_l2_prob', 'avg_l5_prob', 'tmscore'])
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
                                'lddt': np.mean(prediction_result['plddt']),
                                'avg_l2_prob': avg_prob_l2,
                                'avg_l5_prob': avg_prob_l5,
                                'tmscore': float(prediction_result['ptm'])}

                    corr_summary = corr_summary.append(pd.DataFrame(sum_dict, index=[model_count]))
                    target_model_dockq = target_model_dockq.append(pd.DataFrame(sum_dict, index=[model_count]))
                    method_model_dockq = method_model_dockq.append(pd.DataFrame(sum_dict, index=[i]))
                    model_count += 1
                    all_models_path[model] = f"{method}_{i}"
                    method_models_path[model] = f"{method}_{i}"

            if bfail:
                lddt_top1_score = -1
                l2_top1_score = -1
                l5_top1_score = -1
                tmscore_top1_score = -1
                avg_top1_score = -1
                bestof5_score = -1
            else:
                lddt_argmax = method_model_dockq[method_model_dockq.lddt == method_model_dockq.lddt.max()].index
                lddt_top1_score = method_model_dockq.loc[lddt_argmax[0]].dockq

                l2_argmax = method_model_dockq[
                    method_model_dockq.avg_l2_prob == method_model_dockq.avg_l2_prob.max()].index
                l2_top1_score = method_model_dockq.loc[l2_argmax[0]].dockq

                l5_argmax = method_model_dockq[
                    method_model_dockq.avg_l5_prob == method_model_dockq.avg_l5_prob.max()].index
                l5_top1_score = method_model_dockq.loc[l5_argmax[0]].dockq

                tmscore_argmax = method_model_dockq[
                    method_model_dockq.tmscore == method_model_dockq.tmscore.max()].index
                tmscore_top1_score = method_model_dockq.loc[tmscore_argmax[0]].dockq

                clean_dir(tmpdir)
                for model in method_models_path:
                    os.system(f"cp {model} {tmpdir}/{method_models_path[model]}")
                pairwise_pd = pairwise_qa.run(tmpdir)
                pairwise_top1_score = float(
                    method_model_dockq.loc[method_model_dockq.model == pairwise_pd.loc[0].model].dockq)

                bestof5_score = np.max(np.array(dockQ_scores))

                average_pd = method_model_dockq.merge(
                    pairwise_pd,
                    how="inner",  # takes the intersection
                    on="model",  # merges on species identifiers
                    suffixes=("_1", "_2")
                )

                average_scores = []
                for index, row in average_pd.iterrows():
                    average_scores += [(float(row['tmscore']) + float(row['lddt'])/100 + float(row['pairwise_score']))/3]
                average_pd['avg_score'] = average_scores
                # print(average_pd)
                average_pd_sort = average_pd.sort_values(by=['avg_score'], ascending=False, ignore_index=True)
                avg_top1_score = float(
                    method_model_dockq.loc[method_model_dockq.model == average_pd_sort.loc[0].model].dockq)

        ptm_result[f"{method}_top1"] = tmscore_top1_score
        ptm_result[f"{method}_bestof5"] = bestof5_score
        plddt_result[method] = lddt_top1_score
        plddt_result[f"{method}_bestof5"] = bestof5_score
        topl2_result[method] = l2_top1_score
        topl2_result[f"{method}_bestof5"] = bestof5_score
        topl5_result[method] = l5_top1_score
        topl5_result[f"{method}_bestof5"] = bestof5_score
        pairwise_result[method] = pairwise_top1_score
        pairwise_result[f"{method}_bestof5"] = bestof5_score
        average_result[method] = avg_top1_score
        average_result[f"{method}_bestof5"] = bestof5_score

    if model_count > 0:

        sort_lddt_pd = target_model_dockq.sort_values(by=['lddt'], ascending=False, ignore_index=True)
        plddt_result['top1'] = sort_lddt_pd.loc[0].dockq
        plddt_result['bestof5'] = np.max(np.array([sort_lddt_pd.loc[n].dockq for n in range(5)]))

        sort_l2_pd = target_model_dockq.sort_values(by=['avg_l2_prob'], ascending=False, ignore_index=True)
        topl2_result['top1'] = sort_l2_pd.loc[0].dockq
        topl2_result['bestof5'] = np.max(np.array([sort_l2_pd.loc[n].dockq for n in range(5)]))

        sort_l5_pd = target_model_dockq.sort_values(by=['avg_l5_prob'], ascending=False, ignore_index=True)
        topl5_result['top1'] = sort_l5_pd.loc[0].dockq
        topl5_result['bestof5'] = np.max(np.array([sort_l5_pd.loc[n].dockq for n in range(5)]))

        sort_tmscore_pd = target_model_dockq.sort_values(by=['tmscore'], ascending=False, ignore_index=True)
        ptm_result['top1'] = sort_tmscore_pd.loc[0].dockq
        ptm_result['bestof5'] = np.max(np.array([sort_tmscore_pd.loc[n].dockq for n in range(5)]))

        clean_dir(tmpdir)
        for model in all_models_path:
            os.system(f"cp {model} {tmpdir}/{all_models_path[model]}")
        pairwise_pd = pairwise_qa.run(tmpdir)
        pairwise_result['top1'] = float(
            target_model_dockq.loc[target_model_dockq.model == pairwise_pd.loc[0].model].dockq)
        pairwise_result['bestof5'] = np.max(np.array(
            [float(target_model_dockq.loc[target_model_dockq.model == pairwise_pd.loc[n].model].dockq) for n in
             range(5)]))

        average_pd = target_model_dockq.merge(
            pairwise_pd,
            how="inner",  # takes the intersection
            on="model",  # merges on species identifiers
            suffixes=("_1", "_2")
        )
        average_scores = []
        for index, row in average_pd.iterrows():
            average_scores += [(float(row['tmscore']) + float(row['lddt']) / 100 + float(row['pairwise_score'])) / 3]
        average_pd['avg_score'] = average_scores
        average_pd_sort = average_pd.sort_values(by=['avg_score'], ascending=False, ignore_index=True)
        average_result['top1'] = float(
            target_model_dockq.loc[target_model_dockq.model == average_pd_sort.loc[0].model].dockq)
        average_result['bestof5'] = np.max(np.array(
            [float(target_model_dockq.loc[target_model_dockq.model == average_pd_sort.loc[n].model].dockq) for n in
             range(5)]))

        sort_target_model = target_model_dockq.sort_values(by=['dockq'], ascending=False, ignore_index=True)
        plddt_result['best_model'] = float(sort_target_model.loc[0].dockq)
        topl2_result['best_model'] = float(sort_target_model.loc[0].dockq)
        topl5_result['best_model'] = float(sort_target_model.loc[0].dockq)
        ptm_result['best_model'] = float(sort_target_model.loc[0].dockq)
        pairwise_result['best_model'] = float(sort_target_model.loc[0].dockq)
        average_result['best_model'] = float(sort_target_model.loc[0].dockq)
        print(average_pd)

    return {'ptm_result': ptm_result,
            'topl2_result': topl2_result,
            'topl5_result': topl5_result,
            'plddt_result': plddt_result,
            'pairwise_result': pairwise_result,
            'average_result': average_result,
            'corr_summary_pd': corr_summary}


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dimerlist', type=is_file, required=True)
    parser.add_argument('--atom_dir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=is_dir, required=True)
    parser.add_argument('--tmpdir', type=is_dir, required=True)

    args = parser.parse_args()

    dockq_program = '/home/bml_casp15/BML_CASP15/tools/DockQ/DockQ.py'

    pairwise_qa = Pairwise_dockq_qa(dockq_program)

    methods = ['alphafold',
               'uniclust_oxmatch_a3m',
                'pdb_interact_uniref_a3m',
               'species_interact_uniref_a3m',
               'uniprot_distance_uniref_a3m',
                'string_interact_uniref_a3m',
               # 'geno_dist_uniref_a3m',
               'pdb_interact_uniref_sto',
               'species_interact_uniref_sto',
               'uniprot_distance_uniref_sto',
               'string_interact_uniref_sto',
               # 'geno_dist_uniref_sto',
               'pdb_interact_uniprot_sto',
               'species_interact_uniprot_sto',
               'uniprot_distance_uniprot_sto',
               'string_interact_uniprot_sto']

    column_names = []
    for method in methods:
        for top_model in ['top1', 'bestof5']:
            column_names += [f"{method}_{top_model}"]

    selection_methods = ['ptm', 'plddt', 'topl2_prec', 'topl5_prec']

    ptm_results = pd.DataFrame(columns=['dimer', 'top1', 'bestof5', 'best_model'].append(column_names))
    plddt_results = pd.DataFrame(columns=['dimer', 'top1', 'bestof5', 'best_model'].append(column_names))
    topl2_results = pd.DataFrame(columns=['dimer', 'top1', 'bestof5', 'best_model'].append(column_names))
    topl5_results = pd.DataFrame(columns=['dimer', 'top1', 'bestof5', 'best_model'].append(column_names))
    pairwise_results = pd.DataFrame(columns=['dimer', 'top1', 'bestof5', 'best_model'].append(column_names))
    average_results = pd.DataFrame(columns=['dimer', 'top1', 'bestof5', 'best_model'].append(column_names))
    corr_summary = pd.DataFrame(columns=['dockq', 'lddt', 'avg_l2_prob', 'avg_l5_prob', 'tmscore'])

    #dimers = os.listdir(args.output_dir)
    process_list = []
    # for dimer_count in range(len(dimers)):
    #     dimer = dimers[dimer_count]
    #     process_list.append([dimer, args.atom_dir, args.output_dir, args.tmpdir, dockq_program])
    for line in open(args.dimerlist).readlines():
        dimer = line.rstrip('\n')
        process_list.append([dimer, args.atom_dir, args.output_dir, args.tmpdir, dockq_program])

    pool = Pool(processes=10)
    results = pool.map(evaluate_dimer, process_list)
    pool.close()
    pool.join()

    dimer_count = 0
    for result in results:
        ptm_results = ptm_results.append(pd.DataFrame(result['ptm_result'], index=[dimer_count]))
        topl2_results = topl2_results.append(pd.DataFrame(result['topl2_result'], index=[dimer_count]))
        topl5_results = topl5_results.append(pd.DataFrame(result['topl5_result'], index=[dimer_count]))
        plddt_results = plddt_results.append(pd.DataFrame(result['plddt_result'], index=[dimer_count]))
        pairwise_results = pairwise_results.append(pd.DataFrame(result['pairwise_result'], index=[dimer_count]))
        average_results = average_results.append(pd.DataFrame(result['average_result'], index=[dimer_count]))
        corr_summary = corr_summary.append(result['corr_summary_pd'])
        dimer_count += 1

    ptm_results.to_csv('eva_tm.csv')
    topl2_results.to_csv('eva_l2.csv')
    topl5_results.to_csv('eva_l5.csv')
    plddt_results.to_csv('eva_lddt.csv')
    pairwise_results.to_csv('eva_pairwise.csv')
    average_results.to_csv('eva_avg.csv')
    corr_summary.to_csv('eva_summary.csv')
