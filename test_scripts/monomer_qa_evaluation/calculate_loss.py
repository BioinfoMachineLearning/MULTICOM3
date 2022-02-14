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

def cal_scores(tmscore_program, inpdb, nativepdb):
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $3}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[2].rstrip('\n'))
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep GDT-score | awk '{print $3}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    gdtscore = float(tmscore_contents[0].rstrip('\n'))
    return tmscore, gdtscore


def convert_ranking_to_df(infile):
    df = None
    if infile.find('alphafold_ranking.csv') > 0:
        df = pd.read_csv(infile)
        df = df.sort_values(by=['plddt_avg'], ascending=False)
    elif infile.find('enqa_ranking.csv') > 0:
        df = pd.read_csv(infile)
        df = df.sort_values(by=['score'], ascending=False)
        df.reset_index(inplace=True)
    elif infile.find('native_scores.csv') > 0:
        df = pd.read_csv(infile)
        df = df.sort_values(by=['gdtscore'], ascending=False)
        df.reset_index(inplace=True)
    else:
        models = []
        scores = []
        for line in open(infile):
            line = line.rstrip('\n')
            contents = line.split()
            if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or contents[0] == "END":
                continue
            model, score = line.split()
            models += [model]
            scores += [float(score)]
        df = pd.DataFrame({'model': models, 'score': scores})
        df = df.sort_values(by=['score'], ascending=False)
        df.reset_index(inplace=True)
    return df

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--inputdir', type=is_dir, required=True)
    args = parser.parse_args()

    all_af_gdt_corr = []
    all_af_gdt_loss = []
    all_af_lddt_corr = []
    all_af_lddt_loss = []

    all_enqa_gdt_corr = []
    all_enqa_gdt_loss = []
    all_enqa_lddt_corr = []
    all_enqa_lddt_loss = []

    all_pairwise_gdt_corr = []
    all_pairwise_gdt_loss = []
    all_pairwise_lddt_corr = []
    all_pairwise_lddt_loss = []

    targetname = []
    for target in os.listdir(args.inputdir):
        print(f"Processing {target}")
        result_csv = f"{args.inputdir}/{target}/native_scores.csv"
        if not os.path.exists(result_csv):
            print(f"The native scores for {target} is missing: {result_csv}")
            continue

        alphafold_ranking_csv = f"{args.inputdir}/{target}/alphafold_ranking.csv"
        enqa_ranking_csv = f"{args.inputdir}/{target}/enqa_ranking.csv"
        pairwise_ranking_csv = f"{args.inputdir}/{target}/pairwise_ranking.tm"

        if not os.path.exists(alphafold_ranking_csv) or not os.path.exists(enqa_ranking_csv) or not os.path.exists(pairwise_ranking_csv):
            continue

        native_gdt_scores = convert_ranking_to_df(result_csv)

        global_scores = convert_ranking_to_df(alphafold_ranking_csv)
        native_gdt_scores_corr = global_scores.merge(native_gdt_scores, on=['model'], how='inner')

        corr, _ = pearsonr(np.array(native_gdt_scores_corr.gdtscore), np.array(native_gdt_scores_corr.plddt_avg))
        top1_model = global_scores.model[0]
        top1_model_score = float(native_gdt_scores[native_gdt_scores.model == top1_model].gdtscore)
        best_model_score = native_gdt_scores.gdtscore[0]
        gdt_loss = best_model_score - top1_model_score
        all_af_gdt_loss += [gdt_loss]
        all_af_gdt_corr += [corr]

        corr, _ = pearsonr(np.array(native_gdt_scores_corr.lddtscore), np.array(native_gdt_scores_corr.plddt_avg))
        top1_model = global_scores.model[0]
        top1_model_score = float(native_gdt_scores[native_gdt_scores.model == top1_model].lddtscore)
        best_model_score = native_gdt_scores.lddtscore[0]
        lddt_loss = best_model_score - top1_model_score
        all_af_lddt_loss += [lddt_loss]
        all_af_lddt_corr += [corr]


        global_scores = convert_ranking_to_df(enqa_ranking_csv)
        native_gdt_scores_corr = global_scores.merge(native_gdt_scores, on=['model'], how='inner')
        corr, _ = pearsonr(np.array(native_gdt_scores_corr.gdtscore), np.array(native_gdt_scores_corr.score))

        top1_model = global_scores.model[0]
        top1_model_score = float(native_gdt_scores[native_gdt_scores.model == top1_model].gdtscore)
        best_model_score = native_gdt_scores.gdtscore[0]
        gdt_loss = best_model_score - top1_model_score
        all_enqa_loss += [gdt_loss]
        all_enqa_corr += [corr]

        corr, _ = pearsonr(np.array(native_gdt_scores_corr.lddtscore), np.array(native_gdt_scores_corr.score))
        top1_model = global_scores.model[0]
        top1_model_score = float(native_gdt_scores[native_gdt_scores.model == top1_model].lddtscore)
        best_model_score = native_gdt_scores.lddtscore[0]
        lddt_loss = best_model_score - top1_model_score
        all_enqa_lddt_loss += [lddt_loss]
        all_enqa_lddt_corr += [corr]



        global_scores = convert_ranking_to_df(pairwise_ranking_csv)
        native_gdt_scores_corr = global_scores.merge(native_gdt_scores, on=['model'], how='inner')
        corr, _ = pearsonr(np.array(native_gdt_scores_corr.gdtscore), np.array(native_gdt_scores_corr.score))

        top1_model = global_scores.model[0]
        top1_model_score = float(native_gdt_scores[native_gdt_scores.model == top1_model].gdtscore)
        best_model_score = native_gdt_scores.gdtscore[0]
        gdt_loss = best_model_score - top1_model_score
        all_pairwise_loss += [gdt_loss]
        all_pairwise_corr += [corr]

        corr, _ = pearsonr(np.array(native_gdt_scores_corr.lddtscore), np.array(native_gdt_scores_corr.score))
        top1_model = global_scores.model[0]
        top1_model_score = float(native_gdt_scores[native_gdt_scores.model == top1_model].lddtscore)
        best_model_score = native_gdt_scores.lddtscore[0]
        lddt_loss = best_model_score - top1_model_score
        all_pairwise_lddt_loss += [lddt_loss]
        all_pairwise_lddt_corr += [corr]

        targetname += [target]

    print(f"alphafold: \ngdt loss:{np.mean(np.array(all_af_gdt_loss))}\ngdt correlation:{np.mean(np.array(all_af_gdt_corr))}")
    print(f"lddt loss:{np.mean(np.array(all_af_lddt_loss))}\ngdt correlation:{np.mean(np.array(all_af_lddt_corr))}")
    df = pd.DataFrame({'target': targetname,
                       'gdt_loss': all_af_gdt_loss,
                       'gdt_corr': all_af_gdt_corr,
                       'lddt_loss': all_af_lddt_loss,
                       'lddt_corr': all_af_lddt_corr})
    df.to_csv("alphafold_eva_per_target.csv")

    print(f"enqa: \ngdt loss:{np.mean(np.array(all_enqa_gdt_loss))}\ngdt correlation:{np.mean(np.array(all_enqa_gdt_corr))}")
    print(f"lddt loss:{np.mean(np.array(all_enqa_lddt_loss))}\ngdt correlation:{np.mean(np.array(all_enqa_lddt_corr))}")
    df = pd.DataFrame({'target': targetname,
                       'gdt_loss': all_enqa_gdt_loss,
                       'gdt_corr': all_enqa_gdt_corr,
                       'lddt_loss': all_enqa_lddt_loss,
                       'lddt_corr': all_enqa_lddt_corr})
    df.to_csv("enqa_eva_per_target.csv")


    print(f"pairwise: \ngdt loss:{np.mean(np.array(all_pairwise_gdt_loss))}\ngdt correlation:{np.mean(np.array(all_pairwise_gdt_corr))}")
    print(f"lddt loss:{np.mean(np.array(all_pairwise_lddt_loss))}\ngdt correlation:{np.mean(np.array(all_pairwise_lddt_corr))}")
    df = pd.DataFrame({'target': targetname,
                       'gdt_loss': all_pairwise_gdt_loss,
                       'gdt_corr': all_pairwise_gdt_corr,
                       'lddt_loss': all_pairwise_lddt_loss,
                       'lddt_corr': all_pairwise_lddt_corr})
    df.to_csv("pairwise_eva_per_target.csv")