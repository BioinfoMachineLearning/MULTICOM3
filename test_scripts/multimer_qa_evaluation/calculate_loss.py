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
    parser.add_argument('--inputdir', type=is_dir, required=True)
    args = parser.parse_args()

    all_af_plddt_corr = []
    all_af_plddt_loss = []

    all_af_confidence_corr = []
    all_af_confidence_loss = []

    all_af_ptm_corr = []
    all_af_ptm_loss = []

    all_pairwise_corr = []
    all_pairwise_loss = []


    targetname = []
    for target in os.listdir(args.inputdir):
        print(f"Processing {target}")
        result_csv = f"{args.inputdir}/{target}/native_scores.csv"
        if not os.path.exists(result_csv):
            print(f"The native scores for {target} is missing: {result_csv}")
            continue

        alphafold_ranking_csv = f"{args.inputdir}/{target}/alphafold_ranking.csv"
        pairwise_ranking_csv = f"{args.inputdir}/{target}/pairwise_ranking.csv"

        if not os.path.exists(alphafold_ranking_csv) or not os.path.exists(pairwise_ranking_csv):
            continue

        native_scores = pd.read_csv(result_csv).sort_values(by=['dockq_score'], ascending=False)
        # print(native_scores)
        global_scores = pd.read_csv(alphafold_ranking_csv).sort_values(by=['ptm'], ascending=False)
        native_scores_corr = global_scores.merge(native_scores, on=['model'], how='inner')
        corr, _ = pearsonr(np.array(native_scores_corr.dockq_score), np.array(native_scores_corr.ptm))
        top1_model = global_scores.model[0]
        top1_model_score = float(native_scores[native_scores.model == top1_model].dockq_score)
        best_model_score = native_scores.dockq_score[0]
        # print(top1_model_score)
        # print(best_model_score)
        dockq_loss = best_model_score - top1_model_score
        all_af_ptm_loss += [dockq_loss]
        all_af_ptm_corr += [corr]
        ptm_top_1_model = top1_model

        global_scores = pd.read_csv(alphafold_ranking_csv).sort_values(by=['plddt_avg'], ascending=False)
        native_scores_corr = global_scores.merge(native_scores, on=['model'], how='inner')
        corr, _ = pearsonr(np.array(native_scores_corr.dockq_score), np.array(native_scores_corr.plddt_avg))
        top1_model = global_scores.model[0]
        top1_model_score = float(native_scores[native_scores.model == top1_model].dockq_score)
        best_model_score = native_scores.dockq_score[0]
        dockq_loss = best_model_score - top1_model_score
        all_af_plddt_loss += [dockq_loss]
        all_af_plddt_corr += [corr]
        plddt_top1_model = top1_model

        global_scores = pd.read_csv(alphafold_ranking_csv).sort_values(by=['confidence'], ascending=False)
        native_scores_corr = global_scores.merge(native_scores, on=['model'], how='inner')
        corr, _ = pearsonr(np.array(native_scores_corr.dockq_score), np.array(native_scores_corr.confidence))
        top1_model = global_scores.model[0]
        top1_model_score = float(native_scores[native_scores.model == top1_model].dockq_score)
        best_model_score = native_scores.dockq_score[0]
        dockq_loss = best_model_score - top1_model_score
        all_af_confidence_loss += [dockq_loss]
        all_af_confidence_corr += [corr]
        confidence_top1_model = top1_model


        global_scores = pd.read_csv(pairwise_ranking_csv).sort_values(by=['pairwise_score'], ascending=False)
        native_scores_corr = global_scores.merge(native_scores, on=['model'], how='inner')
        corr, _ = pearsonr(np.array(native_scores_corr.dockq_score), np.array(native_scores_corr.pairwise_score))
        top1_model = global_scores.model[0]
        top1_model_score = float(native_scores[native_scores.model == top1_model].dockq_score)
        best_model_score = native_scores.dockq_score[0]
        dockq_loss = best_model_score - top1_model_score
        all_pairwise_loss += [dockq_loss]
        all_pairwise_corr += [corr]

        targetname += [target]

        if ptm_top_1_model == plddt_top1_model == confidence_top1_model:
            print("same alphafold ranking result")

    print(f"alphafold_plddt: \ndockq loss:{np.mean(np.array(all_af_plddt_loss))}\ndockq correlation:{np.mean(np.array(all_af_plddt_corr))}")
    df = pd.DataFrame({'target': targetname,
                       'dockq_loss': all_af_plddt_loss,
                       'dockq_corr': all_af_plddt_corr})
    df.to_csv("alphafold_plddt_eva_per_dimer.csv")

    print(f"alphafold_ptm: \ndockq loss:{np.mean(np.array(all_af_ptm_loss))}\ndockq correlation:{np.mean(np.array(all_af_ptm_corr))}")
    df = pd.DataFrame({'target': targetname,
                       'dockq_loss': all_af_ptm_loss,
                       'dockq_corr': all_af_ptm_corr})
    df.to_csv("alphafold_ptm_eva_per_dimer.csv")

    print(f"alphafold_confidence: \ndockq loss:{np.mean(np.array(all_af_confidence_loss))}\ndockq correlation:{np.mean(np.array(all_af_confidence_corr))}")
    df = pd.DataFrame({'target': targetname,
                       'dockq_loss': all_af_confidence_loss,
                       'dockq_corr': all_af_confidence_corr})
    df.to_csv("alphafold_confidence_eva_per_dimer.csv")

    print(f"pairwise: \ndockq loss:{np.mean(np.array(all_pairwise_loss))}\ndockq correlation:{np.mean(np.array(all_pairwise_corr))}")
    df = pd.DataFrame({'target': targetname,
                       'dockq_loss': all_pairwise_loss,
                       'dockq_corr': all_pairwise_corr})
    df.to_csv("pairwise_eva_per_dimer.csv")