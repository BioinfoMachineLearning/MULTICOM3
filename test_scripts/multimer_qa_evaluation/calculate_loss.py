import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.quaternary_structure_evaluation.pairwise_dockq import *
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir
import pathlib

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


def convert_ranking_to_df(infile, index_col='index', ranked_field="", is_csv=True, ascending=False):
    df = None
    # af=plddt_avg, enqa=score, native_scores='gdtscore
    if is_csv:
        if index_col is None:
            df = pd.read_csv(infile, index_col=index_col)
        else:
            df = pd.read_csv(infile)
        if ranked_field != 'score':
            df['score'] = df[ranked_field]
        df = df.sort_values(by=['score'], ascending=ascending)
    else:
        models = []
        scores = []
        for line in open(infile):
            line = line.rstrip('\n')
            contents = line.split()
            if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or \
                    contents[0] == "END":
                continue
            contents = line.split()
            model = contents[0]
            score = contents[1]
            if model.find('BML_CASP15') > 0:
                model = pathlib.Path(model).name
            if model.find('.pdb') < 0:
                model = model + '.pdb'
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

    gdt_loss_and_corr = {}
    lddt_loss_and_corr = {}

    for method in ['af', 'bfactor', 'multieva', 'af_pairwise_avg', 'bfactor_pairwise_avg']:#, 'foldseek', 'enqa', 'dproq_dockq', 'dproq_evalue']:
        gdt_loss_and_corr[method] = dict(corr=[],
                                         loss=[])
        lddt_loss_and_corr[method] = dict(corr=[],
                                          loss=[])

    targetname = []
    for target in os.listdir(args.inputdir):
        print(f"Processing {target}")
        result_csv = f"{args.inputdir}/{target}/native_scores.csv"
        if not os.path.exists(result_csv):
            print(f"The native scores for {target} is missing: {result_csv}")
            continue

        ranking_csvs = {'af': f"{args.inputdir}/{target}/alphafold_ranking.csv",
                        'bfactor': f"{args.inputdir}/{target}/bfactor_ranking.csv",
                        #'dproq_dockq': f"{args.inputdir}/{target}/dproq_ranking_dockq.csv",
                        #'dproq_evalue': f"{args.inputdir}/{target}/dproq_ranking_evalue.csv",
                        #'foldseek': f"{args.inputdir}/{target}/foldseek_qa.csv",
                        'multieva': f"{args.inputdir}/{target}/multieva.csv",
                        "af_pairwise_avg": f"{args.inputdir}/{target}/pairwise_af_avg.ranking",
                        "bfactor_pairwise_avg": f"{args.inputdir}/{target}/pairwise_bfactor_avg.ranking"}
                        #'enqa': f"{args.inputdir}/{target}/enqa_ranking.csv"}

        find_all_res = True
        for method in ranking_csvs:
            if not os.path.exists(ranking_csvs[method]):
                print(f"cannot find {ranking_csvs[method]}")
                find_all_res = False
                break

        if not find_all_res:
            continue

        native_gdt_scores = pd.read_csv(result_csv)
        native_gdt_scores = native_gdt_scores.sort_values(by=['tmscore'], ascending=False)
        native_gdt_scores.reset_index(inplace=True)


        # print(native_gdt_scores)
        for method in ranking_csvs:
            print(ranking_csvs[method])
            if method == 'af':
                global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
                                                      ranked_field='confidence',
                                                      is_csv=True)
            elif method == 'bfactor':
                global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
                                                      ranked_field='bfactor',
                                                      is_csv=True)
            elif method == 'dproq_dockq' or method == 'dproq_evalue':
                global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
                                                      ranked_field='score',
                                                      is_csv=True)
            elif method == "foldseek":
                global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
                                                      ranked_field='score2',
                                                      is_csv=True)
            elif method == "multieva":
                global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
                                                      ranked_field='MMalign score',
                                                      is_csv=True, index_col=None)
                print(global_scores)
                global_scores['model'] = global_scores['Name'] + '.pdb'
            elif method == 'enqa':
                global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
                                                      ranked_field='score',
                                                      is_csv=True)
            elif method == "af_pairwise_avg" or method == "bfactor_pairwise_avg":
                global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
                                                      ranked_field='avg_score',
                                                      is_csv=True)
            else:
                global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
                                                      is_csv=False)
            # print(global_scores)
            native_gdt_scores_corr = global_scores.merge(native_gdt_scores, on=['model'], how='inner')

            corr, _ = pearsonr(np.array(native_gdt_scores_corr.tmscore), np.array(native_gdt_scores_corr.score))
            gdt_loss_and_corr[method]['corr'] += [corr]
            top1_model = global_scores.model[0]
            # print(top1_model)
            top1_model_score = float(native_gdt_scores[native_gdt_scores.model == top1_model].tmscore)
            best_model_score = native_gdt_scores.tmscore[0]
            gdt_loss = best_model_score - top1_model_score
            gdt_loss_and_corr[method]['loss'] += [gdt_loss]

            # print(corr)
            # print(gdt_loss)



        targetname += [target]

    methods = []
    gdt_losses = []
    gdt_corrs = []
    for method in gdt_loss_and_corr:
        print(
            f"{method}: \n"
            f"gdt loss:{np.mean(np.array(gdt_loss_and_corr[method]['loss']))}\n"
            f"gdt correlation:{np.mean(np.array(gdt_loss_and_corr[method]['corr']))}")
        df = pd.DataFrame({'target': targetname,
                           'gdt_loss': gdt_loss_and_corr[method]['loss'],
                           'gdt_corr': gdt_loss_and_corr[method]['corr']})
        df.to_csv(f"{method}_eva_per_target.csv")
        gdt_losses += [np.mean(np.array(gdt_loss_and_corr[method]['loss']))]
        gdt_corrs += [np.mean(np.array(gdt_loss_and_corr[method]['corr']))]
        methods += [method]

    result_df = pd.DataFrame({'method': methods, 'gdt_loss': gdt_losses, 'gdt_corr': gdt_corrs})
    result_df.to_csv('summary.csv')
    print(pd.DataFrame({'method': methods, 'gdt_loss': gdt_losses, 'gdt_corr': gdt_corrs}))
