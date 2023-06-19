from calendar import c
from curses import raw
import os, sys, argparse, time
from pydoc import doc
from multiprocessing import Pool
from tqdm import tqdm
import random
import numpy as np
from scipy.stats.stats import pearsonr
import pandas as pd

def read_qa_txt_as_df(targetname, infile):
    models = []
    scores = []
    for line in open(infile):
        line = line.rstrip('\n')
        # print(line)
        if len(line) == 0:
            continue
        contents = line.split()
        # if contents[0] == "QMODE":
        #     if float(contents[1]) == 2:
        #         return None

        if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or \
                contents[0] == "END" or contents[0] == "REMARK":
            continue

        contents = line.split()
        
        model, score1 = contents[0], contents[1]
        if score1 == "X":
            continue

        models += [model]
        scores += [float(score1)]
    df = pd.DataFrame({'model': models, 'score': scores})
    df = df.sort_values(by=['score'], ascending=False)
    df.reset_index(inplace=True)
    return df

def read_mmalign(infile):
    # print(infile)
    for line in open(infile):
        line = line.rstrip('\n')
        if len(line) == 0:
            continue
        if line.split()[0] == 'TM-score=':
            # print(line)
            if line.find('Structure_1') >= 0:
                return float(line.split()[1])
    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--tarball', type=str, required=True)
    parser.add_argument('--scoredir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    group_ids = []
    
    for target in sorted(os.listdir(args.tarball)):
        for prediction in os.listdir(args.tarball + '/' + target):
            if prediction.find('QA') < 0:
                continue
            groupid = prediction[prediction.find('QA')+2:prediction.find('_')]
            if groupid not in group_ids:
                group_ids += [groupid]

    group_ids = sorted(group_ids)
    print(group_ids)    
    group_res = {}
    for group_id in group_ids:
        # print(group_id)
        corrs = []
        losses = []
        selected_top1_scores = []
        global_qa = False
        for target in sorted(os.listdir(args.scoredir)):
            prediction = f"{args.tarball}/{target}/{target}QA{group_id}_1"
            print(prediction)
            if not os.path.exists(prediction):
                corrs += ["0"]
                losses += ["10"]
                selected_top1_scores += ["0"]
                continue
            
            print(prediction)
            pred_df = read_qa_txt_as_df(target, prediction)
            # print(pred_df)
            if pred_df is None:
                corrs += ["0"]
                losses += ["11"]
                selected_top1_scores += ["0"]
                continue

            scores_filt = []
            scores_true = []
            true_scores_dict = {}
            for i in range(len(pred_df)):
                model = pred_df.loc[i, 'model']
                scorefile = args.scoredir + '/' + target + '/' + model.rstrip('.pdb') + '_filtered_out'
                if not os.path.exists(scorefile):
                    print(f"Excluding {model} in {target}")
                    continue

                true_score = read_mmalign(scorefile)
                scores_filt += [pred_df.loc[i, 'score']]
                scores_true += [true_score]
                true_scores_dict[model] = true_score

            # print(scores_filt)
            # print(scores_true)

            if len(scores_filt) == 0:
                corrs += ["0"]
                losses += ["12"]
                selected_top1_scores += ["0"]
                continue

            corr = pearsonr(np.array(scores_filt), np.array(scores_true))[0]
            # print(corr)

            top1_model = pred_df.loc[0, 'model']
            # print(top1_model)
            if top1_model not in true_scores_dict:
                corrs += ["0"]
                losses += ["13"]
                selected_top1_scores += ["0"]
                raise Exception(f"Cannot find the {scorefile} for {top1_model}!")

            # print(np.max(np.array(scores_true)))
            # print(true_scores_dict[top1_model])
            loss = float(np.max(np.array(scores_true))) - float(true_scores_dict[top1_model])

            corrs += [str(corr)]
            losses += [str(loss)]
            selected_top1_scores += [str(true_scores_dict[top1_model])]

            global_qa = True

        if global_qa:
            group_res[group_id] = dict(corrs=corrs, losses=losses, selected_scores=selected_top1_scores)

    group_ids = [key for key in group_res]
    print('\t'.join(group_ids))
    
    targets = [target for target in sorted(os.listdir(args.scoredir))]
    print('\t'.join(targets))

    for i in range(len(os.listdir(args.scoredir))):
        contents = []
        for group_id in group_res:
            contents += [group_res[group_id]['corrs'][i]]
            contents += [group_res[group_id]['losses'][i]]
            contents += [group_res[group_id]['selected_scores'][i]]
        print('\t'.join(contents))




    
