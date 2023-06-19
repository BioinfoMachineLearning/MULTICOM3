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
        if contents[0].find(targetname) < 0:
            continue
        
        model, score1, score2 = contents[0], contents[1], contents[2]
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
            if line.find('Chain_2') >= 0:
                return float(line.split()[1])
    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--workdir', type=str, required=True)
    parser.add_argument('--scoredir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    af_groupid = 1
    pairwise_groupid = 2
    average_groupid = 3

    for target in os.listdir(args.scoredir):
        workdir = args.workdir + '/' + target + '/N1_quaternary_structure_evaluation'

        casp_file = f"{workdir}/pairwise_af_avg.ranking"
        if not os.path.exists(casp_file):
            raise Exception(f"Cannot find the result file: {casp_file}")

        if not os.path.exists(f"{args.outdir}/{target}"):
            os.makedirs(f"{args.outdir}/{target}")
            
        df = pd.read_csv(casp_file)
        with open(f"{args.outdir}/{target}/{target}QA{af_groupid}_1", 'w') as fw:
            for model, score in zip(df['model'], df['confidence']):
                fw.write(f"{model} {score} X\n")

        with open(f"{args.outdir}/{target}/{target}QA{pairwise_groupid}_1", 'w') as fw:
            for model, score in zip(df['model'], df['average_MMS']):
                fw.write(f"{model} {score} X\n")
    
        with open(f"{args.outdir}/{target}/{target}QA{average_groupid}_1", 'w') as fw:
            for model, score in zip(df['model'], df['avg_score']):
                fw.write(f"{model} {score} X\n")
    
