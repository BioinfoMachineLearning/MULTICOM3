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

def read_qa_txt_as_df(infile):
    models = []
    scores = []
    for line in open(infile):
        line = line.rstrip('\n')
        contents = line.split()
        if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or \
                contents[0] == "END":
            continue
        model, score = line.split()
        models += [model]
        scores += [float(score)]
    df = pd.DataFrame({'model': models, 'score': scores})
    df = df.sort_values(by=['score'], ascending=False)
    df.reset_index(inplace=True)
    return df

def write_qa_file(outfile, models, scores):
    with open(outfile, 'w') as fw:
        for model, score in zip(models, scores):
            fw.write(f"{model} {score} {score}\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    for target in os.listdir(args.outdir):
        workdir = args.outdir + '/' + target

        alphafold_ranking = workdir + '/alphafold_ranking.csv'

        af_df = pd.read_csv(alphafold_ranking)

        write_qa_file(f"{workdir}/{target}QA1_1", list(af_df['model']), list(af_df['confidence']))
        
        pairwise_ranking = workdir + '/pairwise_ranking.csv'

        pairwise_df = pd.read_csv(pairwise_ranking)

        write_qa_file(f"{workdir}/{target}QA2_1", list(pairwise_df['model']), list(pairwise_df['pairwise_score']))

        avg_ranking_df = pairwise_df.merge(af_df, how="inner", on='model')
        
        for index, weight in enumerate([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]):
            models, scores = [], []
            for i in range(len(avg_ranking_df)):
                pairwise_score = float(avg_ranking_df.loc[i, 'pairwise_score'])
                alphafold_score = float(avg_ranking_df.loc[i, 'confidence'])
                avg_score = weight * alphafold_score + (1-weight) * pairwise_score
                models += [avg_ranking_df.loc[i, 'model']]
                scores += [avg_score]
            write_qa_file(f"{workdir}/{target}QA2{index+1}_1", models, scores)

        
        af_pairwise_ranking = workdir + '/pairwise_af_ranking.csv'
        pairwise_df = pd.read_csv(af_pairwise_ranking)

        write_qa_file(f"{workdir}/{target}QA3_1", list(pairwise_df['model']), list(pairwise_df['pairwise_score']))

    
