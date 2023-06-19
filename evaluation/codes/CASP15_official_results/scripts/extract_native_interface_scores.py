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
    parser.add_argument('--infile', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    targets = []
    df = pd.read_csv(args.infile)
    for trg in df['trg']:
        if trg not in targets:
            targets += [trg]

    for target in targets:
        
        targetname = target.rstrip('.pdb')

        cmd = f"grep {target} {args.infile}"
        print(cmd)
        contents = os.popen(cmd).readlines()

        data_dict = {'model': [], 'qscore': [], 'dockq': [], 'tmscore': []}
        for line in contents:
            tmp = line.split(',')
            model = tmp[1]
            qscore = tmp[3]
            dockq = tmp[9]
            tmscore = tmp[10]
            data_dict['model'] += [model]
            data_dict['qscore'] += [qscore]
            data_dict['dockq'] += [dockq]
            data_dict['tmscore'] += [tmscore]
        pd.DataFrame(data_dict).to_csv(args.outdir + '/' + targetname + '.csv')

