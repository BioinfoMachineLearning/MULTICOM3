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
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    groupnames = []
    targets = []
    for native_score_file in sorted(os.listdir(args.indir)):
        target = native_score_file.rstrip('o.csv')
        targets += [target]
        native_df = pd.read_csv(args.indir + '/' + native_score_file)
        for model in native_df['model']:
            groupname = model[model.find('TS'):model.find('TS')+5]
            if groupname not in groupnames:
                groupnames += [groupname]

        with open(args.outdir + '/' + target + '.results', 'w') as fw:
            for groupname in groupnames:
                write_contents = [groupname]
                for i in range(1, 6):
                    if target[0] == "T":
                        model = f"{target}{groupname}_{i}o"
                    else:
                        model = f"{target}{groupname}_{i}"
                    cmd = f"grep {model} {args.indir}/{native_score_file}"
                    print(cmd)
                    contents = os.popen(cmd).readlines()
                    tmscore = None
                    for line in contents:
                        tmp = line.split(',')
                        tmscore = tmp[4].rstrip('\n')
                        break
                    if tmscore is None:
                        continue
                    tmscore = str(round(float(tmscore),3))
                    write_contents += [tmscore]
                if len(write_contents) == 1:
                    continue
                fw.write('\t'.join(write_contents) + '\n')

    for groupname in groupnames:
        with open(args.outdir + '/' + groupname + '.results', 'w') as fw:
            for target in targets:
                write_contents = [target]
                cmd = f"grep {groupname} {args.outdir}/{target}.results"
                contents = os.popen(cmd).readlines()
                for line in contents:
                    tmp = line.rstrip('\n').split()
                    tmscores = tmp[1:]
                    break
                tmscores_rounded = []
                for tmscore in tmscores:
                    tmscore = str(round(float(tmscore),3))
                    tmscores_rounded += [tmscore]
                write_contents += tmscores_rounded
                print(write_contents)
                fw.write('\t'.join(write_contents) + '\n')                            



