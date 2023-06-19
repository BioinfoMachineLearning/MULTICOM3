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

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', type=str, required=True)
    parser.add_argument('--prefix', type=str, required=True)
    parser.add_argument('--outfile', type=str, required=True)
    args = parser.parse_args()

    infiles = args.infiles.split(',')

    models = []
    for infile in infiles:
        df = pd.read_csv(infile)
        for model in df['model']:
            modelname = model[model.find('TS'):]
            if modelname not in models:
                models += [modelname]
    
    data_dict = {'model': [], 'qscore': [], 'dockq': [], 'tmscore': []}

    for model in models:
        cmd = f"grep {model} {args.prefix}*"
        print(cmd)
        contents = os.popen(cmd).readlines()

        qscores, dockqs, tmscores = [], [], []
        for line in contents:
            tmp = line.split(',')
            qscore = float(tmp[2])
            dockq = float(tmp[3])
            tmscore = float(tmp[4])
            qscores += [qscore]
            dockqs += [dockq]
            tmscores += [tmscore]

        data_dict['model'] += [f"{args.prefix}{model}"]
        data_dict['qscore'] += [str(np.max(np.array(qscores)))]
        data_dict['dockq'] += [str(np.max(np.array(dockqs)))]
        data_dict['tmscore'] += [str(np.max(np.array(tmscores)))]

    pd.DataFrame(data_dict).to_csv(args.outfile)

