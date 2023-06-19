import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir
from bml_casp15.quaternary_structure_evaluation.pairwise_mmalign import *

def run_command(inparams):
    mmalign_program, pdb1, pdb2, outfile = inparams
    cmd = f"{mmalign_program}  {pdb1}  {pdb2} > {outfile} "
    os.system(cmd)

def read_mmalign(infile):
    # print(infile)
    for line in open(infile):
        line = line.rstrip('\n')
        if len(line) == 0:
            continue
        if line.split()[0] == 'TM-score=':
            # print(line)
            if line.find('Structure_2') >= 0:
                return float(line.split()[1])
    return 0


def run_pairwise(mmalign_program, indir, tempdir, outfile, outfile2):

    ranking_pd = pd.DataFrame(columns=['model', 'score'])

    pdbs = os.listdir(indir)

    process_list = []
    for i in range(len(pdbs)):
        for j in range(len(pdbs)):
            pdb1 = pdbs[i]
            pdb2 = pdbs[j]
            if pdb1 == pdb2:
                continue
            tmscore_file = f"{tempdir}/{pdb1}_{pdb2}.mmalign"
            if not os.path.exists(tmscore_file):
                process_list.append([mmalign_program, indir + '/' + pdb1, indir + '/' + pdb2, tmscore_file])

    pool = Pool(processes=40)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()

    for i in range(len(pdbs)):
        pdb1 = pdbs[i]
        ranking = {'model': pdb1}
        scores = []
        for pdb2 in pdbs:
            if pdb1 == pdb2:
                continue

            tmscore_file = f"{tempdir}/{pdb1}_{pdb2}.mmalign"
            score = read_mmalign(tmscore_file)
            scores += [score]

        ranking['score'] = np.mean(np.array(scores))

        ranking_pd = ranking_pd.append(pd.DataFrame(ranking, index=[i]))

    ranking_pd.sort_values(by=['score'], ascending=False, ignore_index=True).to_csv(outfile)

    with open(outfile2, 'w') as fw:
        for model, score in zip(ranking_pd['model'], ranking_pd['score']):
            fw.write(f"{model} {score}\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True)
    parser.add_argument('--mmalign_program', type=is_file, required=True)

    args = parser.parse_args()

    for target in os.listdir(args.indir):
        workdir = args.indir + '/' + target
        tempdir = workdir + '/temp'
        if not os.path.exists(tempdir):
            os.makedirs(tempdir)

        outfile = workdir + '/' + target + '.csv'

        run_pairwise(args.mmalign_program, workdir + '/pdb', tempdir, outfile, f"{args.indir}/{target}/{target}QA2_1")
