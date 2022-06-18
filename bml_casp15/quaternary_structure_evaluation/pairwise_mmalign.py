import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd


def run_command(inparams):
    mmalign_program, input_dir, pdb1, pdb2 = inparams
    cmd = mmalign_program + ' ' + input_dir + '/' + pdb1 + ' ' + input_dir + '/' + pdb2 + " | grep TM-score | awk '{print $2}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[1].rstrip('\n'))
    return pdb1, pdb2, tmscore


class Pairwise_MMalign_qa:

    def __init__(self, mmalign_program):

        self.mmalign_program = mmalign_program

    def run(self, input_dir):

        ranking_pd = pd.DataFrame(columns=['model', 'pairwise_score'])

        pdbs = os.listdir(input_dir)

        process_list = []
        for i in range(len(pdbs)):
            for j in range(len(pdbs)):
                pdb1 = pdbs[i]
                pdb2 = pdbs[j]
                if pdb1 == pdb2:
                    continue
                process_list.append([self.mmalign_program, input_dir, pdb1, pdb2])

        pool = Pool(processes=40)
        results = pool.map(run_command, process_list)
        pool.close()
        pool.join()

        scores_dict = {}
        for result in results:
            pdb1, pdb2, score = result
            scores_dict[f"{pdb1}_{pdb2}"] = score

        for i in range(len(pdbs)):
            pdb1 = pdbs[i]
            ranking = {'model': pdb1}
            scores = []
            for pdb2 in pdbs:
                if pdb1 == pdb2:
                    continue
                if f"{pdb1}_{pdb2}" in scores_dict:
                    score = scores_dict[f"{pdb1}_{pdb2}"]
                scores += [score]

            ranking['pairwise_score'] = np.mean(np.array(scores))

            ranking_pd = ranking_pd.append(pd.DataFrame(ranking, index=[i]))

        return ranking_pd.sort_values(by=['pairwise_score'], ascending=False, ignore_index=True)
