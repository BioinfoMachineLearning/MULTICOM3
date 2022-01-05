import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd


class Pairwise_dockq_qa:

    def __init__(self, dockq_program):

        self.dockq_program = dockq_program

    def run(self, input_dir):

        ranking_pd = pd.DataFrame(columns=['model', 'pairwise_score'])

        pdbs = os.listdir(input_dir)

        scores_dict = {}

        for i in range(len(pdbs)):

            pdb1 = pdbs[i]

            ranking = {'model': pdb1}

            scores = []

            for pdb2 in pdbs:

                if pdb1 == pdb2:
                    continue

                if f"{pdb1}_{pdb2}" in scores_dict:

                    score = scores_dict[f"{pdb1}_{pdb2}"]

                elif f"{pdb2}_{pdb1}" in scores_dict:

                    score = scores_dict[f"{pdb2}_{pdb1}"]

                else:

                    cmd = f"python {self.dockq_program} -short " + input_dir + '/' + pdb1 + " " \
                          + input_dir + '/' + pdb2 + " | grep DockQ | awk '{print $2}'"

                    score = os.popen(cmd).read()

                    score = float(score.rstrip('\n'))

                    scores_dict[f"{pdb1}_{pdb2}"] = score

                scores += [score]

            ranking['pairwise_score'] = np.mean(np.array(scores))

            ranking_pd = ranking_pd.append(pd.DataFrame(ranking, index=[i]))

        return ranking_pd.sort_values(by=['pairwise_score'], ascending=False, ignore_index=True)
