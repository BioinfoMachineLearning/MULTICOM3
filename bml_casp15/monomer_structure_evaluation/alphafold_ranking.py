import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle


class Alphafold_pkl_qa:

    def __init__(self, ranking_methods = ['ptm', 'plddt_avg', 'confidence']):

        self.methods = ranking_methods

    def run(self, input_dir):
        ranking_pd = pd.DataFrame(columns=['model'] + self.methods)
        model_count = 0
        for pkl in os.listdir(input_dir):
            if pkl.index('.pkl') < 0:
                continue
            ranking = {'model': pkl.replace('.pkl', '.pdb')}
            with open(input_dir + '/' + pkl, 'rb') as f:
                prediction_result = pickle.load(f)
                if 'plddt_avg' in self.methods:
                    ranking['plddt_avg'] = np.mean(prediction_result['plddt'])
                if 'ptm' in self.methods:
                    ranking['ptm'] = float(prediction_result['ptm'])
                if 'confidence' in self.methods:
                    ranking['confidence'] = float(prediction_result['ranking_confidence'])
                ranking_pd = ranking_pd.append(pd.DataFrame(ranking, index=[model_count]))
                model_count += 1
        return ranking_pd#.sort_values(by=['confidence'], ascending=False, ignore_index=True)
