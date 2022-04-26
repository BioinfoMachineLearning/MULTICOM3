import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from bml_casp15.common.util import makedir_if_not_exists
import glob
import pathlib

class MultiEva_qa:

    def __init__(self, multieva_program):
        self.multieva_program = multieva_program

    def run(self, chain_id_map, fasta_path, input_dir, stoichiometry,
            alphafold_prediction_dir, outputdir):

        makedir_if_not_exists(outputdir)

        cmd = f"python {self.multieva_program} {fasta_path} {input_dir} {stoichiometry} {alphafold_prediction_dir} 20 {outputdir}"
            
        print(cmd)
        try:
            os.system(cmd)
        except Exception as e:
            print(e)

        result_csv = f"{outputdir}/{pathlib.Path(fasta_path).stem}.csv"
        print(result_csv)
        if not os.path.exists(result_csv):
            return None

        return result_csv

    def run_with_pairwise_ranking(self, input_dir, pkl_dir, pairwise_ranking_file, outputdir):

        input_dir = os.path.abspath(input_dir)
        pkl_dir = os.path.abspath(pkl_dir)
        outputdir = os.path.abspath(outputdir)

        makedir_if_not_exists(outputdir)

        pairwise_ranking = convert_pairwise_ranking_to_df(pairwise_ranking_file)

        makedir_if_not_exists(outputdir + '/alphafold_prediction')

        model_count = 5

        os.chdir(f"{outputdir}/alphafold_prediction")
        for i in range(len(pairwise_ranking)):
            if i >= model_count:
                break
            model = pairwise_ranking.loc[i, 'model']
            os.system(f"ln -s {input_dir}/{model} relaxed_model_{i+1}.pdb")
            os.system(f"ln -s {pkl_dir}/{model.replace('.pdb', '.pkl')} result_model_{i+1}.pkl")

        os.chdir(self.program_path)

        print("Try to run the ensemble model first")
        cmd = f"{self.enqa_program} {input_dir} {outputdir} ensemble {outputdir}/alphafold_prediction"
        if not self.use_gpu:
            cmd += " --cpu"

        print(cmd)
        try:
            os.system(cmd)
        except Exception as e:
            print(e)

        if not os.path.exists(f"{outputdir}/result.txt"):
            print("Try to run the EGNN_esto9 model")
            cmd = f"{self.enqa_program} {input_dir} {outputdir} EGNN_esto9 {outputdir}/alphafold_prediction"
            if not self.use_gpu:
                cmd += " --cpu"
            print(cmd)
            try:
                os.system(cmd)
            except Exception as e:
                print(e)

        models = []
        scores = []
        if os.path.exists(f"{outputdir}/result.txt"):
            for line in open(f"{outputdir}/result.txt"):
                line = line.rstrip('\n')
                contents = line.split()
                if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[
                    0] == "QMODE":
                    continue
                model, score = line.split()
                models += [model]
                scores += [score]
        else:
            pdbs = os.listdir(input_dir)
            for i in range(len(pdbs)):
                models += [pdbs[i]]
                scores += [0.0]

        return pd.DataFrame({'model': models, 'score': scores})
