import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
import pandas as pd
from bml_casp15.quaternary_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from bml_casp15.quaternary_structure_evaluation.pairwise_dockq import Pairwise_dockq_qa


class Quaternary_structure_evaluation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "pairwise"]):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods

        self.pairwise_qa = Pairwise_dockq_qa(params['dockq_program'])
        self.alphafold_qa = Alphafold_pkl_qa()

    def process(self, model_dir, output_dir):

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb'
        makedir_if_not_exists(pdbdir)

        pkldir = output_dir + '/pkl'
        makedir_if_not_exists(pkldir)

        for method in os.listdir(model_dir):
            for i in range(0, 5):
                os.system(f"cp {model_dir}/{method}/ranked_{i}.pdb {pdbdir}/{method}_{i}.pdb")
                #os.system(f"cp {model_dir}/{method}/result_model_{i + 1}_multimer.pkl {pkldir}/{method}_{i}.pkl")

        #if not os.path.exists(output_dir + '/pairwise_ranking.csv'):
        #    pairwise_ranking = self.pairwise_qa.run(pdbdir)
        #    pairwise_ranking.to_csv(output_dir + '/pairwise_ranking.csv')

        #if not os.path.exists(output_dir + '/alphafold_ranking.csv'):
        #    alphafold_ranking = self.alphafold_qa.run(pkldir)
        #    alphafold_ranking.to_csv(output_dir + '/alphafold_ranking.csv')
