import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
import pandas as pd
from bml_casp15.tertiary_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from bml_casp15.tertiary_structure_evaluation.enqa_ranking import En_qa
import copy
import pickle


def extract_monomer_pdb(complex_pdb, chainid, output_pdb):
    chain_start = -1
    chain_end = -1
    residue_count = 0
    prev_num = 0
    with open(output_pdb, 'w') as fw:
        for line in open(complex_pdb, 'r').readlines():
            if not line.startswith('ATOM'):
                continue

            residue_num = int(line[22:26].strip())
            if residue_num != prev_num:
                residue_count += 1
                prev_num = residue_num

            chain_name = line[21]
            if chain_name == chainid:
                # print(line)
                fw.write(line)
                if chain_start == -1:
                    chain_start = residue_count
            elif chain_start >= 0:
                chain_end = residue_count
                break

    if chain_start >= 0 and chain_end == -1:
        chain_end = residue_count + 1

    return chain_start - 1, chain_end - 1


def extract_monomer_pkl(complex_pkl, residue_start, residue_end, output_pkl):
    with open(complex_pkl, 'rb') as f:
        prediction_result = pickle.load(f)
        prediction_result_monomer = copy.deepcopy(prediction_result)
        prediction_result_monomer['plddt'] = prediction_result_monomer['plddt'][residue_start:residue_end, ]
        with open(output_pkl, 'wb') as f:
            pickle.dump(prediction_result_monomer, f, protocol=4)


class Tertiary_structure_evaluation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "apollo", "enQA"], use_gpu=True):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods
        self.alphafold_qa = Alphafold_pkl_qa(ranking_methods=['plddt_avg'])
        self.parwise_qa = params['qscore_program']
        self.tmscore = params['tmscore_program']
        self.enqa = En_qa(program_python=params['enqa_program_python'], program_path=params['enqa_program_path'],
                          program_script=params['enqa_program_script'], use_gpu=use_gpu)

    def process(self, targetname, fasta_file, monomer_model_dir, output_dir, multimer_model_dir="",
                chainid_in_multimer=""):

        output_dir_abs = os.path.abspath(output_dir)

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb'
        makedir_if_not_exists(pdbdir)

        pkldir = output_dir + '/pkl'
        makedir_if_not_exists(pkldir)

        for method in os.listdir(monomer_model_dir):
            for i in range(0, 5):
                os.system(f"cp {monomer_model_dir}/{method}/ranked_{i}.pdb {pdbdir}/{method}_{i}.pdb")
                os.system(f"cp {monomer_model_dir}/{method}/result_model_{i+1}.pkl {pkldir}/{method}_{i}.pkl")

        for method in os.listdir(multimer_model_dir):
            for i in range(0, 5):
                complex_pdb = f"{multimer_model_dir}/{method}/ranked_{i}.pdb"
                if os.path.exists(complex_pdb):
                    residue_start, residue_end = extract_monomer_pdb(
                        complex_pdb=complex_pdb,
                        chainid=chainid_in_multimer,
                        output_pdb=f"{pdbdir}/{method}_{i}.pdb")
                    extract_monomer_pkl(complex_pkl=f"{multimer_model_dir}/{method}/result_model_{i + 1}_multimer.pkl",
                                        residue_start=residue_start,
                                        residue_end=residue_end,
                                        output_pkl=f"{pkldir}/{method}_{i}.pkl")
                # print(f"start: {residue_start}, end:{residue_end}")

        if not os.path.exists(output_dir + '/pairwise_ranking.tm'):
            with open(f"{output_dir}/model.list", 'w') as fw:
                for pdb in os.listdir(pdbdir):
                    fw.write(f"{pdbdir}/{pdb}\n")
            os.system(
                f"{self.parwise_qa} {output_dir}/model.list {fasta_file} {self.tmscore} {output_dir} pairwise_ranking")

        if not os.path.exists(output_dir_abs + '/alphafold_ranking.csv'):
            alphafold_ranking = self.alphafold_qa.run(pkldir)
            alphafold_ranking.to_csv(output_dir_abs + '/alphafold_ranking.csv')

        if not os.path.exists(output_dir_abs + '/enqa_ranking.csv'):
            enqa_ranking = self.enqa.run(input_dir=pdbdir,
                                         alphafold_prediction_dir=f"{monomer_model_dir}/original",
                                         outputdir=output_dir_abs+'/enqa')
            enqa_ranking.to_csv(output_dir_abs + '/enqa_ranking.csv')
