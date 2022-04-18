import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
import pandas as pd
from bml_casp15.monomer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from bml_casp15.monomer_structure_evaluation.enqa_ranking import En_qa
import copy
import pickle


class Monomer_structure_evaluation_human_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "apollo", "enQA", "human"], use_gpu=True):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods
        self.alphafold_qa = Alphafold_pkl_qa(ranking_methods=['plddt_avg'])
        self.parwise_qa = params['qscore_program']
        self.tmscore = params['tmscore_program']
        self.enqa = En_qa(program_python=params['enqa_program_python'], program_path=params['enqa_program_path'],
                          program_script=params['enqa_program_script'], use_gpu=use_gpu)
        self.casp13_human = params['CASP13_human_script']
        self.casp14_human = params['CASP14_human_script']

    def process(self, targetname, fasta_file, monomer_model_dir, output_dir, multimer_model_dir="",
                chainid_in_multimer="", unrelaxed_chainid_in_multimer=""):

        output_dir_abs = os.path.abspath(output_dir)

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb'
        makedir_if_not_exists(pdbdir)

        pkldir = output_dir + '/pkl'
        makedir_if_not_exists(pkldir)

        msadir = output_dir + '/msa'
        makedir_if_not_exists(msadir)

        for method in os.listdir(monomer_model_dir):
            for i in range(0, 5):
                os.system(f"cp {monomer_model_dir}/{method}/ranked_{i}.pdb {pdbdir}/{method}_{i}.pdb")
                extract_pkl(src_pkl=f"{monomer_model_dir}/{method}/result_model_{i + 1}.pkl",
                            output_pkl=f"{pkldir}/{method}_{i}.pkl")
                os.system(f"cp {monomer_model_dir}/{method}/msas/monomer_final.a3m {msadir}/{method}_{i}.a3m")

        if os.path.exists(multimer_model_dir):
            for method in os.listdir(multimer_model_dir):
                for i in range(0, 5):
                    complex_pdb = f"{multimer_model_dir}/{method}/ranked_{i}.pdb"
                    if os.path.exists(complex_pdb):
                        residue_start, residue_end = extract_monomer_pdb(
                            complex_pdb=complex_pdb,
                            chainid=unrelaxed_chainid_in_multimer,
                            output_pdb=f"{pdbdir}/{method}_{i}.pdb")
                        extract_pkl(src_pkl=f"{multimer_model_dir}/{method}/result_model_{i + 1}_multimer.pkl",
                                    residue_start=residue_start,
                                    residue_end=residue_end,
                                    output_pkl=f"{pkldir}/{method}_{i}.pkl")
                        os.system(f"cp {multimer_model_dir}/{method}/msas/{chainid_in_multimer}/monomer_final.a3m "
                                  f"{msadir}/{method}_{i}.a3m")
                    # print(f"start: {residue_start}, end:{residue_end}")

        result_dict = {}
        cwd = os.getcwd()

        if "apollo" in self.run_methods:
            if not os.path.exists(output_dir + '/pairwise_ranking.tm'):
                os.chdir(output_dir)
                with open("model.list", 'w') as fw:
                    for pdb in os.listdir(pdbdir):
                        fw.write(f"pdb/{pdb}\n")
                os.system(
                    f"{self.parwise_qa} model.list {fasta_file} {self.tmscore} . pairwise_ranking")
            result_dict["apollo"] = output_dir + '/pairwise_ranking.tm'

        if "alphafold" in self.run_methods:
            if not os.path.exists(output_dir_abs + '/alphafold_ranking.csv'):
                alphafold_ranking = self.alphafold_qa.run(pkldir)
                alphafold_ranking.to_csv(output_dir_abs + '/alphafold_ranking.csv')
            result_dict["alphafold"] = output_dir_abs + '/alphafold_ranking.csv'

        if "enQA" in self.run_methods:
            if not os.path.exists(output_dir_abs + '/enqa_ranking.csv'):
                # enqa_ranking = self.enqa.run(input_dir=pdbdir,
                #                              alphafold_prediction_dir=f"{monomer_model_dir}/original",
                #                              outputdir=output_dir_abs+'/enqa')
                enqa_ranking = self.enqa.run_with_pairwise_ranking(input_dir=pdbdir,
                                                                   pkl_dir=pkldir,
                                                                   pairwise_ranking_file=output_dir + '/pairwise_ranking.tm',
                                                                   outputdir=output_dir_abs + '/enqa')
                enqa_ranking.to_csv(output_dir_abs + '/enqa_ranking.csv')
            result_dict["enQA"] = output_dir_abs + '/enqa_ranking.csv'

        if "human" in self.run_methods:
            casp13_ranking = output_dir_abs + '/casp13/Full_TS/eva/Human_TS.txt'
            casp14_deeprank_ranking = output_dir_abs + '/casp14/Full_TS/eva/Human_TS.txt'
            casp14_deeprank3_ranking = output_dir_abs + '/casp14/Full_TS/eva_DeepRank3/DeepRank3_Cluster.txt'
            if not os.path.exists(casp13_ranking):
                cmd = f"{self.casp13_script} {targetname} {fasta_file} {output_dir_abs}/casp13 NULL NULL {pdbdir} &> {output_dir_abs}/run_casp13_human.log "    
                os.system(cmd)
            if not os.path.exists(casp14_deeprank_ranking) or not os.path.exists(casp14_deeprank3_ranking):
                cmd = f"{self.casp14_script} {targetname} {fasta_file} {output_dir_abs}/casp14 NULL NULL NULL {pdbdir} &> {output_dir_abs}/run_casp14_human.log "
                os.system(cmd)
            result_dict['casp13'] = casp13_ranking
            result_dict['casp14_deeprank'] = casp14_deeprank_ranking
            result_dict['casp14_deeprank3'] = casp14_deeprank3_ranking

        os.chdir(cwd)

        return result_dict
