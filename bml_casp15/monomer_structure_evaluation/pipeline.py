import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
import pandas as pd
from bml_casp15.monomer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from bml_casp15.monomer_structure_evaluation.enqa_ranking import En_qa
from bml_casp15.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa
from bml_casp15.common.protein import read_qa_txt_as_df
import copy
import pickle
import json


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


def extract_pkl(src_pkl, output_pkl, residue_start=-1, residue_end=-1):
    with open(src_pkl, 'rb') as f:
        prediction_result = pickle.load(f)
        if residue_end != -1 and residue_start != -1:
            prediction_result_monomer = {'plddt': prediction_result['plddt'][residue_start:residue_end, ]}
            distogram_monomer = prediction_result['distogram']
            distogram_monomer['logits'] = distogram_monomer['logits'][residue_start:residue_end,
                                          residue_start:residue_end, :]
            prediction_result_monomer['distogram'] = distogram_monomer

        else:
            prediction_result_monomer = {'plddt': prediction_result['plddt'],
                                         'distogram': prediction_result['distogram']}
        with open(output_pkl, 'wb') as f:
            pickle.dump(prediction_result_monomer, f, protocol=4)


class Monomer_structure_evaluation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "apollo", "enQA", "bfactor"], use_gpu=True):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods
        self.alphafold_qa = Alphafold_pkl_qa(ranking_methods=['plddt_avg'])
        self.parwise_qa = params['qscore_program']
        self.tmscore = params['tmscore_program']
        self.enqa = En_qa(enqa_program=params['enqa_program'], use_gpu=use_gpu)
        self.bfactorqa = Bfactor_qa()

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
            ranking_json_file = f"{monomer_model_dir}/{method}/ranking_debug.json"
            ranking_json = json.loads(open(ranking_json_file).read())
            for i in range(0, 5):
                os.system(f"cp {monomer_model_dir}/{method}/ranked_{i}.pdb {pdbdir}/{method}_{i}.pdb")

                model_name = list(ranking_json["order"])[i]
                extract_pkl(src_pkl=f"{monomer_model_dir}/{method}/result_{model_name}.pkl",
                            output_pkl=f"{pkldir}/{method}_{i}.pkl")
                os.system(f"cp {monomer_model_dir}/{method}/msas/monomer_final.a3m {msadir}/{method}_{i}.a3m")

        if os.path.exists(multimer_model_dir):
            for method in os.listdir(multimer_model_dir):
                ranking_json_file = f"{multimer_model_dir}/{method}/ranking_debug.json"
                if not os.path.exists(ranking_json_file):
                    continue
                ranking_json = json.loads(open(ranking_json_file).read())
                for i in range(0, 5):
                    complex_pdb = f"{multimer_model_dir}/{method}/ranked_{i}.pdb"
                    if os.path.exists(complex_pdb):
                        residue_start, residue_end = extract_monomer_pdb(
                            complex_pdb=complex_pdb,
                            chainid=unrelaxed_chainid_in_multimer,
                            output_pdb=f"{pdbdir}/{method}_{i}.pdb")

                        model_name = list(ranking_json["order"])[i]
                        extract_pkl(src_pkl=f"{multimer_model_dir}/{method}/result_{model_name}.pkl",
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
                print(f"{self.parwise_qa} model.list {fasta_file} {self.tmscore} . pairwise_ranking")
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
        if "bfactor" in self.run_methods:
            if not os.path.exists(output_dir + '/bfactor_ranking.csv'):
                bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
                bfactor_ranking.to_csv(output_dir + '/bfactor_ranking.csv')
            result_dict["bfactor"] = output_dir + '/bfactor_ranking.csv'

        if "apollo" in self.run_methods and "alphafold" in self.run_methods:
            pairwise_ranking_df = read_qa_txt_as_df(result_dict["apollo"])
            ranks = [i+1 for i in range(len(pairwise_ranking_df))]
            print(ranks)
            pairwise_ranking_df['pairwise_rank'] = ranks
            print(pairwise_ranking_df)
            alphafold_ranking_df = pd.read_csv(result_dict['alphafold'])
            ranks = [i + 1 for i in range(len(alphafold_ranking_df))]
            alphafold_ranking_df['alphafold_rank'] = ranks
            avg_ranking_df = pairwise_ranking_df.merge(alphafold_ranking_df, how="inner", on='model')
            avg_scores = []
            avg_rankings = []
            print(avg_ranking_df)
            for i in range(len(avg_ranking_df)):
                pairwise_score = float(avg_ranking_df.loc[i, 'score'])
                alphafold_score = float(avg_ranking_df.loc[i, 'plddt_avg']) / 100
                avg_score = (pairwise_score + alphafold_score) / 2
                avg_scores += [avg_score]
                avg_rank = (int(avg_ranking_df.loc[i, 'pairwise_rank'])+int(avg_ranking_df.loc[i, 'alphafold_rank']))/2
                avg_rankings += [avg_rank]
            avg_ranking_df['avg_score'] = avg_scores
            avg_ranking_df['avg_rank'] = avg_rankings
            avg_ranking_df = avg_ranking_df.sort_values(by=['avg_score'], ascending=False)
            avg_ranking_df.reset_index(inplace=True, drop=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="index"), axis=1, inplace=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="Unnamed"), axis=1, inplace=True)
            avg_ranking_df.to_csv(output_dir + '/pairwise_af_avg.ranking')
            result_dict["pairwise_af_avg"] = output_dir + '/pairwise_af_avg.ranking'

        os.chdir(cwd)

        return result_dict

    def reprocess(self, targetname, fasta_file, output_dir):

        output_dir_abs = os.path.abspath(output_dir)

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb'
        pkldir = output_dir + '/pkl'
        msadir = output_dir + '/msa'

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

        if "bfactor" in self.run_methods:
            if not os.path.exists(output_dir + '/bfactor_ranking.csv'):
                bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
                bfactor_ranking.to_csv(output_dir + '/bfactor_ranking.csv')
            result_dict["bfactor"] = output_dir + '/bfactor_ranking.csv'

        if "apollo" in self.run_methods and "alphafold" in self.run_methods:
            pairwise_ranking_df = read_qa_txt_as_df(result_dict["apollo"])
            alphafold_ranking_df = pd.read_csv(result_dict['alphafold'])
            avg_ranking_df = pairwise_ranking_df.merge(alphafold_ranking_df, how="inner", on='model')
            avg_scores = []
            for i in range(len(avg_ranking_df)):
                pairwise_score = float(avg_ranking_df.loc[i, 'score'])
                alphafold_score = float(avg_ranking_df.loc[i, 'plddt_avg']) / 100
                avg_score = (pairwise_score + alphafold_score) / 2
                avg_scores += [avg_score]
            avg_ranking_df['avg_score'] = avg_scores
            avg_ranking_df = avg_ranking_df.sort_values(by=['avg_score'], ascending=False)
            avg_ranking_df.reset_index(inplace=True, drop=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="index"), axis=1, inplace=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="Unnamed"), axis=1, inplace=True)
            avg_ranking_df.to_csv(output_dir + '/pairwise_af_avg.ranking')

        os.chdir(cwd)

        return result_dict
