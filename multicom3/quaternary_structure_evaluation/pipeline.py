import os, sys, argparse, time, json
from multiprocessing import Pool
from tqdm import tqdm
from multicom3.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
import pandas as pd
from multicom3.monomer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from multicom3.quaternary_structure_evaluation.pairwise_dockq import Pairwise_dockq_qa
from multicom3.quaternary_structure_evaluation.pairwise_mmalign import Pairwise_MMalign_qa
from multicom3.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa
from multicom3.common.protein import complete_result


class Quaternary_structure_evaluation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "bfactor", 'multieva']):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods

        self.multieva_qa = Pairwise_MMalign_qa(params['mmalign_program'])
        self.alphafold_qa = Alphafold_pkl_qa(sort_field='confidence')
        self.bfactorqa = Bfactor_qa()

    def process(self, fasta_path, chain_id_map, model_dir, output_dir, monomer_model_dir="", stoichiometry="", model_count=5):

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb/'
        makedir_if_not_exists(pdbdir)

        pkldir = output_dir + '/pkl/'
        makedir_if_not_exists(pkldir)

        msadir = output_dir + '/msa/'
        makedir_if_not_exists(msadir)

        for method in os.listdir(model_dir):
            print(method)
            ranking_json_file = f"{model_dir}/{method}/ranking_debug.json"
            if not os.path.exists(ranking_json_file):
                continue
            ranking_json = json.loads(open(ranking_json_file).read())

            for i in range(model_count):
                # if not complete_result(model_dir + '/' + method):
                #     continue
                os.system(f"cp {model_dir}/{method}/ranked_{i}.pdb {pdbdir}/{method}_{i}.pdb")

                model_name = list(ranking_json["order"])[i]
                os.system(f"cp {model_dir}/{method}/result_{model_name}.pkl {pkldir}/{method}_{i}.pkl")
                for chain_id in chain_id_map:
                    msa_chain_outdir = msadir + '/' + chain_id_map[chain_id].description
                    makedir_if_not_exists(msa_chain_outdir)
                    os.system(f"cp {model_dir}/{method}/msas/{chain_id}/monomer_final.a3m "
                              f"{msa_chain_outdir}/{method}_{i}.monomer.a3m")
                    os.system(f"cp {model_dir}/{method}/msas/{chain_id_map[chain_id].description}.paired.a3m "
                              f"{msa_chain_outdir}/{method}_{i}.paired.a3m")

        result_dict = {}

        if "alphafold" in self.run_methods:
            if not os.path.exists(output_dir + '/alphafold_ranking.csv'):
                alphafold_ranking = self.alphafold_qa.run(pkldir)
                alphafold_ranking.to_csv(output_dir + '/alphafold_ranking.csv')
            result_dict["alphafold"] = output_dir + '/alphafold_ranking.csv'

        if "bfactor" in self.run_methods:
            if not os.path.exists(output_dir + '/bfactor_ranking.csv'):
                bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
                bfactor_ranking.to_csv(output_dir + '/bfactor_ranking.csv')
            result_dict["bfactor"] = output_dir + '/bfactor_ranking.csv'

        if "multieva" in self.run_methods:
            if not os.path.exists(f"{output_dir}/multieva.csv"):
                workdir = output_dir + '/multieva'
                makedir_if_not_exists(workdir)
        
                multieva_pd = self.multieva.run(input_dir=pdbdir)
                multieva_pd.to_csv(f"{output_dir}/multieva.csv")

            result_dict["multieva"] = output_dir + '/multieva.csv'

        if "multieva" in self.run_methods and "alphafold" in self.run_methods:
            pairwise_ranking_df = pd.read_csv(result_dict["multieva"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
            pairwise_ranking_df['model'] = pairwise_ranking_df['Name'] + '.pdb'
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
                pairwise_score = float(avg_ranking_df.loc[i, 'MMalign score'])
                alphafold_score = float(avg_ranking_df.loc[i, 'confidence'])
                avg_score = (pairwise_score + alphafold_score) / 2
                avg_scores += [avg_score]
                avg_rank = (int(avg_ranking_df.loc[i, 'pairwise_rank']) + int(
                    avg_ranking_df.loc[i, 'alphafold_rank'])) / 2
                avg_rankings += [avg_rank]
            avg_ranking_df['avg_score'] = avg_scores
            avg_ranking_df['avg_rank'] = avg_rankings
            avg_ranking_df = avg_ranking_df.sort_values(by=['avg_score'], ascending=False)
            avg_ranking_df.reset_index(inplace=True, drop=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="index"), axis=1, inplace=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="Unnamed"), axis=1, inplace=True)
            avg_ranking_df.to_csv(output_dir + '/pairwise_af_avg.ranking')
            result_dict["pairwise_af_avg"] = output_dir + '/pairwise_af_avg.ranking'

        if "multieva" in self.run_methods and "bfactor" in self.run_methods:
            pairwise_ranking_df = pd.read_csv(result_dict["multieva"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
            pairwise_ranking_df['model'] = pairwise_ranking_df['Name'] + '.pdb'
            print(ranks)
            pairwise_ranking_df['pairwise_rank'] = ranks
            print(pairwise_ranking_df)
            bfactor_ranking_df = pd.read_csv(result_dict['bfactor'])
            ranks = [i + 1 for i in range(len(bfactor_ranking_df))]
            bfactor_ranking_df['bfactor_rank'] = ranks
            avg_ranking_df = pairwise_ranking_df.merge(bfactor_ranking_df, how="inner", on='model')
            avg_scores = []
            avg_rankings = []
            print(avg_ranking_df)
            for i in range(len(avg_ranking_df)):
                pairwise_score = float(avg_ranking_df.loc[i, 'MMalign score'])
                alphafold_score = float(avg_ranking_df.loc[i, 'bfactor'])/100
                avg_score = (pairwise_score + alphafold_score) / 2
                avg_scores += [avg_score]
                avg_rank = (int(avg_ranking_df.loc[i, 'pairwise_rank']) + int(
                    avg_ranking_df.loc[i, 'bfactor_rank'])) / 2
                avg_rankings += [avg_rank]
            avg_ranking_df['avg_score'] = avg_scores
            avg_ranking_df['avg_rank'] = avg_rankings
            avg_ranking_df = avg_ranking_df.sort_values(by=['avg_score'], ascending=False)
            avg_ranking_df.reset_index(inplace=True, drop=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="index"), axis=1, inplace=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="Unnamed"), axis=1, inplace=True)
            avg_ranking_df.to_csv(output_dir + '/pairwise_bfactor_avg.ranking')
            result_dict["pairwise_bfactor_avg"] = output_dir + '/pairwise_bfactor_avg.ranking'

        return result_dict

    def reprocess(self, fasta_path, chain_id_map, output_dir, monomer_model_dir="", stoichiometry=""):

        self.alphafold_qa = Alphafold_pkl_qa(sort_field='plddt_avg', ranking_methods=['plddt_avg'])

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb/'
        makedir_if_not_exists(pdbdir)

        pkldir = output_dir + '/pkl/'
        makedir_if_not_exists(pkldir)

        msadir = output_dir + '/msa/'
        makedir_if_not_exists(msadir)

        result_dict = {}

        if "alphafold" in self.run_methods:
            if not os.path.exists(output_dir + '/alphafold_ranking.csv'):
                alphafold_ranking = self.alphafold_qa.run(pkldir)
                alphafold_ranking.to_csv(output_dir + '/alphafold_ranking.csv')
            result_dict["alphafold"] = output_dir + '/alphafold_ranking.csv'

        if "bfactor" in self.run_methods:
            if not os.path.exists(output_dir + '/bfactor_ranking.csv'):
                bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
                bfactor_ranking.to_csv(output_dir + '/bfactor_ranking.csv')
            result_dict["bfactor"] = output_dir + '/bfactor_ranking.csv'

        if "multieva" in self.run_methods:
            if not os.path.exists(f"{output_dir}/multieva.csv"):
                workdir = output_dir + '/multieva'
                makedir_if_not_exists(workdir)
        
                multieva_pd = self.multieva.run(input_dir=pdbdir)
                multieva_pd.to_csv(f"{output_dir}/multieva.csv")
            result_dict["multieva"] = output_dir + '/multieva.csv'

        if "multieva" in self.run_methods and "alphafold" in self.run_methods:
            pairwise_ranking_df = pd.read_csv(result_dict["multieva"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
            pairwise_ranking_df['model'] = pairwise_ranking_df['Name'] + '.pdb'
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
                pairwise_score = float(avg_ranking_df.loc[i, 'MMalign score'])
                alphafold_score = float(avg_ranking_df.loc[i, 'plddt_avg'])/100
                avg_score = (pairwise_score + alphafold_score) / 2
                avg_scores += [avg_score]
                avg_rank = (int(avg_ranking_df.loc[i, 'pairwise_rank']) + int(
                    avg_ranking_df.loc[i, 'alphafold_rank'])) / 2
                avg_rankings += [avg_rank]
            avg_ranking_df['avg_score'] = avg_scores
            avg_ranking_df['avg_rank'] = avg_rankings
            avg_ranking_df = avg_ranking_df.sort_values(by=['avg_score'], ascending=False)
            avg_ranking_df.reset_index(inplace=True, drop=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="index"), axis=1, inplace=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="Unnamed"), axis=1, inplace=True)
            avg_ranking_df.to_csv(output_dir + '/pairwise_af_avg.ranking')
            result_dict["pairwise_af_avg"] = output_dir + '/pairwise_af_avg.ranking'

        if "multieva" in self.run_methods and "bfactor" in self.run_methods:
            pairwise_ranking_df = pd.read_csv(result_dict["multieva"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
            pairwise_ranking_df['model'] = pairwise_ranking_df['Name'] + '.pdb'
            print(ranks)
            pairwise_ranking_df['pairwise_rank'] = ranks
            print(pairwise_ranking_df)
            bfactor_ranking_df = pd.read_csv(result_dict['bfactor'])
            ranks = [i + 1 for i in range(len(bfactor_ranking_df))]
            bfactor_ranking_df['bfactor_rank'] = ranks
            avg_ranking_df = pairwise_ranking_df.merge(bfactor_ranking_df, how="inner", on='model')
            avg_scores = []
            avg_rankings = []
            print(avg_ranking_df)
            for i in range(len(avg_ranking_df)):
                pairwise_score = float(avg_ranking_df.loc[i, 'MMalign score'])
                alphafold_score = float(avg_ranking_df.loc[i, 'bfactor'])/100
                avg_score = (pairwise_score + alphafold_score) / 2
                avg_scores += [avg_score]
                avg_rank = (int(avg_ranking_df.loc[i, 'pairwise_rank']) + int(
                    avg_ranking_df.loc[i, 'bfactor_rank'])) / 2
                avg_rankings += [avg_rank]
            avg_ranking_df['avg_score'] = avg_scores
            avg_ranking_df['avg_rank'] = avg_rankings
            avg_ranking_df = avg_ranking_df.sort_values(by=['avg_score'], ascending=False)
            avg_ranking_df.reset_index(inplace=True, drop=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="index"), axis=1, inplace=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="Unnamed"), axis=1, inplace=True)
            avg_ranking_df.to_csv(output_dir + '/pairwise_bfactor_avg.ranking')
            result_dict["pairwise_bfactor_avg"] = output_dir + '/pairwise_bfactor_avg.ranking'

        return result_dict
