import os, sys, argparse, time, json
from multiprocessing import Pool
from tqdm import tqdm
from multicom3.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
import pandas as pd
from multicom3.monomer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from multicom3.multimer_structure_evaluation.pairwise_mmalign import Pairwise_MMalign_qa
from multicom3.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa
from multicom3.common.protein import complete_result


class Multimer_structure_evaluation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "bfactor", 'multieva']):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods

        self.multieva_qa = Pairwise_MMalign_qa(params['mmalign_program'])
        self.alphafold_qa = Alphafold_pkl_qa(sort_field='confidence')
        self.bfactorqa = Bfactor_qa()

    def process(self, fasta_path, chain_id_map, model_dir, output_dir, is_homomer=False, model_count=5):

        makedir_if_not_exists(output_dir)

        pdbdir = os.path.join(output_dir, 'pdb')
        makedir_if_not_exists(pdbdir)

        pkldir = os.path.join(output_dir, 'pkl')
        makedir_if_not_exists(pkldir)

        msadir = os.path.join(output_dir, 'msa')
        makedir_if_not_exists(msadir)

        for method in os.listdir(model_dir):
            print(method)
            ranking_json_file = os.path.join(model_dir, method, "ranking_debug.json")
            if not os.path.exists(ranking_json_file):
                continue
            ranking_json = json.loads(open(ranking_json_file).read())

            for i in range(model_count):
                ranked_pdb = os.path.join(model_dir, method, f"ranked_{i}.pdb")
                trg_pdb = os.path.join(pdbdir, f"{method}_{i}.pdb")
                os.system(f"cp {ranked_pdb} {trg_pdb}")

                model_name = list(ranking_json["order"])[i]
                ranked_pkl = os.path.join(model_dir, method, f"result_{model_name}.pkl")
                trg_pkl = os.path.join(pkldir, f"{method}_{i}.pkl")
                os.system(f"cp {ranked_pkl} {trg_pkl}")
                for chain_id in chain_id_map:
                    msa_chain_outdir = os.path.join(msadir, chain_id)
                    makedir_if_not_exists(msa_chain_outdir)
                    src_monomer_a3m = os.path.join(model_dir, method, 'msas', chain_id, "monomer_final.a3m")
                    trg_monomer_a3m = os.path.join(msa_chain_outdir, f"{method}_{i}.monomer.a3m")
                    os.system(f"cp {src_monomer_a3m} {trg_monomer_a3m}")
                    if not is_homomer:
                        src_paired_a3m = os.path.join(model_dir, method, 'msas', chain_id + ".paired.a3m")
                        trg_paried_a3m = os.path.join(msa_chain_outdir, f"{method}_{i}.paired.a3m")
                        os.system(f"cp {src_paired_a3m} {trg_paried_a3m}")

        result_dict = {}

        if "alphafold" in self.run_methods:
            result_file = os.path.join(output_dir, 'alphafold_ranking.csv')
            if not os.path.exists(result_file):
                alphafold_ranking = self.alphafold_qa.run(pkldir)
                alphafold_ranking.to_csv(result_file)
            result_dict["alphafold"] = result_file

        if "bfactor" in self.run_methods:
            result_file = os.path.join(output_dir, 'bfactor_ranking.csv')
            if not os.path.exists(result_file):
                bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
                bfactor_ranking.to_csv(result_file)
            result_dict["bfactor"] = result_file

        if "multieva" in self.run_methods:
            result_file = os.path.join(output_dir, 'multieva.csv')
            if not os.path.exists(result_file):
                workdir = os.path.join(output_dir, 'multieva')
                makedir_if_not_exists(workdir)
        
                multieva_pd = self.multieva_qa.run(input_dir=pdbdir)
                multieva_pd.to_csv(result_file)

            result_dict["multieva"] = result_file

        if "multieva" in self.run_methods and "alphafold" in self.run_methods:
            pairwise_ranking_df = pd.read_csv(result_dict["multieva"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
            pairwise_ranking_df['model'] = pairwise_ranking_df['Name']
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

            result_file = os.path.join(output_dir, 'pairwise_af_avg.ranking')
            avg_ranking_df.to_csv(result_file)
            result_dict["pairwise_af_avg"] = result_file

        if "multieva" in self.run_methods and "bfactor" in self.run_methods:
            pairwise_ranking_df = pd.read_csv(result_dict["multieva"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
            pairwise_ranking_df['model'] = pairwise_ranking_df['Name']
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

            result_file = os.path.join(output_dir, 'pairwise_bfactor_avg.ranking')
            avg_ranking_df.to_csv(result_file)
            result_dict["pairwise_bfactor_avg"] = result_file

        return result_dict

    def reprocess(self, fasta_path, chain_id_map, output_dir):

        self.alphafold_qa = Alphafold_pkl_qa(sort_field='plddt_avg', ranking_methods=['plddt_avg'])

        makedir_if_not_exists(output_dir)

        pdbdir = os.path.join(output_dir, 'pdb')
        makedir_if_not_exists(pdbdir)

        pkldir = os.path.join(output_dir, 'pkl')
        makedir_if_not_exists(pkldir)

        msadir = os.path.join(output_dir, 'msa')
        makedir_if_not_exists(msadir)

        result_dict = {}

        if "alphafold" in self.run_methods:
            result_file = os.path.join(output_dir, 'alphafold_ranking.csv')
            if not os.path.exists(result_file):
                alphafold_ranking = self.alphafold_qa.run(pkldir)
                alphafold_ranking.to_csv(result_file)
            result_dict["alphafold"] = result_file

        if "bfactor" in self.run_methods:
            result_file = os.path.join(output_dir, 'bfactor_ranking.csv')
            if not os.path.exists(result_file):
                bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
                bfactor_ranking.to_csv(result_file)
            result_dict["bfactor"] = result_file

        if "multieva" in self.run_methods:
            result_file = os.path.join(output_dir, 'multieva.csv')
            if not os.path.exists(result_file):
                workdir = os.path.join(output_dir, 'multieva')
                makedir_if_not_exists(workdir)
        
                multieva_pd = self.multieva_qa.run(input_dir=pdbdir)
                multieva_pd.to_csv(result_file)

            result_dict["multieva"] = result_file

        if "multieva" in self.run_methods and "alphafold" in self.run_methods:
            pairwise_ranking_df = pd.read_csv(result_dict["multieva"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
            pairwise_ranking_df['model'] = pairwise_ranking_df['Name']
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

            result_file = os.path.join(output_dir, 'pairwise_af_avg.ranking')
            avg_ranking_df.to_csv(result_file)
            result_dict["pairwise_af_avg"] = result_file

        if "multieva" in self.run_methods and "bfactor" in self.run_methods:
            pairwise_ranking_df = pd.read_csv(result_dict["multieva"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
            pairwise_ranking_df['model'] = pairwise_ranking_df['Name']
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

            result_file = os.path.join(output_dir, 'pairwise_bfactor_avg.ranking')
            avg_ranking_df.to_csv(result_file)
            result_dict["pairwise_bfactor_avg"] = result_file

        return result_dict
