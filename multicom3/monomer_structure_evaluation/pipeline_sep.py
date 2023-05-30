import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from multicom3.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, clean_dir
import pandas as pd
from multicom3.monomer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from multicom3.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa
from multicom3.common.protein import read_qa_txt_as_df
import copy
import pickle
import json
import pathlib


def get_sequence(contents):
    """Enclosing logic in a function to simplify code"""

    res_codes = [
        # 20 canonical amino acids
        ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
        ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
        ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
        ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
        ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
        # Non-canonical amino acids
        # ('MSE', 'M'), ('SOC', 'C'),
        # Canonical xNA
        ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
        ('  T', 'T'),
    ]

    three_to_one = dict(res_codes)
    # _records = set(['ATOM  ', 'HETATM'])
    _records = set(['ATOM  '])

    sequence = []
    read = set()
    for line in contents:
        line = line.strip()
        if line[0:6] in _records:
            resn = line[17:20]
            resi = line[22:26]
            icode = line[26]
            r_uid = (resn, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
            else:
                continue
            aa_resn = three_to_one.get(resn, 'X')
            sequence.append(aa_resn)

    return ''.join(sequence)


def extract_monomer_pdbs(complex_pdb, sequence, output_prefix):
    chain_contents = {}
    prev_chain_name = ""
    residue_count = 0
    prev_num = 0
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        residue_num = int(line[22:26].strip())
        if residue_num != prev_num:
            residue_count += 1
            prev_num = residue_num

        chain_id = line[21]
        if chain_id not in chain_contents:
            if chain_id != prev_chain_name and prev_chain_name != "":
                chain_contents[prev_chain_name]['res_end'] = residue_count - 1
                chain_contents[prev_chain_name]['sequence'] = get_sequence(chain_contents[prev_chain_name]['content'])
            chain_contents[chain_id] = dict(content=[line],
                                            sequence="",
                                            res_start=residue_count,
                                            res_end=-1)
            prev_chain_name = chain_id
        else:
            chain_contents[chain_id]['content'] += [line]

    if prev_chain_name != "":
        chain_contents[prev_chain_name]['res_end'] = residue_count
        chain_contents[prev_chain_name]['sequence'] = get_sequence(chain_contents[prev_chain_name]['content'])

    result_dict = {}
    for chain_id in chain_contents:
        if chain_contents[chain_id]['sequence'] == sequence:
            with open(f"{output_prefix}_{chain_id}.pdb", 'w') as fw:
                fw.writelines(chain_contents[chain_id]['content'])
            result_dict[chain_id] = dict(pdb=f"{output_prefix}_{chain_id}.pdb",
                                         pdbname=pathlib.Path(f"{output_prefix}_{chain_id}.pdb").name,
                                         chain_start=chain_contents[chain_id]['res_start'] - 1,
                                         chain_end=chain_contents[chain_id]['res_end'] - 1)
    return result_dict


def extract_pkl(src_pkl, output_pkl, residue_start=-1, residue_end=-1):
    has_distogram = False
    with open(src_pkl, 'rb') as f:
        prediction_result = pickle.load(f)
        if residue_end != -1 and residue_start != -1:
            prediction_result_monomer = {'plddt': prediction_result['plddt'][residue_start:residue_end + 1, ]}
            if 'distogram' in prediction_result:
                distogram_monomer = prediction_result['distogram']
                distogram_monomer['logits'] = distogram_monomer['logits'][residue_start:residue_end + 1,
                                              residue_start:residue_end + 1, :]
                prediction_result_monomer['distogram'] = distogram_monomer
                has_distogram = True
        else:
            prediction_result_monomer = {'plddt': prediction_result['plddt']}
            if 'distogram' in prediction_result:
                prediction_result_monomer['distogram'] = prediction_result['distogram']
                has_distogram = True
        with open(output_pkl, 'wb') as f:
            pickle.dump(prediction_result_monomer, f, protocol=4)
    return has_distogram


class Monomer_structure_evaluation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "apollo", "bfactor"], use_gpu=True):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods
        self.alphafold_qa = Alphafold_pkl_qa(ranking_methods=['plddt_avg'])
        self.parwise_qa = params['qscore_program']
        self.tmscore = params['tmscore_program']
        self.bfactorqa = Bfactor_qa()

    def run_qas(self, fasta_file, pdbdir, pkldir, output_dir_abs,
                pdbs_from_monomer, pdbs_from_multimer, pdbs_with_dist):

        result_dict = {}
        cwd = os.getcwd()

        pdb_count = len(os.listdir(pdbdir))

        if "alphafold" in self.run_methods:
            # if os.path.exists(output_dir_abs + '/alphafold_ranking.csv'):
            #     df = pd.read_csv(output_dir_abs + '/alphafold_ranking.csv')
            #     if len(df) != pdb_count:
            #         os.system(f"rm {output_dir_abs}/alphafold_ranking.csv")

            # if not os.path.exists(output_dir_abs + '/alphafold_ranking.csv'):
            alphafold_ranking = self.alphafold_qa.run(pkldir)
            alphafold_ranking.to_csv(output_dir_abs + '/alphafold_ranking.csv')
            result_dict["alphafold"] = output_dir_abs + '/alphafold_ranking.csv'

            # if not os.path.exists(output_dir_abs + '/alphafold_ranking_monomer.csv'):
            monomer_indices = [i for i in range(len(alphafold_ranking)) if
                               alphafold_ranking.loc[i, 'model'] in pdbs_from_monomer]
            alphafold_ranking_monomer = copy.deepcopy(alphafold_ranking.iloc[monomer_indices])
            alphafold_ranking_monomer.reset_index(inplace=True, drop=True)
            alphafold_ranking_monomer.to_csv(output_dir_abs + '/alphafold_ranking_monomer.csv')
            result_dict["alphafold_monomer"] = output_dir_abs + '/alphafold_ranking_monomer.csv'

            # if not os.path.exists(output_dir_abs + '/alphafold_ranking_multimer.csv'):
            monomer_indices = [i for i in range(len(alphafold_ranking)) if
                               alphafold_ranking.loc[i, 'model'] in pdbs_from_multimer]
            alphafold_ranking_multimer = copy.deepcopy(alphafold_ranking.iloc[monomer_indices])
            alphafold_ranking_multimer.reset_index(inplace=True, drop=True)
            alphafold_ranking_multimer.to_csv(output_dir_abs + '/alphafold_ranking_multimer.csv')
            result_dict["alphafold_multimer"] = output_dir_abs + '/alphafold_ranking_multimer.csv'

        if "apollo" in self.run_methods:
            if os.path.exists(output_dir_abs + '/pairwise_ranking.tm'):
                df = read_qa_txt_as_df(output_dir_abs + '/pairwise_ranking.tm')
                if len(df) != pdb_count:
                    os.system(f"rm {output_dir_abs}/pairwise_ranking.tm")

            if not os.path.exists(output_dir_abs + '/pairwise_ranking.tm'):
                os.chdir(output_dir_abs)
                with open("model.list", 'w') as fw:
                    for pdb in os.listdir(pdbdir):
                        fw.write(f"pdb/{pdb}\n")
                print(f"{self.parwise_qa} model.list {fasta_file} {self.tmscore} . pairwise_ranking")
                os.system(
                    f"{self.parwise_qa} model.list {fasta_file} {self.tmscore} . pairwise_ranking")
            result_dict["apollo"] = output_dir_abs + '/pairwise_ranking.tm'
            pairwise_ranking_df = read_qa_txt_as_df(output_dir_abs + '/pairwise_ranking.tm')
            monomer_indices = [i for i in range(len(pairwise_ranking_df)) if
                               pairwise_ranking_df.loc[i, 'model'] in pdbs_from_monomer]
            pairwise_ranking_monomer = copy.deepcopy(pairwise_ranking_df.iloc[monomer_indices])
            pairwise_ranking_monomer.reset_index(inplace=True, drop=True)
            pairwise_ranking_monomer.to_csv(output_dir_abs + '/pairwise_ranking_monomer.csv')
            result_dict["apollo_monomer"] = output_dir_abs + '/pairwise_ranking_monomer.csv'
       
        if "bfactor" in self.run_methods:
            # if os.path.exists(output_dir_abs + '/bfactor_ranking.csv'):
            #     df = pd.read_csv(output_dir_abs + '/bfactor_ranking.csv')
            #     if len(df) != pdb_count:
            #         os.system(f"rm {output_dir_abs}/bfactor_ranking.csv")
            # if not os.path.exists(output_dir_abs + '/bfactor_ranking.csv'):
            bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
            bfactor_ranking.to_csv(output_dir_abs + '/bfactor_ranking.csv')
            result_dict["bfactor"] = output_dir_abs + '/bfactor_ranking.csv'

        if "apollo" in self.run_methods and "alphafold" in self.run_methods:
            pairwise_ranking_df = read_qa_txt_as_df(result_dict["apollo"])
            ranks = [i + 1 for i in range(len(pairwise_ranking_df))]
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
                avg_rank = (int(avg_ranking_df.loc[i, 'pairwise_rank']) + int(
                    avg_ranking_df.loc[i, 'alphafold_rank'])) / 2
                avg_rankings += [avg_rank]
            avg_ranking_df['avg_score'] = avg_scores
            avg_ranking_df['avg_rank'] = avg_rankings
            avg_ranking_df = avg_ranking_df.sort_values(by=['avg_score'], ascending=False)
            avg_ranking_df.reset_index(inplace=True, drop=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="index"), axis=1, inplace=True)
            avg_ranking_df.drop(avg_ranking_df.filter(regex="Unnamed"), axis=1, inplace=True)
            avg_ranking_df.to_csv(output_dir_abs + '/pairwise_af_avg.ranking')
            result_dict["pairwise_af_avg"] = output_dir_abs + '/pairwise_af_avg.ranking'

            monomer_indices = [i for i in range(len(avg_ranking_df)) if
                               avg_ranking_df.loc[i, 'model'] in pdbs_from_monomer]
            avg_ranking_df_monomer = copy.deepcopy(avg_ranking_df.iloc[monomer_indices])
            avg_ranking_df_monomer.reset_index(inplace=True, drop=True)
            avg_ranking_df_monomer.to_csv(output_dir_abs + '/pairwise_af_avg_monomer.ranking')
            result_dict["pairwise_af_avg_monomer"] = output_dir_abs + '/pairwise_af_avg_monomer.ranking'

            monomer_indices = [i for i in range(len(avg_ranking_df)) if
                               avg_ranking_df.loc[i, 'model'] in pdbs_from_multimer]
            avg_ranking_df_multimer = copy.deepcopy(avg_ranking_df.iloc[monomer_indices])
            avg_ranking_df_multimer.reset_index(inplace=True, drop=True)
            avg_ranking_df_multimer.to_csv(output_dir_abs + '/pairwise_af_avg_multimer.ranking')
            result_dict["pairwise_af_avg_multimer"] = output_dir_abs + '/pairwise_af_avg_multimer.ranking'

        os.chdir(cwd)

        return result_dict

    def process(self, targetname, fasta_file, monomer_model_dir, output_dir, multimer_model_dir="", model_count=5):

        query_sequence = open(fasta_file).readlines()[1].rstrip('\n').strip()

        output_dir_abs = os.path.abspath(output_dir)

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb'
        clean_dir(pdbdir)

        pdbdir_monomer = output_dir + '/pdb_monomer'
        makedir_if_not_exists(pdbdir_monomer)

        pdbdir_multimer = output_dir + '/pdb_multimer'
        makedir_if_not_exists(pdbdir_multimer)

        pkldir = output_dir + '/pkl'
        clean_dir(pkldir)

        pkldir_monomer = output_dir + '/pkl_monomer'
        makedir_if_not_exists(pkldir_monomer)

        pkldir_multimer = output_dir + '/pkl_multimer'
        makedir_if_not_exists(pkldir_multimer)

        msadir = output_dir + '/msa'
        clean_dir(msadir)

        msadir_monomer = output_dir + '/msa_monomer'
        makedir_if_not_exists(msadir_monomer)

        msadir_multimer = output_dir + '/msa_multimer'
        makedir_if_not_exists(msadir_multimer)

        pdbs_with_dist = []

        pdbs_from_monomer = []
        if os.path.exists(monomer_model_dir):
            for method in os.listdir(monomer_model_dir):
                ranking_json_file = f"{monomer_model_dir}/{method}/ranking_debug.json"
                ranking_json = json.loads(open(ranking_json_file).read())
                for i in range(model_count):
                    os.system(f"cp {monomer_model_dir}/{method}/ranked_{i}.pdb {pdbdir_monomer}/{method}_{i}.pdb")

                    model_name = list(ranking_json["order"])[i]
                    has_distogram = extract_pkl(src_pkl=f"{monomer_model_dir}/{method}/result_{model_name}.pkl",
                                                output_pkl=f"{pkldir_monomer}/{method}_{i}.pkl")
                    os.system(
                        f"cp {monomer_model_dir}/{method}/msas/monomer_final.a3m {msadir_monomer}/{method}_{i}.a3m")
                    pdbs_from_monomer += [f"{method}_{i}.pdb"]

                    os.system(f"ln -s {pdbdir_monomer}/{method}_{i}.pdb {pdbdir}/{method}_{i}.pdb")
                    os.system(f"ln -s {pkldir_monomer}/{method}_{i}.pkl {pkldir}/{method}_{i}.pkl")
                    os.system(f"ln -s {msadir_monomer}/{method}_{i}.a3m {msadir}/{method}_{i}.a3m")

                    if has_distogram:
                        pdbs_with_dist += [f"{method}_{i}.pdb"]

        pdbs_from_multimer = []
        if os.path.exists(multimer_model_dir):
            for method in os.listdir(multimer_model_dir):
                ranking_json_file = f"{multimer_model_dir}/{method}/ranking_debug.json"
                if not os.path.exists(ranking_json_file):
                    continue
                ranking_json = json.loads(open(ranking_json_file).read())
                for i in range(model_count):
                    complex_pdb = f"{multimer_model_dir}/{method}/ranked_{i}.pdb"
                    if os.path.exists(complex_pdb):
                        chain_pdb_dict = extract_monomer_pdbs(complex_pdb=complex_pdb,
                                                              sequence=query_sequence,
                                                              output_prefix=f"{pdbdir_multimer}/{method}_{i}")
                        print(chain_pdb_dict)
                        for chain_id in chain_pdb_dict:
                            pdbname = chain_pdb_dict[chain_id]['pdbname']
                            print(pdbname)
                            model_name = list(ranking_json["order"])[i]
                            has_distogram = False
                            src_pkl = f"{multimer_model_dir}/{method}/result_{model_name}.pkl"
                            if os.path.exists(src_pkl):
                                has_distogram = extract_pkl(
                                        src_pkl=f"{multimer_model_dir}/{method}/result_{model_name}.pkl",
                                        residue_start=chain_pdb_dict[chain_id]['chain_start'],
                                        residue_end=chain_pdb_dict[chain_id]['chain_end'],
                                        output_pkl=pkldir_multimer + '/' + pdbname.replace('.pdb', '.pkl'))

                            os.system(f"cp {multimer_model_dir}/{method}/msas/{chain_id}/monomer_final.a3m "
                                      f"{msadir_multimer}/{pdbname.replace('.pdb', '.a3m')}")

                            # print(f"start: {residue_start}, end:{residue_end}")
                            pdbs_from_multimer += [pdbname]

                            os.system(f"ln -s {pdbdir_multimer}/{pdbname} {pdbdir}/{pdbname}")
                            os.system(
                                f"ln -s {pkldir_multimer}/{chain_pdb_dict[chain_id]['pdbname'].replace('.pdb', '.pkl')}"
                                f" {pkldir}/{chain_pdb_dict[chain_id]['pdbname'].replace('.pdb', '.pkl')}")
                            os.system(
                                f"ln -s {msadir_multimer}/{chain_pdb_dict[chain_id]['pdbname'].replace('.pdb', '.a3m')} "
                                f"{msadir}/{chain_pdb_dict[chain_id]['pdbname'].replace('.pdb', '.a3m')}")

                            if has_distogram:
                                pdbs_with_dist += [pdbname]

        return self.run_qas(fasta_file=fasta_file, pdbdir=pdbdir, pkldir=pkldir, output_dir_abs=output_dir_abs,
                            pdbs_from_monomer=pdbs_from_monomer, pdbs_from_multimer=pdbs_from_multimer,
                            pdbs_with_dist=pdbs_with_dist)

    def reprocess(self, targetname, fasta_file, output_dir):

        output_dir_abs = os.path.abspath(output_dir)

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb'
        pkldir = output_dir + '/pkl'
        msadir = output_dir + '/msa'

        if not os.path.exists(pdbdir + "_monomer"):
            pdbs_from_monomer = []
        else:
            pdbs_from_monomer = os.listdir(pdbdir + "_monomer")

        if not os.path.exists(pdbdir + "_multimer"):
            pdbs_from_multimer = []
        else:
            pdbs_from_multimer = os.listdir(pdbdir + "_multimer")

        return self.run_qas(fasta_file=fasta_file, pdbdir=pdbdir, pkldir=pkldir, output_dir_abs=output_dir_abs,
                            pdbs_from_monomer=pdbs_from_monomer,
                            pdbs_from_multimer=pdbs_from_multimer,
                            pdbs_with_dist=pdbs_from_monomer)
