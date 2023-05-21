import os, sys, argparse, time, json
from multiprocessing import Pool
from tqdm import tqdm
from multicom3.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
import pandas as pd
from multicom3.monomer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from multicom3.quaternary_structure_evaluation.pairwise_dockq import Pairwise_dockq_qa
from multicom3.quaternary_structure_evaluation.dproq_ranking import DPROQ
from multicom3.quaternary_structure_evaluation.enqa_ranking import En_qa
from multicom3.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa
from multicom3.quaternary_structure_evaluation.multieva_qa import MultiEva_qa
from multicom3.quaternary_structure_evaluation.foldseek_ranking import FoldSeek_qa
from multicom3.common.protein import complete_result
from multicom3.monomer_structure_evaluation.pipeline_sep import extract_monomer_pdbs
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import pickle

PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'


def combine_pdb(pdbs, combine_pdb):
    with open(combine_pdb, 'w') as out:
        for chain_id, pdb in zip(PDB_CHAIN_IDS, pdbs):
            for line in open(pdb):
                if not line.startswith('ATOM'):
                    continue
                out.write(line[:21] + chain_id + line[22:])
            out.write('TER\n')


def extract_multimer_pkl(inpkl, outpkl, chain_starts, chain_ends):
    with open(inpkl, 'rb') as f:
        prediction_result = pickle.load(f)
        indices = []
        for chain_start, chain_end in zip(chain_starts, chain_ends):
            indices += list(range(chain_start, chain_end + 1))
        prediction_result_new = {'plddt': prediction_result['plddt'][indices,]}
        if 'distogram' in prediction_result:
            distogram_new = prediction_result['distogram']
            ttt1 = distogram_new['logits'][indices, :, :]
            distogram_new['logits'] = ttt1[:, indices, :]
            prediction_result_new['distogram'] = distogram_new
        with open(outpkl, 'wb') as f:
            pickle.dump(prediction_result_new, f, protocol=4)


def cal_contact_ratio_between_pdbs(pdb1, pdb2):
    parser = PDBParser(PERMISSIVE=1)
    structure1 = parser.get_structure('', pdb1)
    structure2 = parser.get_structure('', pdb2)

    model1 = structure1[0]
    chain_id1 = list(model1.child_dict.keys())
    xyzPDB1 = model1[chain_id1[0]]

    model2 = structure2[0]
    chain_id2 = list(model2.child_dict.keys())
    xyzPDB2 = model2[chain_id2[0]]

    h_dist_map = np.zeros((len(xyzPDB1), len(xyzPDB2)))

    for i in range(1, len(xyzPDB1) + 1):
        for j in range(1, len(xyzPDB2) + 1):
            res_i = xyzPDB1[i]
            res_j = xyzPDB2[j]
            dist_list = []
            for atom_i in res_i:
                for atom_j in res_j:
                    if ('C' in atom_i.name or 'N' in atom_i.name or 'O' in atom_i.name or 'S' in atom_i.name) and \
                            ('C' in atom_j.name or 'N' in atom_j.name or 'O' in atom_j.name or 'S' in atom_j.name):
                        dist_list.append(atom_i - atom_j)
                    else:
                        continue
            min_dist = np.min(dist_list)
            h_dist_map[i - 1, j - 1] = min_dist

    return (h_dist_map < 6).sum() / (len(xyzPDB1) * len(xyzPDB2))


def pair_chain_ids(src_pair_dict, trg_chain_pdb_dict):
    chain_ids_1 = []
    chain_ids_2 = []
    ratios = []
    for src_pair_id in src_pair_dict:
        for trg_chain_id in trg_chain_pdb_dict:
            max_ratio = 0
            pdb2 = trg_chain_pdb_dict[trg_chain_id]['pdb']
            for src_chain_idx, src_chain_id in enumerate(src_pair_id.split('_')):
                pdb1 = src_pair_dict[src_pair_id]['pdb'][src_chain_idx]
                ratio = cal_contact_ratio_between_pdbs(pdb1, pdb2)
                if ratio > max_ratio:
                    max_ratio = ratio
            chain_ids_1 += [src_pair_id]
            chain_ids_2 += [trg_chain_id]
            ratios += [max_ratio]

    df = pd.DataFrame({'chain1': chain_ids_1, 'chain2': chain_ids_2, 'ratio': ratios})
    df = df.sort_values(by=['ratio'], ascending=False)
    df.reset_index(inplace=True, drop=True)
    print(df)

    processed_ids = []
    res_pair_dict = {}
    for i in range(len(df)):
        chain1 = df.loc[i, 'chain1']
        chain2 = df.loc[i, 'chain2']
        if chain1 in processed_ids or chain2 in processed_ids:
            continue

        res_pair_dict[f"{chain1}_{chain2}"] = {'pdb': [], 'chain_start': [], 'chain_end': []}
        res_pair_dict[f"{chain1}_{chain2}"]['pdb'] += src_pair_dict[chain1]['pdb']
        res_pair_dict[f"{chain1}_{chain2}"]['pdb'] += [trg_chain_pdb_dict[chain2]['pdb']]

        res_pair_dict[f"{chain1}_{chain2}"]['chain_start'] += src_pair_dict[chain1]['chain_start']
        res_pair_dict[f"{chain1}_{chain2}"]['chain_start'] += [trg_chain_pdb_dict[chain2]['chain_start']]

        res_pair_dict[f"{chain1}_{chain2}"]['chain_end'] += src_pair_dict[chain1]['chain_end']
        res_pair_dict[f"{chain1}_{chain2}"]['chain_end'] += [trg_chain_pdb_dict[chain2]['chain_end']]

        processed_ids += [chain1]
        processed_ids += [chain2]

    return res_pair_dict


def extract_multimer_pdbs(chain_id_map, complex_pdb, workdir, complex_pkl,
                          output_pdb_name, output_pkl_name):
    makedir_if_not_exists(workdir)
    multimer_pdbs = {}
    for chain_id in chain_id_map:
        chain_pdb_dict = extract_monomer_pdbs(complex_pdb=complex_pdb,
                                              sequence=chain_id_map[chain_id].sequence,
                                              output_prefix=f"{workdir}/{chain_id}")
        print(chain_pdb_dict)
        multimer_pdbs[chain_id] = chain_pdb_dict

    src_pair_dict = {}
    for chain_idx, chain_id in enumerate(multimer_pdbs):
        if chain_idx == 0:
            for same_chain_id in multimer_pdbs[chain_id]:
                src_pair_dict[same_chain_id] = {'pdb': [multimer_pdbs[chain_id][same_chain_id]['pdb']],
                                                'chain_start': [multimer_pdbs[chain_id][same_chain_id]['chain_start']],
                                                'chain_end': [multimer_pdbs[chain_id][same_chain_id]['chain_end']]}
        else:
            src_pair_dict = pair_chain_ids(src_pair_dict, multimer_pdbs[chain_id])

    print(src_pair_dict)

    for pair_idx, pair in enumerate(src_pair_dict):
        combine_pdb([src_pair_dict[pair]['pdb'][chain_idx] for chain_idx, chain_id in enumerate(pair.split('_'))],
                    output_pdb_name + pair + '.pdb')
        # if os.path.exists(complex_pkl):
        #     extract_multimer_pkl(complex_pkl, output_pkl_name + pair + '.pkl',
        #                 src_pair_dict[pair]['chain_start'], src_pair_dict[pair]['chain_end'])


class Quaternary_structure_evaluation_pipeline_human:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "bfactor", 'multieva']):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods

        self.pairwise_qa = Pairwise_dockq_qa(params['dockq_program'])
        self.alphafold_qa = Alphafold_pkl_qa(sort_field='confidence')
        self.dproq = DPROQ(dproq_program=params['dproq_program'])
        self.enqa = En_qa(enqa_program=params['enqa_program'])
        self.bfactorqa = Bfactor_qa()
        self.multieva = MultiEva_qa(multieva_program=params['multieva_program'])
        self.foldseek_qa = FoldSeek_qa(params=params)

    def process(self, fasta_path, chain_id_map, model_dir, extract_model_dir, output_dir,
                monomer_model_dir="", stoichiometry="", model_count=5):

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb/'
        makedir_if_not_exists(pdbdir)

        pkldir = output_dir + '/pkl/'
        makedir_if_not_exists(pkldir)

        msadir = output_dir + '/msa/'
        makedir_if_not_exists(msadir)

        for method in os.listdir(model_dir):
            ranking_json_file = f"{model_dir}/{method}/ranking_debug.json"
            if not os.path.exists(ranking_json_file):
                continue
            ranking_json = json.loads(open(ranking_json_file).read())

            # if not complete_result(model_dir + '/' + method):
            #     continue

            for i in range(model_count):
                if i >= len(list(ranking_json["order"])):
                    break
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

        if os.path.exists(extract_model_dir):
            for method in os.listdir(extract_model_dir):
                ranking_json_file = f"{extract_model_dir}/{method}/ranking_debug.json"
                if not os.path.exists(ranking_json_file):
                    continue
                ranking_json = json.loads(open(ranking_json_file).read())

                for i in range(model_count):
                    model_name = list(ranking_json["order"])[i]
                    complex_pdb = f"{extract_model_dir}/{method}/ranked_{i}.pdb"
                    complex_pkl = f"{extract_model_dir}/{method}/result_{model_name}.pkl"

                    if os.path.exists(complex_pdb):
                        extract_multimer_pdbs(chain_id_map=chain_id_map,
                                              complex_pdb=complex_pdb,
                                              complex_pkl=complex_pkl,
                                              workdir=f"{extract_model_dir}/{method}/ranked_{i}",
                                              output_pdb_name=f"{pdbdir}/{method}_{i}",
                                              output_pkl_name=f"{pkldir}/{method}_{i}")

        result_dict = {}

        if "pairwise" in self.run_methods:
            if not os.path.exists(output_dir + '/pairwise_ranking.csv'):
                pairwise_ranking = self.pairwise_qa.run(pdbdir)
                pairwise_ranking.to_csv(output_dir + '/pairwise_ranking.csv')
            result_dict["pairwise"] = output_dir + '/pairwise_ranking.csv'

        if "alphafold" in self.run_methods:
            if not os.path.exists(output_dir + '/alphafold_ranking.csv'):
                alphafold_ranking = self.alphafold_qa.run(pkldir)
                alphafold_ranking.to_csv(output_dir + '/alphafold_ranking.csv')
            result_dict["alphafold"] = output_dir + '/alphafold_ranking.csv'

        if "dproq" in self.run_methods:
            if not os.path.exists(output_dir + '/DOCKQ_ranking.csv'):
                dproq_ranking_dockq, dproq_ranking_evalue = self.dproq.run(indir=pdbdir, outdir=output_dir + '/dproq')
                dproq_ranking_dockq.to_csv(output_dir + '/dproq_ranking_dockq.csv')
                dproq_ranking_evalue.to_csv(output_dir + '/dproq_ranking_evalue.csv')
            result_dict["dproq_ranking_dockq"] = output_dir + '/dproq_ranking_dockq.csv'
            result_dict["dproq_ranking_evalue"] = output_dir + '/dproq_ranking_evalue.csv'

        if "enqa" in self.run_methods:
            if not os.path.exists(output_dir + '/enqa_ranking.csv'):
                enqa_ranking = self.enqa.run_with_pairwise_ranking(input_dir=pdbdir,
                                                                   pkl_dir=pkldir,
                                                                   pairwise_ranking_file=output_dir + '/pairwise_ranking.csv',
                                                                   outputdir=output_dir + '/enqa')
                enqa_ranking.to_csv(output_dir + '/enqa_ranking.csv')
            result_dict["enQA"] = output_dir + '/enqa_ranking.csv'

        if "bfactor" in self.run_methods:
            if not os.path.exists(output_dir + '/bfactor_ranking.csv'):
                bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
                bfactor_ranking.to_csv(output_dir + '/bfactor_ranking.csv')
            result_dict["bfactor"] = output_dir + '/bfactor_ranking.csv'

        if "multieva" in self.run_methods:
            if not os.path.exists(f"{output_dir}/multieva.csv"):
                workdir = output_dir + '/multieva'
                makedir_if_not_exists(workdir)
                refdir = workdir + '/monomer_af/'
                makedir_if_not_exists(refdir)

                for chain_id in chain_id_map:
                    default_chain_model = monomer_model_dir + '/' + chain_id_map[
                        chain_id].description + '/pdb/default_0.pdb'
                    os.system(f"cp {default_chain_model} {refdir}/{chain_id_map[chain_id].description}.pdb")

                multieva_csv = self.multieva.run(chain_id_map=chain_id_map,
                                                 fasta_path=fasta_path,
                                                 input_dir=pdbdir, stoichiometry=stoichiometry,
                                                 alphafold_prediction_dir=refdir, outputdir=workdir)

                os.system(f"cp {multieva_csv} {output_dir}/multieva.csv")
            result_dict["multieva"] = output_dir + '/multieva.csv'

        if "foldseek" in self.run_methods:
            if not os.path.exists(f"{output_dir}/foldseek_qa.csv"):
                foldseek_qa_df = self.foldseek_qa.run(chain_id_map=chain_id_map, input_dir=pdbdir,
                                                      outputdir=output_dir + '/foldseek')
                foldseek_qa_df.to_csv(f"{output_dir}/foldseek_qa.csv")
            result_dict["foldseek"] = output_dir + '/foldseek_qa.csv'

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
                pairwise_score = float(avg_ranking_df.loc[i, 'average_MMS'])
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
                pairwise_score = float(avg_ranking_df.loc[i, 'average_MMS'])
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

        if "pairwise" in self.run_methods:
            if not os.path.exists(output_dir + '/pairwise_ranking.csv'):
                pairwise_ranking = self.pairwise_qa.run(pdbdir)
                pairwise_ranking.to_csv(output_dir + '/pairwise_ranking.csv')
            result_dict["pairwise"] = output_dir + '/pairwise_ranking.csv'

        if "alphafold" in self.run_methods:
            if not os.path.exists(output_dir + '/alphafold_ranking.csv'):
                alphafold_ranking = self.alphafold_qa.run(pkldir)
                alphafold_ranking.to_csv(output_dir + '/alphafold_ranking.csv')
            result_dict["alphafold"] = output_dir + '/alphafold_ranking.csv'

        if "dproq" in self.run_methods:
            if not os.path.exists(output_dir + '/DOCKQ_ranking.csv'):
                dproq_ranking_dockq, dproq_ranking_evalue = self.dproq.run(indir=pdbdir, outdir=output_dir + '/dproq')
                dproq_ranking_dockq.to_csv(output_dir + '/dproq_ranking_dockq.csv')
                dproq_ranking_evalue.to_csv(output_dir + '/dproq_ranking_evalue.csv')
            result_dict["dproq_ranking_dockq"] = output_dir + '/dproq_ranking_dockq.csv'
            result_dict["dproq_ranking_evalue"] = output_dir + '/dproq_ranking_evalue.csv'

        if "enqa" in self.run_methods:
            if not os.path.exists(output_dir + '/enqa_ranking.csv'):
                enqa_ranking = self.enqa.run_with_pairwise_ranking(input_dir=pdbdir,
                                                                   pkl_dir=pkldir,
                                                                   pairwise_ranking_file=output_dir + '/pairwise_ranking.csv',
                                                                   outputdir=output_dir + '/enqa')
                enqa_ranking.to_csv(output_dir + '/enqa_ranking.csv')
            result_dict["enQA"] = output_dir + '/enqa_ranking.csv'

        if "bfactor" in self.run_methods:
            if not os.path.exists(output_dir + '/bfactor_ranking.csv'):
                bfactor_ranking = self.bfactorqa.run(input_dir=pdbdir)
                bfactor_ranking.to_csv(output_dir + '/bfactor_ranking.csv')
            result_dict["bfactor"] = output_dir + '/bfactor_ranking.csv'

        if "multieva" in self.run_methods:
            if not os.path.exists(f"{output_dir}/multieva.csv"):
                workdir = output_dir + '/multieva'
                makedir_if_not_exists(workdir)
                refdir = workdir + '/monomer_af/'
                makedir_if_not_exists(refdir)

                for chain_id in chain_id_map:
                    default_chain_model = monomer_model_dir + '/' + chain_id_map[
                        chain_id].description + '/pdb/default_0.pdb'
                    os.system(f"cp {default_chain_model} {refdir}/{chain_id_map[chain_id].description}.pdb")

                multieva_csv = self.multieva.run(chain_id_map=chain_id_map,
                                                 fasta_path=fasta_path,
                                                 input_dir=pdbdir, stoichiometry=stoichiometry,
                                                 alphafold_prediction_dir=refdir, outputdir=workdir)

                os.system(f"cp {multieva_csv} {output_dir}/multieva.csv")
            result_dict["multieva"] = output_dir + '/multieva.csv'

        if "foldseek" in self.run_methods:
            if not os.path.exists(f"{output_dir}/foldseek_qa.csv"):
                foldseek_qa_df = self.foldseek_qa.run(chain_id_map=chain_id_map, input_dir=pdbdir,
                                                      outputdir=output_dir + '/foldseek')
                foldseek_qa_df.to_csv(f"{output_dir}/foldseek_qa.csv")
            result_dict["foldseek"] = output_dir + '/foldseek_qa.csv'

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
                pairwise_score = float(avg_ranking_df.loc[i, 'average_MMS'])
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
                pairwise_score = float(avg_ranking_df.loc[i, 'average_MMS'])
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
