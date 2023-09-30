import copy
import os
import sys
import time, json
from multicom3.common.util import makedir_if_not_exists, check_dirs
import dataclasses
from multicom3.multimer_structure_refinement import iterative_refine_pipeline_heteromer_v1_with_monomer
from multicom3.multimer_structure_refinement import iterative_refine_pipeline_homo_v1
import pandas as pd
import pathlib
import pickle


class refinement_input_multimer:
    def __init__(self, chain_id_map, fasta_path, pdb_path, pkl_path, msa_paths):
        self.chain_id_map = chain_id_map
        self.fasta_path = fasta_path
        self.pdb_path = pdb_path
        self.pkl_path = pkl_path
        self.msa_paths = msa_paths


class Multimer_iterative_refinement_pipeline_server:

    def __init__(self, params):
        self.params = params

    def search(self, refinement_inputs, outdir, stoichiometry):
        result_dirs = []
        for refine_param in refinement_inputs:
            result_dir = ""
            if stoichiometry == "homomer":
                print("refining homomers")
                pipeline_v1 = iterative_refine_pipeline_homo_v1.Multimer_iterative_refinement_pipeline(self.params)
                result_dir = pipeline_v1.search_single(chain_id_map=refine_param.chain_id_map,
                                                            fasta_path=refine_param.fasta_path,
                                                            pdb_path=refine_param.pdb_path,
                                                            pkl_path=refine_param.pkl_path,
                                                            msa_paths=refine_param.msa_paths,
                                                            outdir=os.path.join(outdir, pathlib.Path(refine_param.pdb_path).stem))
            elif stoichiometry == "heteromer":
                print("refining heteromer")
                pipeline_v1_with_monomer = iterative_refine_pipeline_heteromer_v1_with_monomer.Multimer_iterative_refinement_pipeline(
                    self.params)
                result_dir = pipeline_v1_with_monomer.search_single(chain_id_map=refine_param.chain_id_map,
                                                                    fasta_path=refine_param.fasta_path,
                                                                    pdb_path=refine_param.pdb_path,
                                                                    pkl_path=refine_param.pkl_path,
                                                                    msa_paths=refine_param.msa_paths,
                                                                    outdir=os.path.join(outdir, pathlib.Path(refine_param.pdb_path).stem))

            result_dirs += [result_dir]

        return result_dirs


class Multimer_refinement_model_selection:

    def __init__(self, methods=[]):

        self.methods = methods

    def select_v1(self, indir, outdir):
        for pdb in os.listdir(indir):
            start_pdb = os.path.join(indir, pdb, 'iteration1/start.pdb')
            start_pkl = os.path.join(indir, pdb, 'iteration1/start.pkl')
            os.system(f"cp {start_pdb} {os.path.join(outdir, pdb + "_ori.pdb")}")
            os.system(f"cp {start_pkl} {os.path.join(outdir, pdb + "_ori.pkl")}")

            refine_pdb = os.path.join(indir, pdb, 'final/final.pdb')
            refine_pkl = os.path.join(indir, pdb, 'final/final.pkl')
            os.system(f"cp {refine_pdb} " + os.path.join(outdir, pdb + "_ref.pdb"))
            os.system(f"cp {refine_pkl} " + os.path.join(outdir, pdb + "_ref.pkl"))

        pdbs = []
        confidences = []
        for pkl in os.listdir(outdir):
            if pkl.find('.pkl') > 0:
                with open(os.path.join(outdir, pkl), 'rb') as f:
                    prediction_result = pickle.load(f)
                    pdbs += [pkl.replace('.pkl', '.pdb')]
                    confidences += [prediction_result['ranking_confidence']]

        df = pd.DataFrame({'model': pdbs, 'confidence': confidences})
        df = df.sort_values(by=['confidence'], ascending=False)
        df.reset_index(inplace=True)

        for i in range(5):
            pdb_name = df.loc[i, 'model']
            deep_name = f"deep{i+1}"
            os.system("cp " + os.path.join(outdir, pdb_name) + " " + os.path.join(outdir, deep_name + ".pdb"))
            os.system("cp " + os.path.join(outdir, pdb_name.replace('.pdb', '.pkl')) + " " + os.path.join(outdir, deep_name + ".pkl"))

        df.to_csv(os.path.join(outdir, 'final_ranking.csv'))
        return outdir

    def select_v2(self, ranking_df, indir, outdir):
        for i in range(5):
            pdb_name = ranking_df.loc[i, 'model']
            start_pdb = os.path.join(indir, pdb_name, 'iteration1/start.pdb')
            start_pkl = os.path.join(indir, pdb_name, 'iteration1/start.pkl')
            os.system(f"cp {start_pdb} " + os.path.join(outdir, pdb_name + "_ori.pdb"))
            os.system(f"cp {start_pkl} " + os.path.join(outdir, pdb_name + "_ori.pkl"))

            with open(start_pkl, 'rb') as f:
                plddt_start = pickle.load(f)['confidence']

            refine_pdb = os.path.join(indir, pdb_name, 'final/final.pdb')
            refine_pkl = os.path.join(indir, pdb_name, 'final/final.pkl')
            os.system(f"cp {refine_pdb} " + os.path.join(outdir, pdb_name + "_ref.pdb"))
            os.system(f"cp {refine_pkl} " + os.path.join(outdir, pdb_name + "_ref.pkl"))

            with open(refine_pkl, 'rb') as f:
                plddt_ref = pickle.load(f)['confidence']

            casp_pdb_name = f"casp{i+1}"
            if plddt_start > plddt_ref:
                os.system("cp " + os.path.join(outdir, pdb_name + "_ori.pdb") + " " + os.path.join(outdir, casp_pdb_name + ".pdb"))
                os.system("cp " + os.path.join(outdir, pdb_name + "_ori.pkl") + " " + os.path.join(outdir, casp_pdb_name + ".pkl"))
            else:
                os.system("cp " + os.path.join(outdir, pdb_name + "_ref.pdb") + " " + os.path.join(outdir, casp_pdb_name + ".pdb"))
                os.system("cp " + os.path.join(outdir, pdb_name + "_ref.pkl") + " " + os.path.join(outdir, casp_pdb_name + ".pkl"))
        return outdir
