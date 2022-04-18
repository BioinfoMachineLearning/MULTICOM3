import copy
import os
import sys
import time, json
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import dataclasses
from bml_casp15.quaternary_structure_refinement.iterative_refine_pipeline_v1_with_monomer import *
import pandas as pd


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

    def search(self, refinement_inputs, outdir):
        result_dirs = []

        pipeline = Multimer_iterative_refinement_pipeline(self.params)

        for refine_param in refinement_inputs:
            result_dir = pipeline.search_single(chain_id_map=refine_param.chain_id_map,
                                                fasta_path=refine_param.fasta_path, pdb_path=refine_param.pdb_path,
                                                pkl_path=refine_param.pkl_path, msa_paths=refine_param.msa_paths,
                                                outdir=outdir + '/' + pathlib.Path(refine_param.pdb_path).stem)
            result_dirs += [result_dir]

        return result_dirs


class Multimer_refinement_model_selection:

    def __init__(self, methods=[]):

        self.methods = methods

    def select_v1(self, indir, outdir):
        for pdb in os.listdir(indir):
            start_pdb = indir + '/' + pdb + '/iteration1/start.pdb'
            start_pkl = indir + '/' + pdb + '/iteration1/start.pkl'
            os.system(f"cp {start_pdb} {outdir}/{pdb}_ori.pdb")
            os.system(f"cp {start_pkl} {outdir}/{pdb}_ori.pkl")

            refine_pdb = indir + '/' + pdb + '/final/final.pdb'
            refine_pkl = indir + '/' + pdb + '/final/final.pkl'
            os.system(f"cp {refine_pdb} {outdir}/{pdb}_ref.pdb")
            os.system(f"cp {refine_pkl} {outdir}/{pdb}_ref.pdb")

        pdbs = []
        confidences = []
        for pkl in os.listdir(outdir):
            if pkl.find('.pkl') > 0:
                with open(outdir + '/' + pkl, 'rb') as f:
                    prediction_result = pickle.load(f)
                    pdbs += [pkl.replace('.pkl', '.pdb')]
                    confidences += [prediction_result['ranking_confidence']]

        df = pd.DataFrame({'model': pdbs, 'confidence': confidences})
        df = df.sort_values(by=['confidence'], ascending=False)
        df.reset_index(inplace=True)

        for i in range(5):
            pdb_name = df.loc[i, 'model']
            os.system(f"cp {outdir}/{pdb_name} {outdir}/casp{i + 1}.pdb")
            os.system(f"cp {outdir}/{pdb_name.replace('.pdb', '.pkl')} {outdir}/casp{i + 1}.pkl")

        return outdir

    def select_v2(self, ranking_df, indir, outdir):
        for i in range(5):
            pdb_name = ranking_df.loc[i, 'model']
            start_pdb = indir + '/' + pdb_name + '/iteration1/start.pdb'
            start_pkl = indir + '/' + pdb_name + '/iteration1/start.pkl'
            os.system(f"cp {start_pdb} {outdir}/{pdb_name}_ori.pdb")
            os.system(f"cp {start_pkl} {outdir}/{pdb_name}_ori.pkl")
            with open(start_pkl, 'rb') as f:
                plddt_start = pickle.load(f)['confidence']

            refine_pdb = indir + '/' + pdb_name + '/final/final.pdb'
            refine_pkl = indir + '/' + pdb_name + '/final/final.pkl'
            os.system(f"cp {refine_pdb} {outdir}/{pdb_name}_ref.pdb")
            os.system(f"cp {refine_pkl} {outdir}/{pdb_name}_ref.pdb")
            with open(refine_pkl, 'rb') as f:
                plddt_ref = pickle.load(f)['confidence']

            if plddt_start > plddt_ref:
                os.system(f"cp {outdir}/{pdb_name}_ori.pdb {outdir}/casp{i + 1}.pdb")
                os.system(f"cp {outdir}/{pdb_name}_ori.pkl {outdir}/casp{i + 1}.pkl")
            else:
                os.system(f"cp {outdir}/{pdb_name}_ref.pdb {outdir}/casp{i + 1}.pdb")
                os.system(f"cp {outdir}/{pdb_name}_ref.pkl {outdir}/casp{i + 1}.pkl")
        return outdir
