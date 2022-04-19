import copy
import os
import sys
import time, json
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import dataclasses
from bml_casp15.monomer_structure_refinement.iterative_refine_pipeline_v4_50 import *
import pandas as pd


class refinement_input:
    def __init__(self, fasta_path, pdb_path, pkl_path, msa_path):
        self.fasta_path = fasta_path
        self.pdb_path = pdb_path
        self.pkl_path = pkl_path
        self.msa_path = msa_path


class Monomer_iterative_refinement_pipeline_server:

    def __init__(self, params):
        self.params = params

    def search(self, refinement_inputs, outdir):
        result_dirs = []

        pipeline = Monomer_iterative_refinement_pipeline(self.params)

        for refine_param in refinement_inputs:
            result_dir = pipeline.search_single(fasta_path=refine_param.fasta_path, pdb_path=refine_param.pdb_path,
                                                pkl_path=refine_param.pkl_path, msa_path=refine_param.msa_path,
                                                outdir=outdir + '/' + pathlib.Path(refine_param.pdb_path).stem)
            result_dirs += [result_dir]

        return result_dirs


class Monomer_refinement_model_selection:

    def __init__(self, methods=[]):

        self.methods = methods

    def select_v1(self, indir, outdir):
        if os.path.exists(outdir):
            os.system(f"rm -rf {outdir}")
        makedir_if_not_exists(outdir)

        for pdb in os.listdir(indir):
            if not os.path.exists(indir + '/' + pdb + '/iteration1'):
                continue
            start_pdb = indir + '/' + pdb + '/iteration1/start.pdb'
            start_pkl = indir + '/' + pdb + '/iteration1/start.pkl'
            os.system(f"cp {start_pdb} {outdir}/{pdb}_ori.pdb")
            os.system(f"cp {start_pkl} {outdir}/{pdb}_ori.pkl")

            refine_pdb = indir + '/' + pdb + '/final/final.pdb'
            refine_pkl = indir + '/' + pdb + '/final/final.pkl'
            os.system(f"cp {refine_pdb} {outdir}/{pdb}_ref.pdb")
            os.system(f"cp {refine_pkl} {outdir}/{pdb}_ref.pkl")

        pdbs = []
        plddts = []
        for pkl in os.listdir(outdir):
            if pkl.find('.pkl') > 0:
                with open(outdir + '/' + pkl, 'rb') as f:
                    prediction_result = pickle.load(f)
                    pdbs += [pkl.replace('.pkl', '.pdb')]
                    plddts += [np.mean(prediction_result['plddt'])]

        df = pd.DataFrame({'model': pdbs, 'plddt': plddts})
        df = df.sort_values(by=['plddt'], ascending=False)
        df.reset_index(inplace=True, drop=True)
        df.to_csv(outdir + '/final_ranking.csv')

        for i in range(5):
            pdb_name = df.loc[i, 'model']
            os.system(f"cp {outdir}/{pdb_name} {outdir}/casp{i + 1}.pdb")
            os.system(f"cp {outdir}/{pdb_name.replace('.pdb', '.pkl')} {outdir}/casp{i + 1}.pkl")

        return outdir

    def select_v2(self, ranking_df, indir, outdir):
        if os.path.exists(outdir):
            os.system(f"rm -rf {outdir}")
        makedir_if_not_exists(outdir)

        start_pdbs = []
        refine_pdbs = []
        start_plddts = []
        refine_plddts = []
        final_models = []
        for i in range(5):
            pdb_name = ranking_df.loc[i, 'model'].replace('.pdb', '')
            start_pdb = indir + '/' + pdb_name + '/iteration1/start.pdb'
            start_pkl = indir + '/' + pdb_name + '/iteration1/start.pkl'
            os.system(f"cp {start_pdb} {outdir}/{pdb_name}_ori.pdb")
            os.system(f"cp {start_pkl} {outdir}/{pdb_name}_ori.pkl")
            with open(start_pkl, 'rb') as f:
                plddt_start = np.mean(pickle.load(f)['plddt'])

            start_pdbs += [f"{pdb_name}_ori.pdb"]
            start_plddts += [plddt_start]

            refine_pdb = indir + '/' + pdb_name + '/final/final.pdb'
            refine_pkl = indir + '/' + pdb_name + '/final/final.pkl'
            os.system(f"cp {refine_pdb} {outdir}/{pdb_name}_ref.pdb")
            os.system(f"cp {refine_pkl} {outdir}/{pdb_name}_ref.pkl")
            with open(refine_pkl, 'rb') as f:
                plddt_ref = np.mean(pickle.load(f)['plddt'])

            refine_pdbs += [f"{pdb_name}_ref.pdb"]
            refine_plddts += [plddt_ref]

            if plddt_start > plddt_ref:
                os.system(f"cp {outdir}/{pdb_name}_ori.pdb {outdir}/casp{i + 1}.pdb")
                os.system(f"cp {outdir}/{pdb_name}_ori.pkl {outdir}/casp{i + 1}.pkl")
                final_models += [f"{pdb_name}_ori.pdb"]
            else:
                os.system(f"cp {outdir}/{pdb_name}_ref.pdb {outdir}/casp{i + 1}.pdb")
                os.system(f"cp {outdir}/{pdb_name}_ref.pkl {outdir}/casp{i + 1}.pkl")
                final_models += [f"{pdb_name}_ref.pdb"]

        df = pd.DataFrame({'start_model': start_pdbs, 'start_plddt': start_plddts,
                           'refine_model': refine_pdbs, 'refine_plddt': refine_plddts, 'final_model': final_models})
        df.to_csv(outdir + '/final_ranking.csv')

        return outdir
