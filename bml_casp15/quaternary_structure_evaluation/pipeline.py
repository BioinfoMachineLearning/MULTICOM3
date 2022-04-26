import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
import pandas as pd
from bml_casp15.quaternary_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa
from bml_casp15.quaternary_structure_evaluation.pairwise_dockq import Pairwise_dockq_qa
from bml_casp15.quaternary_structure_evaluation.dproq_ranking import DPROQ
from bml_casp15.quaternary_structure_evaluation.enqa_ranking import En_qa
from bml_casp15.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa
from bml_casp15.quaternary_structure_evaluation.multieva_qa import MultiEva_qa
from bml_casp15.quaternary_structure_evaluation.foldseek_ranking import FoldSeek_qa
from bml_casp15.common.protein import complete_result


class Quaternary_structure_evaluation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self, params, run_methods=["alphafold", "pairwise", "enqa", "dproq", "bfactor", 'multieva', 'foldseek']):
        """Initializes the data pipeline."""

        self.params = params
        self.run_methods = run_methods

        self.pairwise_qa = Pairwise_dockq_qa(params['dockq_program'])
        self.alphafold_qa = Alphafold_pkl_qa()
        self.dproq = DPROQ(dproq_program=params['dproq_program'])
        self.enqa = En_qa(enqa_program=params['enqa_program'])
        self.bfactorqa = Bfactor_qa()
        self.multieva = MultiEva_qa(multieva_program=params['multieva_program'])
        self.foldseek_qa = FoldSeek_qa(params=params)

    def process(self, fasta_path, chain_id_map, model_dir, output_dir, monomer_model_dir = "", stoichiometry = ""):

        makedir_if_not_exists(output_dir)

        pdbdir = output_dir + '/pdb/'
        makedir_if_not_exists(pdbdir)

        pkldir = output_dir + '/pkl/'
        makedir_if_not_exists(pkldir)

        msadir = output_dir + '/msa/'
        makedir_if_not_exists(msadir)

        for method in os.listdir(model_dir):
            for i in range(0, 5):
                if not complete_result(model_dir + '/' + method):
                    continue
                os.system(f"cp {model_dir}/{method}/ranked_{i}.pdb {pdbdir}/{method}_{i}.pdb")
                os.system(f"cp {model_dir}/{method}/result_model_{i + 1}_multimer.pkl {pkldir}/{method}_{i}.pkl")
                for chain_id in chain_id_map:
                    msa_chain_outdir = msadir + '/' + chain_id_map[chain_id].description
                    makedir_if_not_exists(msa_chain_outdir)
                    os.system(f"cp {model_dir}/{method}/msas/{chain_id}/monomer_final.a3m "
                              f"{msa_chain_outdir}/{method}_{i}.monomer.a3m")
                    os.system(f"cp {model_dir}/{method}/msas/{chain_id_map[chain_id].description}.paired.a3m "
                              f"{msa_chain_outdir}/{method}_{i}.paired.a3m")

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
                dproq_ranking_dockq, dproq_ranking_evalue = self.dproq.run(indir=pdbdir, outdir=output_dir+'/dproq')
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
            
        return result_dict
