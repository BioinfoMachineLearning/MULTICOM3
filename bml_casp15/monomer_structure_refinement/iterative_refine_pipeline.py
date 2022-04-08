import copy
import os
import sys
import time, json
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
from bml_casp15.tool.foldseek import *
import pickle
import numpy as np
from bml_casp15.monomer_structure_refinement.iterative_refine_pipeline_v4_50 import *


class Monomer_iterative_refinement_pipeline:

    def __init__(self, params):

        self.params = params

    def search(self, fasta_file, input_pdb_dir, outdir, native_pdb=""):

        input_pdb_dir = os.path.abspath(input_pdb_dir)

        fasta_file = os.path.abspath(fasta_file)

        query_sequence = ""
        for line in open(fasta_file):
            line = line.rstrip('\n')
            if line.startswith('>'):
                continue
            else:
                query_sequence = line

        targetname = pathlib.Path(fasta_file).stem

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        iteration_scores = {}

        true_tm_scores = {}

        iteration_result_all = {'targetname': [],
                                'model': [],
                                'start_lddt': [],
                                'end_lddt': [],
                                'start_tmscore': [],
                                'end_tmscore': []}

        iteration_result_avg = {'targetname': [targetname], 'start_lddt': [], 'end_lddt': [], 'start_tmscore': [],
                                'end_tmscore': []}

        iteration_result_max = {'targetname': [targetname], 'start_lddt': [], 'end_lddt': [], 'start_tmscore': [],
                                'end_tmscore': []}

        cwd = os.getcwd()

        for i in range(0, 5):
            model_outdir = f"{outdir}/ranked_{i}"
            makedir_if_not_exists(model_outdir)

            current_ref_dir = input_pdb_dir
            ref_start_pdb = f"ranked_{i}.pdb"
            ref_start_ranking_json_file = f"ranking_debug.json"

            model_iteration_scores = []
            model_iteration_tmscores = []

            print(f"Start to refine {ref_start_pdb}")

            for num_iteration in range(self.max_iteration):
                os.chdir(cwd)
                current_work_dir = f"{model_outdir}/iteration{num_iteration + 1}"
                makedir_if_not_exists(current_work_dir)

                start_pdb = f"{current_work_dir}/start.pdb"
                start_msa = f"{current_work_dir}/start.a3m"
                start_ranking_json_file = f"{current_work_dir}/start_ranking.json"

                os.system(f"cp {current_ref_dir}/{ref_start_pdb} {start_pdb}")
                os.system(f"cp {current_ref_dir}/{ref_start_ranking_json_file} {start_ranking_json_file}")
                os.system(f"cp {current_ref_dir}/msas/final.a3m {start_msa}")

                ranking_json = json.loads(open(start_ranking_json_file).read())

                if num_iteration == 0:
                    ref_avg_lddt = ranking_json["plddts"][list(ranking_json["order"])[i]]
                else:
                    ref_avg_lddt = ranking_json["plddts"][list(ranking_json["order"])[0]]

                ref_tmscore = 0
                if os.path.exists(native_pdb):
                    ref_tmscore, _ = _cal_tmscore(self.params['tmscore_program'], start_pdb, native_pdb, current_work_dir + '/tmp')

                model_iteration_scores += [ref_avg_lddt]
                model_iteration_tmscores += [ref_tmscore]

                out_model_dir = f"{current_work_dir}/alphafold"
                if not _complete_result(out_model_dir):

                    foldseek_res = self.search_templates(inpdb=start_pdb, outdir=current_work_dir + '/foldseek')

                    if not self.check_and_rank_templates(foldseek_res, f"{current_work_dir}/structure_templates.csv", query_sequence):
                        print(f"Cannot find any templates in iteration {num_iteration + 1}")
                        break

                    self.generate_msa_from_templates(fasta_file=fasta_file,
                                                     template_file=f"{current_work_dir}/structure_templates.csv",
                                                     start_msa=start_msa,
                                                     outfile=f"{current_work_dir}/iteration{num_iteration + 1}.a3m")

                    out_template_dir = f"{current_work_dir}/template_pdbs"
                    makedir_if_not_exists(out_template_dir)
                    self.copy_atoms_and_unzip(template_csv=f"{current_work_dir}/structure_templates.csv",
                                              outdir=out_template_dir)

                    makedir_if_not_exists(out_model_dir)
                    cmd = f"python run_alphafold_custom_sim.py " \
                          f"--fasta_path {fasta_file} " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir']} " \
                          f"--custom_msa {current_work_dir}/iteration{num_iteration + 1}.a3m " \
                          f"--temp_struct_csv {current_work_dir}/structure_templates.csv " \
                          f"--struct_atom_dir {out_template_dir} " \
                          f"--output_dir {out_model_dir}"

                    try:
                        os.chdir(self.params['alphafold_program_dir'])
                        os.system(cmd)
                    except Exception as e:
                        print(e)

                new_ranking_json_file = f"{out_model_dir}/ranking_debug.json"
                new_ranking_json = json.loads(open(new_ranking_json_file).read())
                max_lddt_score = new_ranking_json["plddts"][list(new_ranking_json["order"])[0]]

                print(f'#########Iteration: {num_iteration + 1}#############')
                print(f"plddt before: {ref_avg_lddt}")
                print(f"plddt after: {max_lddt_score}")
                if max_lddt_score > ref_avg_lddt:
                    print("Continue to refine")
                    current_ref_dir = out_model_dir
                    ref_start_pdb = f"ranked_0.pdb"
                    ref_start_ranking_json_file = f"ranking_debug.json"
                    print('##################################################')
                    if num_iteration + 1 >= self.max_iteration:
                        print("Reach maximum iteration")
                        ranking_json = json.loads(open(out_model_dir + '/ranking_debug.json').read())
                        ref_avg_lddt = ranking_json["plddts"][list(ranking_json["order"])[0]]

                        ref_tmscore = 0
                        if os.path.exists(native_pdb):
                            ref_tmscore, _ = _cal_tmscore(self.params['tmscore_program'],
                                                          out_model_dir + '/' + ref_start_pdb, native_pdb,
                                                          current_work_dir + '/tmp')
                        model_iteration_scores += [ref_avg_lddt]
                        model_iteration_tmscores += [ref_tmscore]
                else:
                    # keep the models in iteration 1 even through the plddt score decreases
                    if num_iteration == 0:

                        ranking_json = json.loads(open(out_model_dir + '/ranking_debug.json').read())
                        ref_avg_lddt = ranking_json["plddts"][list(ranking_json["order"])[0]]

                        ref_tmscore = 0
                        if os.path.exists(native_pdb):
                            ref_tmscore, _ = _cal_tmscore(self.params['tmscore_program'],
                                                       out_model_dir + '/' + ref_start_pdb, native_pdb, current_work_dir + '/tmp')
                        model_iteration_scores += [ref_avg_lddt]
                        model_iteration_tmscores += [ref_tmscore]
                    break

            # model_iteration_scores += [max_lddt_score]

            if len(model_iteration_scores) > 0:
                iteration_result_all['targetname'] += [targetname]
                iteration_result_all['model'] += [i]
                iteration_result_all['start_lddt'] += [model_iteration_scores[0]]
                iteration_result_all['end_lddt'] += [model_iteration_scores[len(model_iteration_scores) - 1]]
                iteration_result_all['start_tmscore'] += [model_iteration_tmscores[0]]
                iteration_result_all['end_tmscore'] += [model_iteration_tmscores[len(model_iteration_tmscores) - 1]]

            while len(model_iteration_scores) <= self.max_iteration:
                model_iteration_scores += [0]

            while len(model_iteration_tmscores) <= self.max_iteration:
                model_iteration_tmscores += [0]

            iteration_scores[f'model{i + 1}'] = model_iteration_scores
            true_tm_scores[f'model{i + 1}'] = model_iteration_tmscores

        iteration_result_avg['start_lddt'] = [np.mean(np.array(iteration_result_all['start_lddt']))]
        iteration_result_avg['end_lddt'] = [np.mean(np.array(iteration_result_all['end_lddt']))]
        iteration_result_avg['start_tmscore'] = [np.mean(np.array(iteration_result_all['start_tmscore']))]
        iteration_result_avg['end_tmscore'] = [np.mean(np.array(iteration_result_all['end_tmscore']))]

        iteration_result_max['start_lddt'] = [np.max(np.array(iteration_result_all['start_lddt']))]
        iteration_result_max['end_lddt'] = [np.max(np.array(iteration_result_all['end_lddt']))]
        iteration_result_max['start_tmscore'] = [np.max(np.array(iteration_result_all['start_tmscore']))]
        iteration_result_max['end_tmscore'] = [np.max(np.array(iteration_result_all['end_tmscore']))]

        print(iteration_scores)
        df = pd.DataFrame(iteration_scores)
        df.to_csv(outdir + '/summary.csv')

        df = pd.DataFrame(true_tm_scores)
        df.to_csv(outdir + '/tmscores.csv')

        df = pd.DataFrame(iteration_result_avg)
        df.to_csv(outdir + '/iteration_result_avg.csv')

        df = pd.DataFrame(iteration_result_all)
        df.to_csv(outdir + '/iteration_result_all.csv')

        df = pd.DataFrame(iteration_result_max)
        df.to_csv(outdir + '/iteration_result_max.csv')

        os.chdir(cwd)

        return iteration_result_all, iteration_result_avg, iteration_result_max
