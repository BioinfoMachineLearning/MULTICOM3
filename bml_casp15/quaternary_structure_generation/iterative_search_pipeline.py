import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
from bml_casp15.tool.foldseek import *
import pickle
import numpy as np


def _combine_a3ms(infiles, outfile):
    descriptions = []
    seqs = []
    for infile in infiles:
        for line in open(infile):
            line = line.rstrip('\n')
            if line.startswith('>'):
                descriptions += [line]
            else:
                seqs += [line]

    with open(outfile, 'w') as fw:
        for (desc, seq) in zip(descriptions, seqs):
            fw.write(f"{desc}\n{seq}\n")


def _convert_taln_seq_to_a3m(query_non_gaps, aln):
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, aln):
        if is_query_res_non_gap:
            yield sequence_res


def _complete_result(outputdir):
    complete = True
    for i in range(0, 5):
        model = f'{outputdir}/ranked_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def _cal_tmscore(mmalign_program, inpdb, nativepdb):
    cmd = mmalign_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $3}' "
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[1].rstrip('\n'))
    return tmscore, gdtscore


class Monomer_iterative_generation_pipeline:

    def __init__(self, params, max_template_count=50):

        self.params = params

        self.max_iteration = 5

        self.max_template_count = max_template_count

    def search_templates(self, inpdb, outdir):
        makedir_if_not_exists(outdir)
        foldseek_program = self.params['foldseek_program']
        foldseek_pdb_database = self.params['foldseek_pdb_database']
        foldseek_af_database = self.params['foldseek_af_database']
        foldseek_runner = Foldseek(binary_path=foldseek_program,
                                   databases=[foldseek_pdb_database, foldseek_af_database])
        return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive=True)

    def check_and_rank_templates(self, template_file, outfile):
        templates = pd.read_csv(template_file, sep='\t')
        sort_indices = []
        for i in range(len(templates)):
            target = templates.loc[i, 'target']
            evalue = float(templates.loc[i, 'evalue'])
            if target.find('.atom.gz') > 0 and evalue < 1e-10:
                sort_indices += [i]
        for i in range(len(templates)):
            if len(sort_indices) >= self.max_template_count:
                break
            if i in sort_indices:
                continue
            sort_indices += [i]
        if len(sort_indices) == 0:
            os.system(f"cp {template_file} {outfile}")
            return False

        templates_sorted = copy.deepcopy(templates.iloc[sort_indices])
        templates_sorted.drop(templates_sorted.filter(regex="Unnamed"), axis=1, inplace=True)
        templates_sorted.reset_index(inplace=True, drop=True)
        templates_sorted.to_csv(outfile, sep='\t')
        return True

    def generate_msa_from_templates(self, fasta_file, start_msa, template_file, outfile):
        targetname = None
        seq = None
        for line in open(fasta_file):
            line = line.rstrip('\n')
            if line.startswith('>'):
                targetname = line[1:]
            else:
                seq = line

        templates = pd.read_csv(template_file, sep='\t')

        alignments = {targetname: seq}
        seen_seq = []
        for i in range(len(templates)):
            target = templates.loc[i, 'target']
            qaln = templates.loc[i, 'qaln']
            qstart = int(templates.loc[i, 'qstart'])
            qend = int(templates.loc[i, 'qend'])
            taln = templates.loc[i, 'taln']
            tstart = templates.loc[i, 'tstart']
            tend = templates.loc[i, 'tend']

            query_non_gaps = [res != '-' for res in qaln]
            out_sequence = ''.join(_convert_taln_seq_to_a3m(query_non_gaps, taln))

            aln_full = ['-'] * len(seq)
            aln_full[qstart - 1:qend] = out_sequence
            taln_full_seq = ''.join(aln_full)
            if taln_full_seq in seen_seq:
                continue
            alignments[target] = taln_full_seq
            seen_seq += [taln_full_seq]

        fasta_chunks = (f">{k}\n{alignments[k]}" for k in alignments)

        with open(outfile + '.temp', 'w') as fw:
            fw.write('\n'.join(fasta_chunks) + '\n')

        _combine_a3ms([start_msa, f"{outfile}.temp"], f"{outfile}.comb")

        cmd = f"{self.params['hhfilter_program']} -diff 50000 -i {outfile}.comb -o {outfile} -id 90"

        os.system(cmd)

    def copy_atoms_and_unzip(self, template_csv, outdir):
        os.chdir(outdir)
        templates = pd.read_csv(template_csv, sep='\t')
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'target']
            if template_pdb.find('.pdb.gz') > 0:
                os.system(f"cp {self.params['foldseek_af_database_dir']}/{template_pdb} {outdir}")
            else:
                os.system(f"cp {self.params['foldseek_pdb_database_dir']}/{template_pdb} {outdir}")
            os.system(f"gunzip -f {template_pdb}")

    def split_pdb(complex_pdb, outdir):
        pre_chain = None
        fw = None
        for line in open(complex_pdb, 'r').readlines():
            if not line.startswith('ATOM'):
                continue
            chain_name = line[21]
            if pre_chain is None:
                pre_chain = chain_name
                fw = open(outdir + '/' + chain_name + '.pdb', 'w')
                fw.write(line)
            elif chain_name == pre_chain:
                fw.write(line)
            else:
                fw.close()
                i = i + 1
                fw = open(outdir + '/' + chain_name + '.pdb', 'w')
                fw.write(line)
                pre_chain = chain_name
        fw.close()

    def search(self, fasta_file, input_pdb_dir, outdir, native_pdb=""):

        input_pdb_dir = os.path.abspath(input_pdb_dir)

        fasta_file = os.path.abspath(fasta_file)

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

        iteration_result_avg = {'targetname': [targetname], 'start_lddt': [], 'end_lddt': [], 'start_tmscore': [], 'end_tmscore': []}

        cwd = os.getcwd()

        for i in range(0, 5):
            model_outdir = f"{outdir}/ranked_{i}"
            makedir_if_not_exists(model_outdir)

            current_ref_dir = input_pdb_dir
            ref_start_pdb = f"ranked_{i}.pdb"
            ref_start_pkl = f"result_model_{i + 1}.pkl"

            model_iteration_scores = []
            model_iteration_tmscores = []

            print(f"Start to refine {ref_start_pdb}")

            for num_iteration in range(self.max_iteration):
                os.chdir(cwd)
                current_work_dir = f"{model_outdir}/iteration{num_iteration + 1}"
                makedir_if_not_exists(current_work_dir)

                start_pdb = f"{current_work_dir}/start.pdb"
                start_msa = f"{current_work_dir}/start.a3m"
                start_pkl = f"{current_work_dir}/start.pkl"

                os.system(f"cp {current_ref_dir}/{ref_start_pdb} {start_pdb}")
                os.system(f"cp {current_ref_dir}/{ref_start_pkl} {start_pkl}")
                os.system(f"cp {current_ref_dir}/msas/final.a3m {start_msa}")
                ref_avg_lddt = 0
                with open(start_pkl, 'rb') as f:
                    prediction_result = pickle.load(f)
                    ref_avg_lddt = np.mean(prediction_result['plddt'])

                ref_tmscore = 0
                if os.path.exists(native_pdb):
                    ref_tmscore, _ = _cal_tmscore(self.params['mmalign_program'], start_pdb, native_pdb)

                model_iteration_scores += [ref_avg_lddt]
                model_iteration_tmscores += [ref_tmscore]

                out_model_dir = f"{current_work_dir}/alphafold"
                if not _complete_result(out_model_dir):


                    foldseek_res = self.search_templates(inpdb=start_pdb, outdir=current_work_dir + '/foldseek')

                    if not self.check_and_rank_templates(foldseek_res, f"{current_work_dir}/structure_templates.csv"):
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

                max_lddt_score = 0
                max_index = -1
                for j in range(0, 5):
                    new_pkl = f"{out_model_dir}/result_model_{j + 1}.pkl"
                    with open(new_pkl, 'rb') as f:
                        new_prediction_result = pickle.load(f)
                        new_avg_lddt = np.mean(new_prediction_result['plddt'])
                        if new_avg_lddt > max_lddt_score:
                            max_lddt_score = new_avg_lddt
                            max_index = j

                print(f'#########Iteration: {num_iteration + 1}#############')
                print(f"plddt before: {ref_avg_lddt}")
                print(f"plddt after: {max_lddt_score}")
                if max_lddt_score > ref_avg_lddt:
                    print("Continue to refine")
                    current_ref_dir = out_model_dir
                    ref_start_pdb = f"ranked_{max_index}.pdb"
                    ref_start_pkl = f"result_model_{max_index + 1}.pkl"
                    print('##################################################')
                else:
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

        print(iteration_scores)
        df = pd.DataFrame(iteration_scores)
        df.to_csv(outdir + '/summary.csv')

        df = pd.DataFrame(true_tm_scores)
        df.to_csv(outdir + '/tmscores.csv')

        os.chdir(cwd)

        return iteration_result_all, iteration_result_avg
