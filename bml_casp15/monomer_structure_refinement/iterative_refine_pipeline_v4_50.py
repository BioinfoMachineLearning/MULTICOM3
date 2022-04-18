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
from bml_casp15.complex_templates_search.sequence_based_pipeline import assess_hhsearch_hit, PrefilterError
from bml_casp15.complex_templates_search.parsers import TemplateHit


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


def _cal_tmscore(tmscore_program, inpdb, nativepdb, tmpdir):

    cwd = os.getcwd()

    makedir_if_not_exists(tmpdir)

    os.chdir(tmpdir)

    os.system(f"cp {inpdb} inpdb.pdb")
    os.system(f"cp {nativepdb} native.pdb")

    print(inpdb)
    print(nativepdb)
    inpdb = "inpdb.pdb"
    nativepdb = "native.pdb"

    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $3}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[2].rstrip('\n'))
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep GDT-score | awk '{print $3}' "
    tmscore_contents = os.popen(cmd).read().split('\n')
    gdtscore = float(tmscore_contents[0].rstrip('\n'))

    os.chdir(cwd)

    os.system("rm -rf " + tmpdir)

    return tmscore, gdtscore


def build_alignment_indices(sequence, start_index):
    indices_list = []
    counter = start_index
    for symbol in sequence:
        if symbol == '-':
            indices_list.append(-1)
        else:
            indices_list.append(counter)
            counter += 1
    return indices_list


class Monomer_iterative_refinement_pipeline:

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
        return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive_threshold=2000)

    def check_and_rank_templates(self, template_result, outfile, query_sequence):

        evalue_keep_indices = []
        for i in range(len(template_result['local_alignment'])):
            hit = TemplateHit(index=i,
                              name=template_result['local_alignment'].loc[i, 'target'].split('.')[0],
                              aligned_cols=int(template_result['local_alignment'].loc[i, 'alnlen']),
                              query=template_result['local_alignment'].loc[i, 'qaln'],
                              hit_sequence=template_result['local_alignment'].loc[i, 'taln'],
                              indices_query=build_alignment_indices(template_result['local_alignment'].loc[i, 'qaln'],
                                                                    template_result['local_alignment'].loc[
                                                                        i, 'qstart']),
                              indices_hit=build_alignment_indices(template_result['local_alignment'].loc[i, 'taln'],
                                                                  template_result['local_alignment'].loc[i, 'tstart']),
                              sum_probs=0.0)
            try:
                assess_hhsearch_hit(hit=hit, query_sequence=query_sequence)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            evalue_keep_indices += [i]

        tmscore_keep_indices = []
        for i in range(len(template_result['global_alignment'])):
            hit = TemplateHit(index=i,
                              name=template_result['global_alignment'].loc[i, 'target'].split('.')[0],
                              aligned_cols=int(template_result['global_alignment'].loc[i, 'alnlen']),
                              query=template_result['global_alignment'].loc[i, 'qaln'],
                              hit_sequence=template_result['global_alignment'].loc[i, 'taln'],
                              indices_query=build_alignment_indices(template_result['global_alignment'].loc[i, 'qaln'],
                                                                    template_result['global_alignment'].loc[
                                                                        i, 'qstart']),
                              indices_hit=build_alignment_indices(template_result['global_alignment'].loc[i, 'taln'],
                                                                  template_result['global_alignment'].loc[i, 'tstart']),
                              sum_probs=0.0)
            try:
                assess_hhsearch_hit(hit=hit, query_sequence=query_sequence)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            tmscore_keep_indices += [i]

        if len(evalue_keep_indices) == 0 and len(tmscore_keep_indices) == 0:
            return False

        evalue_thresholds = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
        tmscore_thresholds = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3]

        templates_sorted = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])

        evalue_af_indices = []
        evalue_pdb_indices = []
        tmscore_af_indices = []
        tmscore_pdb_indices = []
        for evalue_threshold, tmscore_threshold in zip(evalue_thresholds, tmscore_thresholds):
            evalue_af_indices = []
            evalue_pdb_indices = []
            for i in evalue_keep_indices:
                target = template_result['local_alignment'].loc[i, 'target']
                evalue = float(template_result['local_alignment'].loc[i, 'evalue'])
                if evalue < evalue_threshold:
                    if target.find('.atom.gz') > 0:
                        evalue_pdb_indices += [i]
                    else:
                        evalue_af_indices += [i]

            tmscore_af_indices = []
            tmscore_pdb_indices = []
            for i in tmscore_keep_indices:
                target = template_result['global_alignment'].loc[i, 'target']
                evalue = float(template_result['global_alignment'].loc[i, 'evalue'])
                if evalue > tmscore_threshold:
                    if target.find('.atom.gz') > 0:
                        tmscore_pdb_indices += [i]
                    else:
                        tmscore_af_indices += [i]

            if len(evalue_af_indices) + len(evalue_pdb_indices) \
                    + len(tmscore_af_indices) + len(tmscore_pdb_indices) >= self.max_template_count:
                break

        templates_sorted = copy.deepcopy(template_result['local_alignment'].iloc[evalue_pdb_indices])
        templates_sorted = templates_sorted.append(
            copy.deepcopy(template_result['global_alignment'].iloc[tmscore_pdb_indices]))
        templates_sorted = templates_sorted.append(
                copy.deepcopy(template_result['local_alignment'].iloc[evalue_af_indices]))
        templates_sorted = templates_sorted.append(
                copy.deepcopy(template_result['global_alignment'].iloc[tmscore_af_indices]))

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

    def search_single(self, fasta_path, pdb_path, pkl_path, msa_path, outdir):

        query_sequence = ""
        for line in open(fasta_path):
            line = line.rstrip('\n')
            if line.startswith('>'):
                continue
            else:
                query_sequence = line

        makedir_if_not_exists(outdir)

        cwd = os.getcwd()

        ref_start_pdb = pdb_path
        ref_start_pkl = pkl_path
        ref_start_msa = msa_path

        model_iteration_scores = []

        print(f"Start to refine {pdb_path}")

        for num_iteration in range(self.max_iteration):
            os.chdir(cwd)
            current_work_dir = f"{outdir}/iteration{num_iteration + 1}"
            makedir_if_not_exists(current_work_dir)

            start_pdb = f"{current_work_dir}/start.pdb"
            start_msa = f"{current_work_dir}/start.a3m"
            start_pkl = f"{current_work_dir}/start.pkl"

            os.system(f"cp {ref_start_pdb} {start_pdb}")
            os.system(f"cp {ref_start_msa} {start_msa}")
            os.system(f"cp {ref_start_pkl} {start_pkl}")

            with open(ref_start_pkl, 'rb') as f:
                ref_avg_lddt = np.mean(pickle.load(f)['plddt'])

            model_iteration_scores += [ref_avg_lddt]

            out_model_dir = f"{current_work_dir}/alphafold"
            if not _complete_result(out_model_dir):

                foldseek_res = self.search_templates(inpdb=start_pdb, outdir=current_work_dir + '/foldseek')

                if not self.check_and_rank_templates(foldseek_res, f"{current_work_dir}/structure_templates.csv", query_sequence):
                    print(f"Cannot find any templates in iteration {num_iteration + 1}")
                    break

                self.generate_msa_from_templates(fasta_file=fasta_path,
                                                 template_file=f"{current_work_dir}/structure_templates.csv",
                                                 start_msa=start_msa,
                                                 outfile=f"{current_work_dir}/iteration{num_iteration + 1}.a3m")

                out_template_dir = f"{current_work_dir}/template_pdbs"
                makedir_if_not_exists(out_template_dir)
                self.copy_atoms_and_unzip(template_csv=f"{current_work_dir}/structure_templates.csv",
                                          outdir=out_template_dir)

                makedir_if_not_exists(out_model_dir)
                cmd = f"python run_alphafold_custom_sim.py " \
                      f"--fasta_path {fasta_path} " \
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
                ref_start_pdb = f"{out_model_dir}/ranked_0.pdb"
                model_num = list(new_ranking_json["order"])[0].split('_')[1]
                ref_start_pkl = f"{out_model_dir}/result_model_{model_num}.pkl"
                ref_start_msa = f"{out_model_dir}/msas/monomer_final.a3m"
                print('##################################################')
                if num_iteration + 1 >= self.max_iteration:
                    print("Reach maximum iteration")
                    model_iteration_scores += [max_lddt_score]
            else:
                # keep the models in iteration 1 even through the plddt score decreases
                if num_iteration == 0:
                    ref_start_pdb = f"{out_model_dir}/ranked_0.pdb"
                    model_num = list(new_ranking_json["order"])[0].split('_')[1]
                    ref_start_pkl = f"{out_model_dir}/result_model_{model_num}.pkl"
                    ref_start_msa = f"{out_model_dir}/msas/monomer_final.a3m"
                    model_iteration_scores += [max_lddt_score]
                break

        while len(model_iteration_scores) <= self.max_iteration:
            model_iteration_scores += [0]

        print(model_iteration_scores)
        df = pd.DataFrame(model_iteration_scores)
        df.to_csv(outdir + '/summary.csv')

        final_model_dir = outdir + '/final'

        makedir_if_not_exists(final_model_dir)

        os.system(f"cp {ref_start_pdb} {final_model_dir}/final.pdb")
        os.system(f"cp {ref_start_pkl} {final_model_dir}/final.pkl")
        os.system(f"cp {ref_start_msa} {final_model_dir}/final.a3m")

        os.chdir(cwd)

        return final_model_dir

