import copy
import os
import sys
import time, json
from multicom3.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
from multicom3.tool.foldseek import *
import pickle
import numpy as np
from multicom3.monomer_templates_concatenation import parsers
from multicom3.monomer_structure_refinement.util import *
from multicom3.common.protein import complete_result


# Prefilter exceptions.
class PrefilterError(Exception):
    """A base class for template prefilter exceptions."""


class DateError(PrefilterError):
    """An error indicating that the hit date was after the max allowed date."""


class AlignRatioError(PrefilterError):
    """An error indicating that the hit align ratio to the query was too small."""


class DuplicateError(PrefilterError):
    """An error indicating that the hit was an exact subsequence of the query."""


class LengthError(PrefilterError):
    """An error indicating that the hit was too short."""

def assess_foldseek_hit(
        hit: parsers.TemplateHit,
        query_sequence: str,
        max_subsequence_ratio: float = 0.95,
        min_align_ratio: float = 0.1) -> bool:

    aligned_cols = hit.aligned_cols
    align_ratio = aligned_cols / len(query_sequence)

    template_sequence = hit.hit_sequence.replace('-', '')
    length_ratio = float(len(template_sequence)) / len(query_sequence)

    duplicate = (template_sequence in query_sequence and
                 length_ratio > max_subsequence_ratio)

    if align_ratio <= min_align_ratio:
        raise AlignRatioError('Proportion of residues aligned to query too small. '
                              f'Align ratio: {align_ratio}.')

    if duplicate:
        raise DuplicateError('Template is an exact subsequence of query with large '
                             f'coverage. Length ratio: {length_ratio}.')

    if len(template_sequence) < 10:
        raise LengthError(f'Template too short. Length: {len(template_sequence)}.')

    return True

class Monomer_iterative_refinement_pipeline:

    def __init__(self, params, max_template_count=50):

        self.params = params

        self.max_iteration = 5

        self.max_template_count = max_template_count

        release_date_df = pd.read_csv(params['pdb_release_date_file'])
        self._release_dates = dict(zip(release_date_df['pdbcode'], release_date_df['release_date']))
        self._max_template_date = datetime.datetime.strptime(params['max_template_date'], '%Y-%m-%d')

    def search_templates(self, inpdb, outdir):
        makedir_if_not_exists(outdir)
        foldseek_program = self.params['foldseek_program']
        foldseek_pdb_database = self.params['foldseek_pdb_database']
        foldseek_af_database = self.params['foldseek_af_database']
        foldseek_runner = Foldseek(binary_path=foldseek_program, pdb_database=foldseek_pdb_database,
                                   max_template_date=self._max_template_date, release_dates=self._release_dates,
                                   other_databases=[foldseek_af_database])
        return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive_threshold=2000)

    def check_and_rank_templates(self, template_result, outfile, query_sequence):

        evalue_keep_indices = []
        for i in range(len(template_result['local_alignment'])):
            hit = parsers.TemplateHit(index=i,
                              name=template_result['local_alignment'].loc[i, 'target'].split('.')[0],
                              aligned_cols=int(template_result['local_alignment'].loc[i, 'alnlen']),
                              query=template_result['local_alignment'].loc[i, 'qaln'],
                              hit_sequence=template_result['local_alignment'].loc[i, 'taln'],
                              indices_query=build_alignment_indices(template_result['local_alignment'].loc[i, 'qaln'],
                                                                    template_result['local_alignment'].loc[i, 'qstart']),
                              indices_hit=build_alignment_indices(template_result['local_alignment'].loc[i, 'taln'],
                                                                  template_result['local_alignment'].loc[i, 'tstart']),
                              sum_probs=0.0)
            try:
                assess_foldseek_hit(hit=hit, query_sequence=query_sequence)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            evalue_keep_indices += [i]

        tmscore_keep_indices = []
        for i in range(len(template_result['global_alignment'])):
            hit = parsers.TemplateHit(index=i,
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
                assess_foldseek_hit(hit=hit, query_sequence=query_sequence)
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
            out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, taln))

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

        combine_a3ms([start_msa, f"{outfile}.temp"], f"{outfile}.comb")

        cmd = f"{self.params['hhfilter_program']} -diff 50000 -i {outfile}.comb -o {outfile} -id 90"

        os.system(cmd)

    def copy_atoms_and_unzip(self, template_csv, outdir):
        os.chdir(outdir)
        templates = pd.read_csv(template_csv, sep='\t')
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'target']
            if template_pdb.find('.pdb.gz') > 0:
                template_path = os.path.join(self.params['foldseek_af_database_dir'], template_pdb)
                os.system(f"cp {template_path} {outdir}")
            else:
                template_path = os.path.join(self.params['foldseek_pdb_database_dir'], template_pdb)
                os.system(f"cp {template_path} {outdir}")
            os.system(f"gunzip -f {template_pdb}")

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
            current_work_dir = os.path.join(outdir, f"iteration{num_iteration + 1}")
            makedir_if_not_exists(current_work_dir)

            start_pdb = os.path.join(current_work_dir, "start.pdb")
            start_msa = os.path.join(current_work_dir, "start.a3m")
            start_pkl = os.path.join(current_work_dir, "start.pkl")

            os.system(f"cp {ref_start_pdb} {start_pdb}")
            os.system(f"cp {ref_start_msa} {start_msa}")
            os.system(f"cp {ref_start_pkl} {start_pkl}")

            with open(ref_start_pkl, 'rb') as f:
                ref_avg_lddt = np.mean(pickle.load(f)['plddt'])

            model_iteration_scores += [ref_avg_lddt]

            out_model_dir = os.path.join(current_work_dir, "alphafold")
            if not complete_result(out_model_dir, 5 * int(self.params['num_monomer_predictions_per_model'])):

                foldseek_res = self.search_templates(inpdb=start_pdb, outdir=os.path.join(current_work_dir, 'foldseek'))

                if not self.check_and_rank_templates(foldseek_res, os.path.join(current_work_dir, "structure_templates.csv"),
                                                     query_sequence):
                    print(f"Cannot find any templates in iteration {num_iteration + 1}")
                    break

                self.generate_msa_from_templates(fasta_file=fasta_path,
                                                 template_file=os.path.join(current_work_dir, "structure_templates.csv"),
                                                 start_msa=start_msa,
                                                 outfile=os.path.join(current_work_dir, f"iteration{num_iteration + 1}.a3m"))

                out_template_dir = os.path.join(current_work_dir, "template_pdbs")
                makedir_if_not_exists(out_template_dir)
                self.copy_atoms_and_unzip(template_csv=os.path.join(current_work_dir, "structure_templates.csv"),
                                          outdir=out_template_dir)

                makedir_if_not_exists(out_model_dir)
                custom_msa = os.path.join(current_work_dir, f"iteration{num_iteration + 1}.a3m")
                temp_struct_csv = os.path.join(current_work_dir, "structure_templates.csv")
                cmd = f"python {self.params['alphafold_program']} " \
                      f"--fasta_path={fasta_path} " \
                      f"--env_dir={self.params['alphafold_env_dir']} " \
                      f"--database_dir={self.params['alphafold_database_dir']} " \
                      f"--monomer_num_ensemble={self.params['monomer_num_ensemble']} " \
                      f"--monomer_num_recycle={self.params['monomer_num_recycle']} " \
                      f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                      f"--model_preset={self.params['monomer_model_preset']} " \
                      f"--benchmark={self.params['alphafold_benchmark']} " \
                      f"--use_gpu_relax={self.params['use_gpu_relax']} " \
                      f"--models_to_relax={self.params['models_to_relax']} " \
                      f"--max_template_date={self.params['max_template_date']} " \
                      f"--custom_msa={custom_msa} " \
                      f"--temp_struct_csv={temp_struct_csv} " \
                      f"--struct_atom_dir={out_template_dir} " \
                      f"--output_dir={out_model_dir}"

                try:
                    os.chdir(self.params['alphafold_program_dir'])
                    os.system(cmd)
                except Exception as e:
                    print(e)

            new_ranking_json_file = os.path.join(out_model_dir, "ranking_debug.json")
            new_ranking_json = json.loads(open(new_ranking_json_file).read())
            max_lddt_score = new_ranking_json["plddts"][list(new_ranking_json["order"])[0]]

            print(f'#########Iteration: {num_iteration + 1}#############')
            print(f"plddt before: {ref_avg_lddt}")
            print(f"plddt after: {max_lddt_score}")
            if max_lddt_score > ref_avg_lddt:
                print("Continue to refine")
                ref_start_pdb = os.path.join(out_model_dir, "ranked_0.pdb")
                model_name = list(new_ranking_json["order"])[0]
                ref_start_pkl = os.path.join(out_model_dir, f"result_{model_name}.pkl")
                ref_start_msa = os.path.join(out_model_dir, 'msas', "monomer_final.a3m")
                print('##################################################')
                if num_iteration + 1 >= self.max_iteration:
                    print("Reach maximum iteration")
                    model_iteration_scores += [max_lddt_score]
            else:
                # keep the models in iteration 1 even through the plddt score decreases
                if num_iteration == 0:
                    ref_start_pdb = os.path.join(out_model_dir, "ranked_0.pdb")
                    model_name = list(new_ranking_json["order"])[0]
                    ref_start_pkl = os.path.join(out_model_dir, f"result_{model_name}.pkl")
                    ref_start_msa = os.path.join(out_model_dir, 'msas', "monomer_final.a3m")
                    model_iteration_scores += [max_lddt_score]
                break

        while len(model_iteration_scores) <= self.max_iteration:
            model_iteration_scores += [0]

        print(model_iteration_scores)
        df = pd.DataFrame(model_iteration_scores)
        df.to_csv(os.path.join(outdir, 'summary.csv'))

        final_model_dir = os.path.join(outdir, 'final')

        makedir_if_not_exists(final_model_dir)

        os.system("cp " + ref_start_pdb + " " + os.path.join(final_model_dir, "final.pdb"))
        os.system("cp " + ref_start_pkl + " " + os.path.join(final_model_dir, "final.pkl"))
        os.system("cp " + ref_start_msa + " " + os.path.join(final_model_dir, "final.a3m"))
        os.chdir(cwd)

        return final_model_dir
