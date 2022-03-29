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
from bml_casp15.complex_templates_search.sequence_based_pipeline import assess_hhsearch_hit
from bml_casp15.complex_templates_search.parsers import TemplateHit
from bml_casp15.tertiary_structure_generation.iterative_search_pipeline import build_alignment_indices, PrefilterError
from bml_casp15.quaternary_structure_generation.iterative_search_pipeline import *
from bml_casp15.monomer_alignment_generation.alignment import read_a3m


class Multimer_iterative_generation_pipeline_monomer:

    def __init__(self, params, max_template_count=50):

        self.params = params

        self.max_iteration = 5

        self.max_template_count = max_template_count

    def concatenate_msa_and_templates(self,
                                      chain_id_map,
                                      template_files,
                                      monomer_a3ms,
                                      outpath,
                                      template_path,
                                      rank_templates_by_monomers=True):

        prev_df = None
        for i, chain_id in enumerate(chain_id_map):
            templates = pd.read_csv(template_files[i], sep='\t')
            curr_df = create_template_df(templates)
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, how="inner", on='tpdbcode', suffixes=(str(i), str(i + 1)))

        min_evalues = []
        for i in range(len(prev_df)):
            evalues = []
            for j, chain_id in enumerate(chain_id_map):
                evalue = float(prev_df.loc[i, f'evalue{j + 1}'])
                evalues += [evalue]
            min_evalues += [np.min(np.array(evalues))]
        prev_df['min_evalue'] = min_evalues
        prev_df = prev_df.sort_values(by=['min_evalue'])

        keep_indices = []
        chain_template_msas = {}
        for chain_id in chain_id_map:
            chain_template_msas[chain_id] = {'desc': [chain_id_map[chain_id].description],
                                             'seq': [chain_id_map[chain_id].sequence]}

        print(prev_df)
        for i in range(len(prev_df)):
            template_infos = []
            for j, chain_id in enumerate(chain_id_map):
                template = prev_df.loc[i, f'template{j + 1}']
                qaln = prev_df.loc[i, f'aln_query{j + 1}']
                qstart = int(prev_df.loc[i, f'qstart{j + 1}'])
                qend = int(prev_df.loc[i, f'qend{j + 1}'])
                taln = prev_df.loc[i, f'aln_temp{j + 1}']
                tstart = prev_df.loc[i, f'tstart{j + 1}']
                tend = prev_df.loc[i, f'tend{j + 1}']
                evalue = float(prev_df.loc[i, f'evalue{j + 1}'])
                row_dict = dict(chainid=chain_id,
                                template=template,
                                tpdbcode=template[0:4],
                                aln_temp=taln,
                                tstart=tstart,
                                tend=tend,
                                aln_query=qaln,
                                qstart=qstart,
                                qend=qend,
                                evalue=evalue)
                template_infos += [row_dict]

            if not assess_complex_templates(chain_id_map, template_infos, template_path):
                continue

            if len(keep_indices) >= self.max_template_count:
                break

            keep_indices += [i]
            for j, chain_id in enumerate(chain_id_map):
                query_non_gaps = [res != '-' for res in prev_df.loc[i, f'aln_query{j + 1}']]
                out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, prev_df.loc[i, f'aln_temp{j + 1}']))
                aln_full = ['-'] * len(chain_id_map[chain_id].sequence)

                qstart = int(prev_df.loc[i, f'qstart{j + 1}'])
                qend = int(prev_df.loc[i, f'qend{j + 1}'])
                aln_full[qstart - 1:qend] = out_sequence
                taln_full_seq = ''.join(aln_full)
                chain_template_msas[chain_id]['desc'] += [prev_df.loc[i, f'template{j + 1}']]
                chain_template_msas[chain_id]['seq'] += [taln_full_seq]

        msa_out_path = outpath  # + '/msas'
        makedir_if_not_exists(msa_out_path)

        interact_template_count = len(keep_indices)
        num_monomer = len(template_files)
        max_msa_count = 50000
        num_msa_per_monomer = int((max_msa_count - interact_template_count) / num_monomer)
        print(num_msa_per_monomer)

        out_msas = []

        msa_per_monomer = {}
        for chain_id in chain_id_map:
            msa_per_monomer[chain_id] = {'desc': [], 'seq': []}

        for chain_idx, chain_id in enumerate(chain_id_map):
            fasta_chunks = (f">{chain_template_msas[chain_id]['desc'][i]}\n{chain_template_msas[chain_id]['seq'][i]}"
                            for i in range(len(chain_template_msas[chain_id]['desc'])))
            with open(msa_out_path + '/' + chain_id_map[chain_id].description + '.temp.interact', 'w') as fw:
                fw.write('\n'.join(fasta_chunks) + '\n')

            templates = pd.read_csv(template_files[chain_idx], sep='\t')
            seen_sequences = [chain_template_msas[chain_id]['seq'][i] for i in
                              range(len(chain_template_msas[chain_id]['desc']))]
            monomer_msas = {'desc': [], 'seq': []}
            print(templates)

            for j in range(len(templates)):
                if len(monomer_msas['seq']) > num_msa_per_monomer:
                    break
                query_non_gaps = [res != '-' for res in templates.loc[j, f'qaln']]
                out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, templates.loc[j, f'taln']))
                aln_full = ['-'] * len(chain_id_map[chain_id].sequence)

                qstart = int(templates.loc[j, f'qstart'])
                qend = int(templates.loc[j, f'qend'])
                aln_full[qstart - 1:qend] = out_sequence
                taln_full_seq = ''.join(aln_full)

                if taln_full_seq not in seen_sequences:
                    monomer_msas['desc'] += [templates.loc[j, f'target']]
                    monomer_msas['seq'] += [taln_full_seq]
                    seen_sequences += [taln_full_seq]

            print(f"monomer template alignment depth: {len(monomer_msas['seq'])}")

            if len(monomer_msas['seq']) < num_msa_per_monomer:
                seqs = None
                with open(monomer_a3ms[chain_idx]) as f:
                    seqs = read_a3m(f)
                for header in seqs:
                    if len(monomer_msas['seq']) > num_msa_per_monomer:
                        break
                    if header == chain_id_map[chain_id].description:
                        continue
                    if seqs[header] not in seen_sequences:
                        monomer_msas['desc'] += [header]
                        monomer_msas['seq'] += [seqs[header]]
                        seen_sequences += [seqs[header]]

            print(f"monomer msa alignment depth: {len(monomer_msas['seq'])}")

            msa_per_monomer[chain_id]['desc'] += monomer_msas['desc']
            msa_per_monomer[chain_id]['seq'] += monomer_msas['seq']

            print(f"final monomer msa alignment depth: {len(msa_per_monomer[chain_id]['seq'])}")

            for other_chain_id in chain_id_map:
                if other_chain_id == chain_id:
                    continue
                msa_per_monomer[other_chain_id]['desc'] += ['placeholder'] * len(monomer_msas['desc'])
                other_residue_num = len(chain_id_map[other_chain_id].sequence)
                msa_per_monomer[other_chain_id]['seq'] += [''.join(['-'] * other_residue_num)] * len(
                    monomer_msas['desc'])

        for chain_id in msa_per_monomer:
            fasta_chunks = (f">{msa_per_monomer[chain_id]['desc'][i]}\n{msa_per_monomer[chain_id]['seq'][i]}"
                            for i in range(len(msa_per_monomer[chain_id]['desc'])))
            with open(msa_out_path + '/' + chain_id_map[chain_id].description + '.temp.monomer', 'w') as fw:
                fw.write('\n'.join(fasta_chunks) + '\n')

            combine_a3ms([msa_out_path + '/' + chain_id_map[chain_id].description + '.temp.interact',
                          msa_out_path + '/' + chain_id_map[chain_id].description + '.temp.monomer'],
                         f"{msa_out_path}/{chain_id_map[chain_id].description}.iteration.a3m")
            out_msas += [f"{msa_out_path}/{chain_id_map[chain_id].description}.iteration.a3m"]

        interact_dict = {}
        msa_len = -1
        for i in range(0, len(out_msas)):
            msa_sequences, msa_descriptions = parse_fasta(out_msas[i])
            current_len = len(msa_descriptions)
            if msa_len == -1:
                msa_len = current_len
            elif current_len != msa_len:
                raise Exception(f"The length of each msas are not equal! {out_msas}")
            interact_dict[f'index_{i + 1}'] = [j for j in range(msa_len)]
            
        interact_df = pd.DataFrame(interact_dict)
        interact_csv = outpath + f'/interaction.iteration.csv'
        interact_df.to_csv(interact_csv)

        if rank_templates_by_monomers:
            top_template_files = []
            for template_file, chain_id in zip(template_files, chain_id_map):
                templates = pd.read_csv(template_file, sep='\t')
                keep_indices = []
                for i in range(len(templates)):
                    hit = TemplateHit(index=i,
                                      name=templates.loc[i, 'target'].split('.')[0],
                                      aligned_cols=int(templates.loc[i, 'alnlen']),
                                      query=templates.loc[i, 'qaln'],
                                      hit_sequence=templates.loc[i, 'taln'],
                                      indices_query=build_alignment_indices(templates.loc[i, 'qaln'],
                                                                            templates.loc[i, 'qstart']),
                                      indices_hit=build_alignment_indices(templates.loc[i, 'taln'],
                                                                          templates.loc[i, 'tstart']),
                                      sum_probs=0.0)
                    try:
                        assess_hhsearch_hit(hit=hit, query_sequence=chain_id_map[chain_id].sequence)
                    except PrefilterError as e:
                        msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                        print(msg)
                        continue
                    keep_indices += [i]
                templates_sorted = copy.deepcopy(templates.iloc[keep_indices])
                templates_sorted.drop(templates_sorted.filter(regex="Unnamed"), axis=1, inplace=True)
                templates_sorted.reset_index(inplace=True, drop=True)
                templates_sorted.to_csv(template_file + f'.top{self.max_template_count}', sep='\t')
                top_template_files += [template_file + f'.top{self.max_template_count}']
            return top_template_files, out_msas, interact_csv

        prev_df.iloc[keep_indices].to_csv(outpath + '/complex_templates.csv')
        return [outpath + '/complex_templates.csv'], out_msas, interact_csv

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

    def search(self, fasta_file, monomer_pdb_dir, outdir, native_pdb_dir=""):

        fasta_file = os.path.abspath(fasta_file)

        targetname = pathlib.Path(fasta_file).stem

        print(f"Processing {targetname}")

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        sequences, descriptions = parse_fasta(fasta_file)

        monomer_abs_dirs = {}
        chain_id_map = {}
        for chain_id, sequence, description in zip(PDB_CHAIN_IDS_UNRELAX, sequences, descriptions):
            chain_id_map[chain_id] = FastaChain(sequence=sequence, description=description)

            if not os.path.exists(monomer_pdb_dir + '/' + description):
                print(f"Cannot find monomer models for {description}: {monomer_pdb_dir}/{description}")
                return
            monomer_abs_dirs[chain_id] = os.path.abspath(monomer_pdb_dir + '/' + description)

        native_pdb = ""
        if os.path.exists(native_pdb_dir):
            native_pdb = outdir + '/' + '_'.join(
                [chain_id_map[chain_id].description for chain_id in chain_id_map]) + '.atom'
            combine_pdb(
                [native_pdb_dir + '/' + chain_id_map[chain_id].description + '.atom'
                 if os.path.exists(native_pdb_dir + '/' + chain_id_map[chain_id].description + '.atom')
                 else native_pdb_dir + '/' + chain_id_map[chain_id].description + '.pdb'
                 for chain_id in chain_id_map],
                native_pdb)

        cwd = os.getcwd()

        # outdir = f"{outdir}/{'_'.join(descriptions)}"

        makedir_if_not_exists(outdir)

        out_model_dir = outdir + '/alphafold'
        makedir_if_not_exists(out_model_dir)

        if not complete_result(out_model_dir):
            out_template_dir = outdir + '/templates'

            makedir_if_not_exists(out_template_dir)

            template_files = []
            alphafold_a3ms = []

            for chain_id in chain_id_map:

                monomer_work_dir = outdir + '/' + chain_id_map[chain_id].description

                makedir_if_not_exists(monomer_work_dir)

                chain_final_a3m = monomer_abs_dirs[chain_id] + '/msas/final.a3m'

                if not os.path.exists(chain_final_a3m):
                    raise Exception(f"Cannot find the final a3m in {monomer_abs_dirs[chain_id]}")

                os.system(f"cp {chain_final_a3m} {outdir}/{chain_id_map[chain_id].description}.alphafold.a3m")

                alphafold_a3ms += [f"{outdir}/{chain_id_map[chain_id].description}.alphafold.a3m"]

                chain_pdb = monomer_abs_dirs[chain_id] + '/ranked_0.pdb'

                os.system(f"cp {chain_pdb} {monomer_work_dir}/{chain_id_map[chain_id].description}.pdb")

                foldseek_res = search_templates_foldseek(
                    foldseek_program=self.params['foldseek_program'],
                    databases=[self.params['foldseek_pdb_database'], self.params['foldseek_af_database']],
                    inpdb=f"{monomer_work_dir}/{chain_id_map[chain_id].description}.pdb",
                    outdir=monomer_work_dir + '/foldseek')

                if not check_and_rank_templates(foldseek_res, f"{monomer_work_dir}/structure_templates.csv"):
                    print(
                        f"Cannot find any templates for {chain_id_map[chain_id].description} in iteration {num_iteration + 1}")
                    break

                template_files += [f"{monomer_work_dir}/structure_templates.csv"]

                self.copy_atoms_and_unzip(template_csv=f"{monomer_work_dir}/structure_templates.csv",
                                          outdir=out_template_dir)

            if len(template_files) != len(chain_id_map):
                return

            template_files, msa_files, msa_pair_file = self.concatenate_msa_and_templates(chain_id_map=chain_id_map,
                                                                                          template_files=template_files,
                                                                                          monomer_a3ms=alphafold_a3ms,
                                                                                          template_path=out_template_dir,
                                                                                          outpath=outdir)
            
            if len(template_files) == 1:
                cmd = f"python run_alphafold_multimer_custom_sim.py " \
                      f"--fasta_path {fasta_file} " \
                      f"--env_dir {self.params['alphafold_env_dir']} " \
                      f"--database_dir {self.params['alphafold_database_dir']} " \
                      f"--a3ms {','.join(msa_files)} " \
                      f"--msa_pair_file {msa_pair_file} " \
                      f"--temp_struct_csv {template_files[0]} " \
                      f"--struct_atom_dir {out_template_dir} " \
                      f"--output_dir {out_model_dir}"
            else:
                cmd = f"python run_alphafold_multimer_custom_sim.py " \
                      f"--fasta_path {fasta_file} " \
                      f"--env_dir {self.params['alphafold_env_dir']} " \
                      f"--database_dir {self.params['alphafold_database_dir']} " \
                      f"--a3ms {','.join(msa_files)} " \
                      f"--msa_pair_file {msa_pair_file} " \
                      f"--monomer_temp_csvs {','.join(template_files)} " \
                      f"--struct_atom_dir {out_template_dir} " \
                      f"--output_dir {out_model_dir}"

            try:
                os.chdir(self.params['alphafold_program_dir'])
                print(cmd)
                os.system(cmd)
            except Exception as e:
                print(e)

        targets = [targetname] * 5
        tmscores = [0] * 5
        tmaligns = [0] * 5
        if os.path.exists(out_model_dir + '/ranked_0.pdb'):
            for i in range(0, 5):
                inpdb = f"{out_model_dir}/ranked_{i}.pdb"
                if os.path.exists(native_pdb):
                    tmscores[i] = cal_tmscore(self.params['mmalign_program'],
                                              inpdb, native_pdb)
                    tmaligns[i] = cal_tmalign(self.params['tmalign_program'],
                                              inpdb, native_pdb,
                                              out_model_dir + '/tmp')

        df_all = {'targetname': targets, 'tmscore': tmscores, 'tmalign': tmaligns}
        df_max = {'targetname': [targetname],
                  'tmscore': [np.max(np.array(tmscores))],
                  'tmalign': [np.max(np.array(tmaligns))]}
        return df_all, df_max
