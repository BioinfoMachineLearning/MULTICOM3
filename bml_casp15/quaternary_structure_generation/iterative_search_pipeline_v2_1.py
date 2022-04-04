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
from bml_casp15.quaternary_structure_generation.util import *


def search_templates_foldseek(foldseek_program, databases, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_runner = Foldseek(binary_path=foldseek_program, databases=databases)
    return foldseek_runner.query(pdb=inpdb, outdir=outdir)


class Multimer_iterative_generation_pipeline:

    def __init__(self, params, max_template_count=50):

        self.params = params

        self.max_iteration = 5

        self.max_template_count = max_template_count

    def concatenate_msa_and_templates(self,
                                      chain_id_map,
                                      template_results,
                                      start_msa_path,
                                      outpath,
                                      template_path,
                                      iteration):

        prev_df = None
        for i, chain_id in enumerate(chain_id_map):
            templates = template_results[i]['all_alignment']
            curr_df = create_template_df(templates)
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, how="inner", on='tpdbcode', suffixes=(str(i), str(i + 1)))

        keep_indices = []
        chain_template_msas = {}
        for chain_id in chain_id_map:
            chain_template_msas[chain_id] = {'desc': [chain_id_map[chain_id].description],
                                             'seq': [chain_id_map[chain_id].sequence]}

        print(prev_df)

        seen_complex_seq = []
        seen_complex_seq += ["".join([chain_template_msas[chain_id]['seq'][0] for chain_id in chain_template_msas])]
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

            if not assess_complex_templates(chain_id_map, template_infos):
                continue

            monomer_template_seqs = {}
            for j, chain_id in enumerate(chain_id_map):
                query_non_gaps = [res != '-' for res in prev_df.loc[i, f'aln_query{j + 1}']]
                out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, prev_df.loc[i, f'aln_temp{j + 1}']))
                aln_full = ['-'] * len(chain_id_map[chain_id].sequence)
                qstart = int(prev_df.loc[i, f'qstart{j + 1}'])
                qend = int(prev_df.loc[i, f'qend{j + 1}'])
                aln_full[qstart - 1:qend] = out_sequence
                taln_full_seq = ''.join(aln_full)
                monomer_template_dict = {'desc': prev_df.loc[i, f'template{j + 1}'], 'seq': taln_full_seq}
                monomer_template_seqs[chain_id] = monomer_template_dict

            complex_template_seq = "".join(
                [monomer_template_seqs[chain_id]['seq'] for chain_id in monomer_template_seqs])
            if complex_template_seq not in seen_complex_seq:
                for chain_id in monomer_template_seqs:
                    chain_template_msas[chain_id]['desc'] += [monomer_template_seqs[chain_id]['desc']]
                    chain_template_msas[chain_id]['seq'] += [monomer_template_seqs[chain_id]['seq']]
                seen_complex_seq += [complex_template_seq]
                keep_indices += [i]

        msa_out_path = outpath  # + '/msas'
        makedir_if_not_exists(msa_out_path)

        out_msas = []
        for chain_id in chain_id_map:
            start_msa = start_msa_path + '/' + chain_id_map[chain_id].description + '.start.a3m'
            fasta_chunks = (f">{chain_template_msas[chain_id]['desc'][i]}\n{chain_template_msas[chain_id]['seq'][i]}"
                            for i in range(len(chain_template_msas[chain_id]['desc'])))
            with open(start_msa + '.temp', 'w') as fw:
                fw.write('\n'.join(fasta_chunks) + '\n')
            combine_a3ms([start_msa, f"{start_msa}.temp"],
                         f"{msa_out_path}/{chain_id_map[chain_id].description}.iteration{iteration}.a3m")
            out_msas += [f"{msa_out_path}/{chain_id_map[chain_id].description}.iteration{iteration}.a3m"]

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
        interact_csv = outpath + f'/interaction.iteration{iteration}.csv'
        interact_df.to_csv(interact_csv)

        top_template_files = []
        for template_result, chain_id in zip(template_results, chain_id_map):
            check_and_rank_monomer_templates_local_or_global(template_result=template_result,
                                                             outfile=f"{outpath}/{chain_id_map[chain_id].description}.top{self.max_template_count}",
                                                             query_sequence=chain_id_map[chain_id].sequence,
                                                             max_template_count=self.max_template_count)
            top_template_files += [f"{outpath}/{chain_id_map[chain_id].description}.top{self.max_template_count}"]
        return top_template_files, out_msas, interact_csv
        #
        # prev_df.iloc[keep_indices].to_csv(outpath + '/complex_templates.csv')
        # return [outpath + '/complex_templates.csv'], out_msas, interact_csv

    def copy_atoms_and_unzip(self, templates, outdir):
        os.chdir(outdir)
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'target']
            if template_pdb.find('.pdb.gz') > 0:
                os.system(f"cp {self.params['foldseek_af_database_dir']}/{template_pdb} {outdir}")
            else:
                os.system(f"cp {self.params['foldseek_pdb_database_dir']}/{template_pdb} {outdir}")
            os.system(f"gunzip -f {template_pdb}")

    def search(self, fasta_file, input_pdb_dir, outdir, native_pdb_dir=""):

        input_pdb_dir = os.path.abspath(input_pdb_dir)

        fasta_file = os.path.abspath(fasta_file)

        targetname = pathlib.Path(fasta_file).stem

        print(f"Processing {targetname}")

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        sequences, descriptions = parse_fasta(fasta_file)

        chain_id_map = {}
        for chain_id, sequence, description in zip(PDB_CHAIN_IDS_UNRELAX, sequences, descriptions):
            chain_id_map[chain_id] = FastaChain(sequence=sequence, description=description)

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

        iteration_scores = {}

        true_tm_scores = {}

        iteration_result_all = {'targetname': [],
                                'model': [],
                                'start_lddt': [],
                                'end_lddt': [],
                                'start_tmscore': [],
                                'end_tmscore': [],
                                'start_tmalign': [],
                                'end_tmalign': []
                                }

        iteration_result_avg = {'targetname': [targetname], 'start_lddt': [], 'end_lddt': [],
                                'start_tmscore': [], 'end_tmscore': [], 'start_tmalign': [],
                                'end_tmalign': []}

        iteration_result_max = {'targetname': [targetname], 'start_lddt': [], 'end_lddt': [],
                                'start_tmscore': [], 'end_tmscore': [], 'start_tmalign': [],
                                'end_tmalign': []}

        cwd = os.getcwd()

        for i in range(0, 5):
            model_outdir = f"{outdir}/ranked_{i}"
            makedir_if_not_exists(model_outdir)

            current_ref_dir = input_pdb_dir
            ref_start_pdb = f"ranked_{i}.pdb"
            ref_start_ranking_json_file = f"ranking_debug.json"

            model_iteration_scores = []
            model_iteration_tmscores = []
            model_iteration_tmaligns = []

            print(f"Start to refine {ref_start_pdb}")

            for num_iteration in range(self.max_iteration):
                os.chdir(cwd)
                current_work_dir = f"{model_outdir}/iteration{num_iteration + 1}"
                makedir_if_not_exists(current_work_dir)

                start_pdb = f"{current_work_dir}/start.pdb"
                start_msa_path = f"{current_work_dir}/start_msas"
                start_ranking_json_file = f"{current_work_dir}/start_ranking.json"

                os.system(f"cp {current_ref_dir}/{ref_start_pdb} {start_pdb}")
                os.system(f"cp {current_ref_dir}/{ref_start_ranking_json_file} {start_ranking_json_file}")
                if os.path.exists(start_msa_path):
                    os.system(f"rm -rf {start_msa_path}")
                makedir_if_not_exists(start_msa_path)

                for chain_id in chain_id_map:
                    os.system(f"cp {current_ref_dir}/msas/{chain_id_map[chain_id].description}.paired.a3m "
                              f"{start_msa_path}/{chain_id_map[chain_id].description}.start.a3m")

                ranking_json = json.loads(open(start_ranking_json_file).read())

                if num_iteration == 0:
                    ref_avg_lddt = ranking_json["iptm+ptm"][list(ranking_json["order"])[i]]
                else:
                    ref_avg_lddt = ranking_json["iptm+ptm"][list(ranking_json["order"])[0]]

                ref_tmscore = 0
                ref_tmalign = 0
                if os.path.exists(native_pdb):
                    ref_tmscore = cal_tmscore(self.params['mmalign_program'], start_pdb, native_pdb)
                    ref_tmalign = cal_tmalign(self.params['tmalign_program'], start_pdb, native_pdb,
                                              current_work_dir + '/tmp')

                model_iteration_scores += [ref_avg_lddt]
                model_iteration_tmscores += [ref_tmscore]
                model_iteration_tmaligns += [ref_tmalign]

                out_model_dir = f"{current_work_dir}/alphafold"

                if not complete_result(out_model_dir):

                    chain_pdbs = split_pdb(start_pdb, current_work_dir)

                    template_results = []

                    out_template_dir = current_work_dir + '/templates'
                    makedir_if_not_exists(out_template_dir)

                    for chain_id in chain_pdbs:

                        if chain_id not in chain_id_map:
                            raise Exception("Multimer fasta file and model doesn't match!")

                        monomer_work_dir = current_work_dir + '/' + chain_id_map[chain_id].description
                        makedir_if_not_exists(monomer_work_dir)
                        os.system(
                            f"mv {chain_pdbs[chain_id]} {monomer_work_dir}/{chain_id_map[chain_id].description}.pdb")
                        foldseek_res = search_templates_foldseek(
                            foldseek_program=self.params['foldseek_program'],
                            databases=[self.params['foldseek_pdb_database'], self.params['foldseek_af_database']],
                            inpdb=f"{monomer_work_dir}/{chain_id_map[chain_id].description}.pdb",
                            outdir=monomer_work_dir + '/foldseek')

                        if len(foldseek_res['all_alignment']) == 0:
                            print(
                                f"Cannot find any templates for {chain_id_map[chain_id].description} in iteration {num_iteration + 1}")
                            break

                        # os.system(f"cp {foldseek_res} {monomer_work_dir}/structure_templates.csv")

                        template_results += [foldseek_res]

                        self.copy_atoms_and_unzip(templates=foldseek_res['all_alignment'],
                                                  outdir=out_template_dir)

                    if len(template_results) != len(chain_id_map):
                        break

                    template_files, msa_files, msa_pair_file = self.concatenate_msa_and_templates(
                        chain_id_map=chain_id_map,
                        template_results=template_results,
                        start_msa_path=start_msa_path,
                        template_path=out_template_dir,
                        outpath=current_work_dir,
                        iteration=num_iteration + 1)

                    find_templates = True
                    for chain_id, template_file in zip(chain_id_map, template_files):
                        if len(pd.read_csv(template_file, sep='\t')) == 0:
                            print(
                                f"Cannot find any templates for {chain_id_map[chain_id].description} in iteration {num_iteration + 1}")
                            find_templates = False

                    if not find_templates:
                        break

                    makedir_if_not_exists(out_model_dir)

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

                new_ranking_json_file = f"{out_model_dir}/ranking_debug.json"
                new_ranking_json = json.loads(open(new_ranking_json_file).read())
                max_lddt_score = new_ranking_json["iptm+ptm"][list(new_ranking_json["order"])[0]]

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
                        ref_avg_lddt = ranking_json["iptm+ptm"][list(ranking_json["order"])[0]]

                        ref_tmscore = 0
                        if os.path.exists(native_pdb):
                            ref_tmscore = cal_tmscore(self.params['mmalign_program'],
                                                      out_model_dir + '/' + ref_start_pdb, native_pdb)
                            ref_tmalign = cal_tmalign(self.params['tmalign_program'],
                                                      out_model_dir + '/' + ref_start_pdb, native_pdb,
                                                      out_model_dir + '/tmp')
                        model_iteration_scores += [ref_avg_lddt]
                        model_iteration_tmscores += [ref_tmscore]
                        model_iteration_tmaligns += [ref_tmalign]
                else:
                    # keep the models in iteration 1 even through the plddt score decreases
                    if num_iteration == 0:
                        ranking_json = json.loads(open(out_model_dir + '/ranking_debug.json').read())
                        ref_avg_lddt = ranking_json["iptm+ptm"][list(ranking_json["order"])[0]]

                        ref_tmscore = 0
                        if os.path.exists(native_pdb):
                            ref_tmscore = cal_tmscore(self.params['mmalign_program'],
                                                      out_model_dir + '/' + ref_start_pdb, native_pdb)
                            ref_tmalign = cal_tmalign(self.params['tmalign_program'],
                                                      out_model_dir + '/' + ref_start_pdb, native_pdb,
                                                      out_model_dir + '/tmp')
                        model_iteration_scores += [ref_avg_lddt]
                        model_iteration_tmscores += [ref_tmscore]
                        model_iteration_tmaligns += [ref_tmalign]
                    break

            # model_iteration_scores += [max_lddt_score]

            if len(model_iteration_scores) > 0:
                iteration_result_all['targetname'] += [targetname]
                iteration_result_all['model'] += [i]
                iteration_result_all['start_lddt'] += [model_iteration_scores[0]]
                iteration_result_all['end_lddt'] += [model_iteration_scores[len(model_iteration_scores) - 1]]
                iteration_result_all['start_tmscore'] += [model_iteration_tmscores[0]]
                iteration_result_all['end_tmscore'] += [model_iteration_tmscores[len(model_iteration_tmscores) - 1]]
                iteration_result_all['start_tmalign'] += [model_iteration_tmaligns[0]]
                iteration_result_all['end_tmalign'] += [model_iteration_tmaligns[len(model_iteration_tmaligns) - 1]]

            while len(model_iteration_scores) <= self.max_iteration:
                model_iteration_scores += [0]

            while len(model_iteration_tmscores) <= self.max_iteration:
                model_iteration_tmscores += [0]

            while len(model_iteration_tmaligns) <= self.max_iteration:
                model_iteration_tmaligns += [0]

            iteration_scores[f'model{i + 1}'] = model_iteration_scores
            true_tm_scores[f'model{i + 1}'] = model_iteration_tmscores

        iteration_result_avg['start_lddt'] = [np.mean(np.array(iteration_result_all['start_lddt']))]
        iteration_result_avg['end_lddt'] = [np.mean(np.array(iteration_result_all['end_lddt']))]
        iteration_result_avg['start_tmscore'] = [np.mean(np.array(iteration_result_all['start_tmscore']))]
        iteration_result_avg['end_tmscore'] = [np.mean(np.array(iteration_result_all['end_tmscore']))]
        iteration_result_avg['start_tmalign'] = [np.mean(np.array(iteration_result_all['start_tmalign']))]
        iteration_result_avg['end_tmalign'] = [np.mean(np.array(iteration_result_all['end_tmalign']))]

        iteration_result_max['start_lddt'] = [np.max(np.array(iteration_result_all['start_lddt']))]
        iteration_result_max['end_lddt'] = [np.max(np.array(iteration_result_all['end_lddt']))]
        iteration_result_max['start_tmscore'] = [np.max(np.array(iteration_result_all['start_tmscore']))]
        iteration_result_max['end_tmscore'] = [np.max(np.array(iteration_result_all['end_tmscore']))]
        iteration_result_max['start_tmalign'] = [np.max(np.array(iteration_result_all['start_tmalign']))]
        iteration_result_max['end_tmalign'] = [np.max(np.array(iteration_result_all['end_tmalign']))]

        print(iteration_scores)
        df = pd.DataFrame(iteration_scores)
        df.to_csv(outdir + '/summary.csv')

        df = pd.DataFrame(true_tm_scores)
        df.to_csv(outdir + '/tmscores.csv')

        print(iteration_result_avg)
        df = pd.DataFrame(iteration_result_avg)
        df.to_csv(outdir + '/iteration_result_avg.csv')

        df = pd.DataFrame(iteration_result_all)
        df.to_csv(outdir + '/iteration_result_all.csv')

        df = pd.DataFrame(iteration_result_max)
        df.to_csv(outdir + '/iteration_result_max.csv')

        os.chdir(cwd)

        return iteration_result_all, iteration_result_avg, iteration_result_max
