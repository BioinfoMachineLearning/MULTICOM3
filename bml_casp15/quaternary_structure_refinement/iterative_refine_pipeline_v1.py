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
from bml_casp15.quaternary_structure_refinement.util import *


def search_templates_foldseek(foldseek_program, databases, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_runner = Foldseek(binary_path=foldseek_program, databases=databases)
    return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive_threshold=2000, maxseq=300)


class Multimer_iterative_refinement_pipeline:

    def __init__(self, params, max_template_count=50):

        self.params = params

        self.max_iteration = 5

        self.max_template_count = max_template_count

    def concatenate_msa_and_templates(self,
                                      chain_id_map,
                                      template_results,
                                      start_msa_path,
                                      outpath,
                                      iteration):

        chain_template_msas = {}
        evalue_thresholds = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
        tmscore_thresholds = [0.8, 0.7, 0.6, 0.5, 0.4]
        complex_templates_df_filtered = None
        for evalue_threshold, tmscore_threshold in zip(evalue_thresholds, tmscore_thresholds):
            chain_template_msas = {}
            for chain_id in chain_id_map:
                chain_template_msas[chain_id] = {'desc': [chain_id_map[chain_id].description],
                                                 'seq': [chain_id_map[chain_id].sequence]}

            complex_templates_df = None
            for chain_idx, (chain_id, template_result) in enumerate(zip(chain_id_map, template_results)):
                evalue_keep_indices = []
                for i in range(len(template_result['local_alignment'])):
                    if template_result['local_alignment'].loc[i, 'evalue'] < evalue_threshold:
                        evalue_keep_indices += [i]

                tmscore_keep_indices = []
                for i in range(len(template_result['global_alignment'])):
                    if template_result['global_alignment'].loc[i, 'evalue'] > tmscore_threshold:
                        tmscore_keep_indices += [i]

                templates_filtered = copy.deepcopy(template_result['local_alignment'].iloc[evalue_keep_indices])
                templates_filtered = templates_filtered.append(
                    copy.deepcopy(template_result['global_alignment'].iloc[tmscore_keep_indices]))
                templates_filtered.drop(templates_filtered.filter(regex="Unnamed"), axis=1, inplace=True)
                templates_filtered.reset_index(inplace=True, drop=True)

                curr_df = create_template_df(templates_filtered)
                # print(curr_df)
                if complex_templates_df is None:
                    complex_templates_df = curr_df
                else:
                    complex_templates_df = complex_templates_df.merge(curr_df, how="inner", on='tpdbcode',
                                                                      suffixes=(str(chain_idx), str(chain_idx + 1)))

            keep_indices = []
            seen_complex_seq = []
            seen_complex_seq += ["".join([chain_template_msas[chain_id]['seq'][0] for chain_id in chain_template_msas])]
            for i in range(len(complex_templates_df)):
                if len(keep_indices) > self.max_template_count:
                    break
                template_infos = []
                for j, chain_id in enumerate(chain_id_map):
                    template = complex_templates_df.loc[i, f'template{j + 1}']
                    qaln = complex_templates_df.loc[i, f'aln_query{j + 1}']
                    qstart = int(complex_templates_df.loc[i, f'qstart{j + 1}'])
                    qend = int(complex_templates_df.loc[i, f'qend{j + 1}'])
                    taln = complex_templates_df.loc[i, f'aln_temp{j + 1}']
                    tstart = complex_templates_df.loc[i, f'tstart{j + 1}']
                    tend = complex_templates_df.loc[i, f'tend{j + 1}']
                    evalue = float(complex_templates_df.loc[i, f'evalue{j + 1}'])
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
                    query_non_gaps = [res != '-' for res in complex_templates_df.loc[i, f'aln_query{j + 1}']]
                    out_sequence = ''.join(
                        convert_taln_seq_to_a3m(query_non_gaps, complex_templates_df.loc[i, f'aln_temp{j + 1}']))
                    aln_full = ['-'] * len(chain_id_map[chain_id].sequence)
                    qstart = int(complex_templates_df.loc[i, f'qstart{j + 1}'])
                    qend = int(complex_templates_df.loc[i, f'qend{j + 1}'])
                    aln_full[qstart - 1:qend] = out_sequence
                    taln_full_seq = ''.join(aln_full)
                    monomer_template_dict = {'desc': complex_templates_df.loc[i, f'template{j + 1}'],
                                             'seq': taln_full_seq}
                    monomer_template_seqs[chain_id] = monomer_template_dict

                complex_template_seq = "".join(
                    [monomer_template_seqs[chain_id]['seq'] for chain_id in monomer_template_seqs])
                if complex_template_seq not in seen_complex_seq:
                    for chainid in monomer_template_seqs:
                        chain_template_msas[chainid]['desc'] += [monomer_template_seqs[chainid]['desc']]
                        chain_template_msas[chainid]['seq'] += [monomer_template_seqs[chainid]['seq']]
                    seen_complex_seq += [complex_template_seq]
                    keep_indices += [i]

            complex_templates_df_filtered = copy.deepcopy(complex_templates_df.iloc[keep_indices])
            complex_templates_df_filtered.drop(complex_templates_df_filtered.filter(regex="Unnamed"), axis=1,
                                               inplace=True)
            complex_templates_df_filtered.reset_index(inplace=True, drop=True)

            if len(complex_templates_df_filtered) > self.max_template_count:
                break

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
            check_and_rank_monomer_templates_local_and_global(template_result=template_result,
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

        finaldir = outdir + '/final'
        makedir_if_not_exists(finaldir)

        for i in range(0, 5):
            model_outdir = f"{outdir}/ranked_{i}"
            makedir_if_not_exists(model_outdir)

            current_ref_dir = input_pdb_dir
            ref_start_pdb = f"ranked_{i}.pdb"
            ref_start_ranking_json_file = f"ranking_debug.json"

            new_ranking_json = json.loads(open(current_ref_dir + '/' + ref_start_ranking_json_file).read())
            model_num = list(new_ranking_json["order"])[i].split('_')[1]
            ref_start_pkl = f"result_model_{model_num}_multimer.pkl"

            model_iteration_scores = []
            model_iteration_tmscores = []
            model_iteration_tmaligns = []

            print(f"Start to refine {ref_start_pdb}")

            os.system(f"cp {current_ref_dir}/{ref_start_pdb} {finaldir}/ranked_{i}_start.pdb")
            os.system(f"cp {current_ref_dir}/{ref_start_pkl} {finaldir}/ranked_{i}_start.pkl")

            for num_iteration in range(self.max_iteration):
                os.chdir(cwd)
                current_work_dir = f"{model_outdir}/iteration{num_iteration + 1}"
                makedir_if_not_exists(current_work_dir)

                start_pdb = f"{current_work_dir}/start.pdb"
                start_msa_path = f"{current_work_dir}/start_msas"
                start_pkl = f"{current_work_dir}/start.pkl"

                os.system(f"cp {current_ref_dir}/{ref_start_pdb} {start_pdb}")
                os.system(f"cp {current_ref_dir}/{ref_start_pkl} {start_pkl}")
                if os.path.exists(start_msa_path):
                    os.system(f"rm -rf {start_msa_path}")
                makedir_if_not_exists(start_msa_path)

                for chain_id in chain_id_map:
                    os.system(f"cp {current_ref_dir}/msas/{chain_id_map[chain_id].description}.paired.a3m "
                              f"{start_msa_path}/{chain_id_map[chain_id].description}.start.a3m")

                with open(start_pkl, 'rb') as f:
                    ref_avg_lddt = np.mean(pickle.load(f)['ranking_confidence'])

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
                    model_num = list(new_ranking_json["order"])[0].split('_')[1]
                    ref_start_pkl = f"result_model_{model_num}_multimer.pkl"
                    print('##################################################')
                    if num_iteration + 1 >= self.max_iteration:
                        print("Reach maximum iteration")
                        ref_tmscore = 0
                        if os.path.exists(native_pdb):
                            ref_tmscore = cal_tmscore(self.params['mmalign_program'],
                                                      out_model_dir + '/' + ref_start_pdb, native_pdb)
                            ref_tmalign = cal_tmalign(self.params['tmalign_program'],
                                                      out_model_dir + '/' + ref_start_pdb, native_pdb,
                                                      out_model_dir + '/tmp')
                        model_iteration_scores += [max_lddt_score]
                        model_iteration_tmscores += [ref_tmscore]
                        model_iteration_tmaligns += [ref_tmalign]
                else:
                    # keep the models in iteration 1 even through the plddt score decreases
                    if num_iteration == 0:
                        current_ref_dir = out_model_dir
                        ref_start_pdb = f"ranked_0.pdb"
                        model_num = list(new_ranking_json["order"])[0].split('_')[1]
                        ref_start_pkl = f"result_model_{model_num}_multimer.pkl"
                        ref_tmscore = 0
                        if os.path.exists(native_pdb):
                            ref_tmscore = cal_tmscore(self.params['mmalign_program'],
                                                      out_model_dir + '/' + ref_start_pdb, native_pdb)
                            ref_tmalign = cal_tmalign(self.params['tmalign_program'],
                                                      out_model_dir + '/' + ref_start_pdb, native_pdb,
                                                      out_model_dir + '/tmp')
                        model_iteration_scores += [max_lddt_score]
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

            os.system(f"cp {current_ref_dir}/{ref_start_pdb} {finaldir}/ranked_{i}_ref.pdb")
            os.system(f"cp {current_ref_dir}/{ref_start_pkl} {finaldir}/ranked_{i}_ref.pkl")

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

        pdbs = []
        plddts = []
        for pkl in os.listdir(finaldir):
            if pkl.find('.pkl') > 0:
                with open(finaldir + '/' + pkl, 'rb') as f:
                    prediction_result = pickle.load(f)
                    pdbs += [pkl.replace('.pkl', '.pdb')]
                    plddts += [prediction_result['ranking_confidence']]

        df = pd.DataFrame({'model': pdbs, 'plddt': plddts})
        df = df.sort_values(by=['plddt'], ascending=False)
        df.reset_index(inplace=True, drop=True)
        df.to_csv(outdir + '/final_ranking_v1.csv')

        v1_scores_all = []
        v1_scores_avg = 0
        v1_scores_max = 0
        for i in range(5):
            pdb_name = df.loc[i, 'model']
            ref_tmscore = cal_tmscore(self.params['mmalign_program'], finaldir + '/' + pdb_name, native_pdb)
            v1_scores_all += [ref_tmscore]

        v1_scores_avg = np.mean(np.array(v1_scores_all))
        v1_scores_max = np.max(np.array(v1_scores_all))

        v2_scores_all = []
        v2_scores_avg = 0
        v2_scores_max = 0
        for i in range(5):
            start_pdb = f"{finaldir}/ranked_{i}_start.pdb"
            start_pkl = f"{finaldir}/ranked_{i}_start.pkl"

            with open(start_pkl, 'rb') as f:
                plddt_start = pickle.load(f)['ranking_confidence']

            refine_pdb = f"{finaldir}/ranked_{i}_ref.pdb"
            refine_pkl = f"{finaldir}/ranked_{i}_ref.pkl"

            with open(refine_pkl, 'rb') as f:
                plddt_ref = pickle.load(f)['ranking_confidence']

            select_pdb = refine_pdb
            if plddt_start > plddt_ref:
                select_pdb = start_pdb

            ref_tmscore = cal_tmscore(self.params['mmalign_program'], select_pdb, native_pdb)
            v2_scores_all += [ref_tmscore]

        v2_scores_avg = np.mean(np.array(v2_scores_all))
        v2_scores_max = np.max(np.array(v2_scores_all))

        select_result_all_v1 = {'targetname': [targetname] * 5,
                                'model': [i + 1 for i in range(5)],
                                'tmscore': v1_scores_all}

        select_result_avg_v1 = {'targetname': [targetname], 'tmscore': [v1_scores_avg]}

        select_result_max_v1 = {'targetname': [targetname], 'tmscore': [v1_scores_max]}

        select_result_all_v2 = {'targetname': [targetname] * 5,
                                'model': [i + 1 for i in range(5)],
                                'tmscore': v2_scores_all}

        select_result_avg_v2 = {'targetname': [targetname], 'tmscore': [v2_scores_avg]}

        select_result_max_v2 = {'targetname': [targetname], 'tmscore': [v2_scores_max]}

        os.chdir(cwd)

        print(select_result_all_v1)

        return iteration_result_all, iteration_result_avg, iteration_result_max, \
               select_result_all_v1, select_result_avg_v1, select_result_max_v1, \
               select_result_all_v2, select_result_avg_v2, select_result_max_v2

    def search_single(self, chain_id_map, fasta_path, pdb_path, pkl_path, msa_paths, outdir):

        fasta_path = os.path.abspath(fasta_path)

        targetname = pathlib.Path(fasta_path).stem

        print(f"Processing {targetname}")

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        cwd = os.getcwd()

        ref_start_pdb = pdb_path
        ref_start_pkl = pkl_path
        ref_start_msa_paths = msa_paths

        model_iteration_scores = []

        print(f"Start to refine {ref_start_pdb}")

        for num_iteration in range(self.max_iteration):
            os.chdir(cwd)
            current_work_dir = f"{outdir}/iteration{num_iteration + 1}"
            makedir_if_not_exists(current_work_dir)

            start_pdb = f"{current_work_dir}/start.pdb"
            start_pkl = f"{current_work_dir}/start.pkl"
            start_msa_path = f"{current_work_dir}/start_msas"
            if os.path.exists(start_msa_path):
                os.system(f"rm -rf {start_msa_path}")
            makedir_if_not_exists(start_msa_path)

            with open(ref_start_pkl, 'rb') as f:
                ref_avg_lddt = float(pickle.load(f)['ranking_confidence'])

            for chain_id in chain_id_map:
                os.system(f"cp {ref_start_msa_paths[chain_id]['paired_msa']} "
                          f"{start_msa_path}/{chain_id_map[chain_id].description}.start.a3m")

            os.system(f"cp {ref_start_pdb} {start_pdb}")
            os.system(f"cp {ref_start_pkl} {start_pkl}")

            model_iteration_scores += [ref_avg_lddt]

            out_model_dir = f"{current_work_dir}/alphafold"

            if not complete_result(out_model_dir):

                chain_pdbs = split_pdb_unrelax2relax(start_pdb, current_work_dir)

                template_results = []

                out_template_dir = current_work_dir + '/templates'
                makedir_if_not_exists(out_template_dir)

                for chain_id in chain_pdbs:
                    print(chain_id)
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
                          f"--fasta_path {fasta_path} " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir']} " \
                          f"--multimer_a3ms {','.join(msa_files)} " \
                          f"--monomer_a3ms {','.join(msa_files)} " \
                          f"--msa_pair_file {msa_pair_file} " \
                          f"--temp_struct_csv {template_files[0]} " \
                          f"--struct_atom_dir {out_template_dir} " \
                          f"--output_dir {out_model_dir}"
                else:
                    cmd = f"python run_alphafold_multimer_custom_sim.py " \
                          f"--fasta_path {fasta_path} " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir']} " \
                          f"--multimer_a3ms {','.join(msa_files)} " \
                          f"--monomer_a3ms {','.join(msa_files)} " \
                          f"--msa_pair_file {msa_pair_file} " \
                          f"--monomer_temp_csvs {','.join(template_files)} " \
                          f"--struct_atom_dir {out_template_dir} " \
                          f"--output_dir {out_model_dir}"

                try:
                    os.chdir(self.params['alphafold_program_dir_v2'])
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
                ref_start_pdb = f"{out_model_dir}/ranked_0.pdb"
                model_num = list(new_ranking_json["order"])[0].split('_')[1]
                ref_start_pkl = f"{out_model_dir}/result_model_{model_num}_multimer.pkl"
                ref_start_msa_paths = {}
                for chain_id in chain_id_map:
                    ref_start_msa_paths[chain_id] = dict(
                        paired_msa=f"{out_model_dir}/msas/{chain_id_map[chain_id].description}.paired.a3m")
                    # monomer_msa=f"{out_model_dir}/msas/{chain_id}/monomer_final.a3m")

                print('##################################################')
                if num_iteration + 1 >= self.max_iteration:
                    print("Reach maximum iteration")
                    model_iteration_scores += [max_lddt_score]
            else:
                # keep the models in iteration 1 even through the plddt score decreases
                if num_iteration == 0:
                    ref_start_pdb = f"{out_model_dir}/ranked_0.pdb"
                    model_num = list(new_ranking_json["order"])[0].split('_')[1]
                    ref_start_pkl = f"{out_model_dir}/result_model_{model_num}_multimer.pkl"
                    ref_start_msa_paths = {}
                    for chain_id in chain_id_map:
                        ref_start_msa_paths[chain_id] = dict(
                            paired_msa=f"{out_model_dir}/msas/{chain_id_map[chain_id].description}.paired.a3m")
                        # monomer_msa=f"{out_model_dir}/msas/{chain_id}/monomer_final.a3m")
                    model_iteration_scores += [max_lddt_score]
                break

        # model_iteration_scores += [max_lddt_score]
        while len(model_iteration_scores) <= self.max_iteration:
            model_iteration_scores += [0]

        print(model_iteration_scores)
        df = pd.DataFrame(model_iteration_scores)
        df.to_csv(outdir + '/summary.csv')

        final_model_dir = outdir + '/final'
        makedir_if_not_exists(final_model_dir)
        os.system(f"cp {ref_start_pdb} {final_model_dir}/final.pdb")
        os.system(f"cp {ref_start_pkl} {final_model_dir}/final.pkl")

        for chain_id in chain_id_map:
            os.system(f"cp {ref_start_msa_paths[chain_id]['paired_msa']} "
                      f"{final_model_dir}/{chain_id_map[chain_id].description}.paired.a3m")
            # os.system(f"cp {ref_start_msa_paths[chain_id]['monomer_msa']} "
            #           f"{final_model_dir}/{chain_id_map[chain_id].description}.monomer.a3m")

        os.chdir(cwd)

        return final_model_dir

    def search_single_homo(self, chain_id_map, fasta_path, pdb_path, pkl_path, msa_paths, outdir):

        fasta_path = os.path.abspath(fasta_path)

        targetname = pathlib.Path(fasta_path).stem

        print(f"Processing {targetname}")

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        cwd = os.getcwd()

        ref_start_pdb = pdb_path
        ref_start_pkl = pkl_path
        ref_start_msa_paths = msa_paths

        model_iteration_scores = []

        print(f"Start to refine {ref_start_pdb}")

        for num_iteration in range(self.max_iteration):
            os.chdir(cwd)
            current_work_dir = f"{outdir}/iteration{num_iteration + 1}"
            makedir_if_not_exists(current_work_dir)

            start_pdb = f"{current_work_dir}/start.pdb"
            start_pkl = f"{current_work_dir}/start.pkl"
            start_msa_path = f"{current_work_dir}/start_msas"
            if os.path.exists(start_msa_path):
                os.system(f"rm -rf {start_msa_path}")
            makedir_if_not_exists(start_msa_path)

            with open(ref_start_pkl, 'rb') as f:
                ref_avg_lddt = float(pickle.load(f)['ranking_confidence'])

            for chain_id in chain_id_map:
                os.system(f"cp {ref_start_msa_paths[chain_id]['monomer_msa']} "
                          f"{start_msa_path}/{chain_id_map[chain_id].description}.start.a3m")

            os.system(f"cp {ref_start_pdb} {start_pdb}")
            os.system(f"cp {ref_start_pkl} {start_pkl}")

            model_iteration_scores += [ref_avg_lddt]

            out_model_dir = f"{current_work_dir}/alphafold"

            if not complete_result(out_model_dir):

                chain_pdbs = split_pdb_unrelax2relax(start_pdb, current_work_dir)

                template_results = []

                out_template_dir = current_work_dir + '/templates'
                makedir_if_not_exists(out_template_dir)

                for chain_id in chain_pdbs:
                    print(chain_id)
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
                          f"--fasta_path {fasta_path} " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir']} " \
                          f"--multimer_a3ms {','.join(msa_files)} " \
                          f"--monomer_a3ms {','.join(msa_files)} " \
                          f"--temp_struct_csv {template_files[0]} " \
                          f"--struct_atom_dir {out_template_dir} " \
                          f"--output_dir {out_model_dir}"
                else:
                    cmd = f"python run_alphafold_multimer_custom_sim.py " \
                          f"--fasta_path {fasta_path} " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir']} " \
                          f"--multimer_a3ms {','.join(msa_files)} " \
                          f"--monomer_a3ms {','.join(msa_files)} " \
                          f"--monomer_temp_csvs {','.join(template_files)} " \
                          f"--struct_atom_dir {out_template_dir} " \
                          f"--output_dir {out_model_dir}"

                try:
                    os.chdir(self.params['alphafold_program_dir_v2'])
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
                ref_start_pdb = f"{out_model_dir}/ranked_0.pdb"
                model_num = list(new_ranking_json["order"])[0].split('_')[1]
                ref_start_pkl = f"{out_model_dir}/result_model_{model_num}_multimer.pkl"
                ref_start_msa_paths = {}
                for chain_id in chain_id_map:
                    ref_start_msa_paths[chain_id] = dict(
                        # paired_msa=f"{out_model_dir}/msas/{chain_id_map[chain_id].description}.paired.a3m",
                        monomer_msa=f"{out_model_dir}/msas/{chain_id}/monomer_final.a3m")

                print('##################################################')
                if num_iteration + 1 >= self.max_iteration:
                    print("Reach maximum iteration")
                    model_iteration_scores += [max_lddt_score]
            else:
                # keep the models in iteration 1 even through the plddt score decreases
                if num_iteration == 0:
                    ref_start_pdb = f"{out_model_dir}/ranked_0.pdb"
                    model_num = list(new_ranking_json["order"])[0].split('_')[1]
                    ref_start_pkl = f"{out_model_dir}/result_model_{model_num}_multimer.pkl"
                    ref_start_msa_paths = {}
                    for chain_id in chain_id_map:
                        ref_start_msa_paths[chain_id] = dict(
                            # paired_msa=f"{out_model_dir}/msas/{chain_id_map[chain_id].description}.paired.a3m",
                            monomer_msa=f"{out_model_dir}/msas/{chain_id}/monomer_final.a3m")
                    model_iteration_scores += [max_lddt_score]
                break

        # model_iteration_scores += [max_lddt_score]
        while len(model_iteration_scores) <= self.max_iteration:
            model_iteration_scores += [0]

        print(model_iteration_scores)
        df = pd.DataFrame(model_iteration_scores)
        df.to_csv(outdir + '/summary.csv')

        final_model_dir = outdir + '/final'
        makedir_if_not_exists(final_model_dir)
        os.system(f"cp {ref_start_pdb} {final_model_dir}/final.pdb")
        os.system(f"cp {ref_start_pkl} {final_model_dir}/final.pkl")

        for chain_id in chain_id_map:
            # os.system(f"cp {ref_start_msa_paths[chain_id]['paired_msa']} "
            #           f"{final_model_dir}/{chain_id_map[chain_id].description}.paired.a3m")
            os.system(f"cp {ref_start_msa_paths[chain_id]['monomer_msa']} "
                      f"{final_model_dir}/{chain_id_map[chain_id].description}.monomer.a3m")

        os.chdir(cwd)

        return final_model_dir
