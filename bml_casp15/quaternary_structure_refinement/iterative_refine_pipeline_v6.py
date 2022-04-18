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
from bml_casp15.complex_templates_search.sequence_based_pipeline import assess_hhsearch_hit
from bml_casp15.complex_templates_search.parsers import TemplateHit
from bml_casp15.quaternary_structure_refinement.util import *
import itertools


def search_templates_foldseek(foldseek_program, databases, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_runner = Foldseek(binary_path=foldseek_program, databases=databases)
    return foldseek_runner.query_only_local(pdb=inpdb, outdir=outdir)


def reindex_pdb(inpdb):
    resCounter = 0
    atomCounter = 0
    prevrNum = "XX"
    prevchain = "XX"
    contents = []
    chain_sep = {}
    for line in open(inpdb):
        if not line.startswith('ATOM'):
            continue
        rnum = line[22:27]
        chain = line[21:22]
        if prevchain != chain:
            prevchain = chain
            prevrNum = rnum
            resCounter += 1
            chain_sep[chain] = resCounter
        elif prevrNum != rnum:
            prevrNum = rnum
            resCounter += 1
        atomCounter += 1
        rnum_string = "{:>4}".format(resCounter)
        anum_string = "{:>5}".format(atomCounter)
        row = f"{line[:6]}{anum_string}{line[11:21]} {rnum_string}{line[27:]}"
        contents += [row]
    with open(inpdb, 'w') as fw:
        fw.writelines("".join(contents))


def combine_pdbs_to_single_chain(inpdbs, outpdb):
    with open(outpdb, 'w') as fw:
        for inpdb in inpdbs:
            for line in open(inpdb):
                fw.write(line[:21] + ' ' + line[22:])
    reindex_pdb(outpdb)


class complex_template:
    def __init__(self, monomer_templates):
        self.monomer_templates = monomer_templates


class monomer_template:

    def __init__(self,
                 name,
                 template,
                 tstart,
                 tend,
                 qstart,
                 qend,
                 taln,
                 qaln):
        self.name = name
        self.template = template
        self.tstart = tstart
        self.tend = tend
        self.qstart = qstart
        self.qend = qend
        self.qaln = qaln
        self.taln = taln


class Multimer_iterative_refinement_pipeline:

    def __init__(self, params, max_template_count=50):

        self.params = params

        self.max_iteration = 5

        self.max_template_count = max_template_count

    def concatenate_msa_and_templates(self,
                                      chain_id_map,
                                      pair_template_files,
                                      pair_work_dirs,
                                      start_msa_path,
                                      outpath,
                                      iteration):

        # H0957A_H0957B, H0957B_H0957A

        complex_templates_reorder = []

        for pair in pair_template_files:
            templates = pd.read_csv(pair_template_files[pair], sep='\t')

            seq_len_dict = {}
            for chain_id in chain_id_map:
                seq_len_dict[chain_id_map[chain_id].description] = len(chain_id_map[chain_id].sequence)

            print(seq_len_dict)
            # A_B_C, A=60, B=60, C=120
            order = pair.split('_')
            sep_indices = []
            # start=1, end=60
            # start=61, end=60+60=120
            # start=121, end=240
            chain_tempate_dicts = {}
            for i, member in enumerate(order):
                if len(sep_indices) == 0:
                    sep_indices += [dict(start=1,
                                         end=seq_len_dict[member])]
                else:
                    sep_indices += [dict(start=sep_indices[i - 1]['end'] + 1,
                                         end=sep_indices[i - 1]['end'] + seq_len_dict[member])]
                chain_tempate_dicts[member] = []

            # sep_indices = [60, 120]
            # print(sep_indices)
            if len(templates) == 0:
                continue

            for i in range(len(templates)):
                template = templates.loc[i, 'target']
                qaln = templates.loc[i, 'qaln']
                qstart = int(templates.loc[i, 'qstart'])
                qend = int(templates.loc[i, 'qend'])
                taln = templates.loc[i, 'taln']
                tstart = templates.loc[i, 'tstart']
                tend = templates.loc[i, 'tend']
                evalue = float(templates.loc[i, 'evalue'])
                aligned_length = int(templates.loc[i, 'alnlen'])

                # check whether the query sequence contains all the chains
                contain_all_chains = qstart < sep_indices[0]['end'] and qend > sep_indices[len(sep_indices) - 1][
                    'start']
                if not contain_all_chains:
                    continue

                # cut query alignment based on seperate indices
                res_counter = qstart - 1
                chain_qalns = []
                chain_qaln = []
                sep_index = 0
                for char in qaln:
                    if char != '-':
                        res_counter += 1
                    if res_counter > sep_indices[sep_index]['end']:
                        sep_index += 1
                        chain_qalns += [''.join(chain_qaln)]
                        chain_qaln = []
                    chain_qaln += [char]
                chain_qalns += [''.join(chain_qaln)]

                if len(chain_qalns) != len(order):
                    raise Exception(
                        f"The length of seperated query alignments and members in order doesn't match: {len(chain_qalns)} and {len(order)}")

                # print(chain_qalns)

                # qstart = 2, qend = 110
                # qstart = 2, qend = 60
                # qstart = 61, qend = 110 -> qstart = 1, qend = 50
                curr_chain_qstart = qstart
                curr_chain_tstart = tstart
                taln_counter = 0
                for j, member in enumerate(order):
                    chain_qaln = chain_qalns[j]
                    chain_qstart = curr_chain_qstart
                    chain_qend = curr_chain_qstart + len([char for char in chain_qaln if char != '-']) - 1

                    chain_taln = taln[taln_counter:taln_counter + len(chain_qaln)]
                    chain_tstart = curr_chain_tstart
                    chain_tend = curr_chain_tstart + len([char for char in chain_taln if char != '-']) - 1

                    chain_evalue = evalue
                    chain_aligned_length = len(chain_qaln)

                    curr_chain_qstart = chain_qend + 1 - seq_len_dict[member]
                    curr_chain_tstart = chain_tend + 1
                    taln_counter += len(chain_qaln)

                    row_dict = dict(merge_id=i,
                                    template=template.replace('.atom', '') + 'X',
                                    aln_temp=chain_taln,
                                    tstart=chain_tstart,
                                    tend=chain_tend,
                                    aln_query=chain_qaln,
                                    qstart=chain_qstart,
                                    qend=chain_qend,
                                    evalue=chain_evalue,
                                    aligned_length=chain_aligned_length)

                    chain_tempate_dicts[member] += [row_dict]

            prev_df_reorder = None
            for chaind_index, chain_id in enumerate(chain_id_map):
                member_work_dir = pair_work_dirs[pair] + '/' + chain_id_map[chain_id].description
                makedir_if_not_exists(member_work_dir)
                curr_df_reorder = pd.DataFrame(chain_tempate_dicts[chain_id_map[chain_id].description])
                if len(curr_df_reorder) == 0:
                    prev_df_reorder = None
                    break
                curr_df_reorder.to_csv(member_work_dir + '/structure_templates.csv')
                if prev_df_reorder is None:
                    prev_df_reorder = curr_df_reorder
                else:
                    prev_df_reorder = prev_df_reorder.merge(curr_df_reorder, how="inner", on='merge_id',
                                                            suffixes=(str(chaind_index), str(chaind_index + 1)))

            if prev_df_reorder is not None:
                prev_df_reorder.to_csv(pair_work_dirs[pair] + '/structure_templates_reorder.csv')

                complex_templates_reorder += [pair_work_dirs[pair] + '/structure_templates_reorder.csv']

        complex_df = None
        for complex_template_reorder in complex_templates_reorder:
            if complex_df is None:
                complex_df = pd.read_csv(complex_template_reorder)
            else:
                complex_df = complex_df.append(pd.read_csv(complex_template_reorder))
        if complex_df is None:
            return None, None, None

        complex_df.sort_values(by=['evalue1'])
        complex_df.reset_index(inplace=True, drop=True)

        keep_indices = []
        chain_template_msas = {}
        for chain_id in chain_id_map:
            chain_template_msas[chain_id] = {'desc': [chain_id_map[chain_id].description],
                                             'seq': [chain_id_map[chain_id].sequence]}

        seen_complex_seq = []
        seen_complex_seq += ["".join([chain_template_msas[chain_id]['seq'][0] for chain_id in chain_template_msas])]
        for i in range(len(complex_df)):

            if len(keep_indices) >= self.max_template_count:
                break

            monomer_template_seqs = {}
            for j, chain_id in enumerate(chain_id_map):
                query_non_gaps = [res != '-' for res in complex_df.loc[i, f'aln_query{j + 1}']]
                out_sequence = ''.join(convert_taln_seq_to_a3m(query_non_gaps, complex_df.loc[i, f'aln_temp{j + 1}']))
                aln_full = ['-'] * len(chain_id_map[chain_id].sequence)
                qstart = int(complex_df.loc[i, f'qstart{j + 1}'])
                qend = int(complex_df.loc[i, f'qend{j + 1}'])
                aln_full[qstart - 1:qend] = out_sequence
                taln_full_seq = ''.join(aln_full)
                monomer_template_dict = {'desc': complex_df.loc[i, f'template{j + 1}'], 'seq': taln_full_seq}
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
        complex_df.iloc[keep_indices].to_csv(outpath + '/complex_templates.csv')
        return [outpath + '/complex_templates.csv'], out_msas, interact_csv

    def copy_atoms_and_unzip(self, template_csv, outdir):
        os.chdir(outdir)
        templates = pd.read_csv(template_csv, sep='\t')
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'target']
            pdbcode = template_pdb[:4]
            os.system(f"cp {self.params['foldseek_complex_comb_database_dir']}/{template_pdb} {outdir}/{pdbcode}X.atom")

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

                    # B, C
                    chain_pdbs = split_pdb(start_pdb, current_work_dir)

                    # (B, C) or (C, B)
                    combination_pairs = itertools.permutations([chain_id for chain_id in chain_id_map],
                                                               len(chain_id_map))

                    out_template_dir = current_work_dir + '/templates'
                    makedir_if_not_exists(out_template_dir)

                    pair_template_files = {}
                    pair_work_dirs = {}
                    for combination_pair in combination_pairs:
                        pair_fullname = '_'.join([chain_id_map[member].description for member in combination_pair])

                        comb_res_dir = current_work_dir + '/' + pair_fullname
                        makedir_if_not_exists(comb_res_dir)

                        combine_start_pdb = pair_fullname + '.pdb'

                        combine_pdbs_to_single_chain([chain_pdbs[chain_id] for chain_id in combination_pair],
                                                     comb_res_dir + '/' + combine_start_pdb)

                        foldseek_res = search_templates_foldseek(
                            foldseek_program=self.params['foldseek_program'],
                            databases=[self.params['foldseek_complex_comb_database']],
                            inpdb=comb_res_dir + '/' + combine_start_pdb,
                            outdir=comb_res_dir + '/foldseek')

                        query_sequence = ''.join([chain_id_map[chain_id].sequence for chain_id in combination_pair])
                        check_and_rank_monomer_templates_local_or_global(template_result=foldseek_res,
                                                                         outfile=f"{comb_res_dir}/structure_templates.csv",
                                                                         query_sequence=query_sequence,
                                                                         max_template_count=2000)

                        pair_template_files[pair_fullname] = comb_res_dir + '/structure_templates.csv'
                        pair_work_dirs[pair_fullname] = comb_res_dir

                        self.copy_atoms_and_unzip(template_csv=f"{comb_res_dir}/structure_templates.csv",
                                                  outdir=out_template_dir)

                    # if len(pair_template_files) != len(combination_pairs):
                    #     break

                    template_files, msa_files, msa_pair_file = self.concatenate_msa_and_templates(
                        chain_id_map=chain_id_map,
                        pair_template_files=pair_template_files,
                        pair_work_dirs=pair_work_dirs,
                        start_msa_path=start_msa_path,
                        outpath=current_work_dir,
                        iteration=num_iteration + 1)

                    if template_files is None:
                        print(
                            f"Cannot find any templates for {chain_id_map[chain_id].description} in iteration {num_iteration + 1}")
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
