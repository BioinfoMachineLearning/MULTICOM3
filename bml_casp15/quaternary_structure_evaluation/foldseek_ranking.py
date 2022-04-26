import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from bml_casp15.common.util import makedir_if_not_exists
from bml_casp15.common.protein import extract_pdb
import glob
import itertools
from bml_casp15.tool.foldseek import *
import copy
from bml_casp15.quaternary_structure_refinement.util import cal_tmscore


def search_templates_foldseek(foldseek_program, databases, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_runner = Foldseek(binary_path=foldseek_program, databases=databases)
    return foldseek_runner.query(pdb=inpdb, outdir=outdir, maxseq=300)


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


def split_pdb(complex_pdb, outdir):
    pre_chain = None
    chain_pdbs = {}
    resCounter = 0
    prevrNum = "XX"
    prevchain = "XX"
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue

        rnum = line[22:27]
        chain_name = line[21:22]
        if prevchain != chain_name:
            prevchain = chain_name
            prevrNum = rnum
            resCounter += 1
        elif prevrNum != rnum:
            prevrNum = rnum
            resCounter += 1

        if pre_chain is None:
            pre_chain = chain_name
            fw = open(outdir + '/' + pre_chain + '.pdb', 'w')
            fw.write(line)
            chain_pdbs[pre_chain] = outdir + '/' + pre_chain + '.pdb'
        elif chain_name == pre_chain:
            fw.write(line)
        else:
            fw.close()
            pre_chain = chain_name
            fw = open(outdir + '/' + pre_chain + '.pdb', 'w')
            fw.write(line)
            chain_pdbs[pre_chain] = outdir + '/' + pre_chain + '.pdb'
    fw.close()
    chain_pdbs[pre_chain] = outdir + '/' + pre_chain + '.pdb'
    return chain_pdbs, resCounter


def combine_pdbs_to_single_chain(inpdbs, outpdb):
    with open(outpdb, 'w') as fw:
        for inpdb in inpdbs:
            for line in open(inpdb):
                fw.write(line[:21] + ' ' + line[22:])
    reindex_pdb(outpdb)


def merge_chains(inpdb, outpdb):
    contents = []
    for line in open(inpdb):
        if line.startswith('ATOM'):
            contents += [line[:21] + ' ' + line[22:]]
    with open(outpdb, 'w') as fw:
        fw.writelines(contents)
    reindex_pdb(outpdb)


def cal_mmalign(inparams):
    mmalign_program, src_pdb, trg_pdb, alnlen, total_length = inparams
    tmscore = cal_tmscore(mmalign_program, src_pdb, trg_pdb)
    tmscore = tmscore * alnlen / total_length
    return tmscore


class FoldSeek_qa:

    def __init__(self, params):
        self.params = params

    def run(self, chain_id_map, input_dir, outputdir):

        cwd = os.getcwd()

        input_dir = os.path.abspath(input_dir)

        outputdir = os.path.abspath(outputdir)

        if os.path.exists(outputdir + '/foldseek_qa.csv'):
            return pd.read_csv(outputdir + '/foldseek_qa.csv')

        makedir_if_not_exists(outputdir)

        pdb_foldseek_dir = outputdir + '/pdb_foldseek'
        makedir_if_not_exists(pdb_foldseek_dir)

        for pdb in os.listdir(input_dir):
            merge_chains(input_dir + '/' + pdb, outputdir + '/pdb_foldseek/' + pdb)

        os.chdir(outputdir)

        os.system(f"{self.params['foldseek_program']} createdb {outputdir}/pdb_foldseek pdbs")

        models = []
        scores_0 = []
        scores_2 = []
        for pdb in os.listdir(input_dir):

            pdbdir = outputdir + '/' + pdb
            makedir_if_not_exists(pdbdir)
            chain_pdbs, model_length = split_pdb(input_dir + '/' + pdb, pdbdir)

            combination_pairs = itertools.permutations([chain_id for chain_id in chain_pdbs], len(chain_pdbs))

            pair_template_files = {}
            pair_work_dirs = {}
            for combination_pair in combination_pairs:
                pair_fullname = '_'.join([chain_id for chain_id in combination_pair])

                comb_res_dir = pdbdir + '/' + pair_fullname

                makedir_if_not_exists(comb_res_dir)

                combine_start_pdb = pair_fullname + '.pdb'

                combine_pdbs_to_single_chain([chain_pdbs[chain_id] for chain_id in combination_pair],
                                             comb_res_dir + '/' + combine_start_pdb)

                foldseek_res = search_templates_foldseek(
                    foldseek_program=self.params['foldseek_program'],
                    databases=['pdbs'],
                    inpdb=comb_res_dir + '/' + combine_start_pdb,
                    outdir=comb_res_dir + '/foldseek')

                foldseek_res['all_alignment'].to_csv(comb_res_dir + '/structure_templates.csv')
                pair_template_files[pair_fullname] = comb_res_dir + '/structure_templates.csv'
                pair_work_dirs[pair_fullname] = comb_res_dir

            keep_indices = []
            result_df = None
            for pair_fullname in pair_template_files:
                pair_df = pd.read_csv(pair_template_files[pair_fullname])
                if result_df is None:
                    result_df = pair_df
                else:
                    result_df = result_df.append(pair_df)

            result_df = result_df.sort_values(by=['target'])
            result_df.reset_index(inplace=True, drop=True)
            result_df.to_csv(pdbdir + '/combined_templates.csv')

            prev_template_name = result_df.loc[0, 'target']
            min_evalue = result_df.loc[0, 'evalue']
            min_index = 0
            for i in range(1, len(result_df)):
                if result_df.loc[i, 'target'] != prev_template_name:
                    keep_indices += [min_index]
                    prev_template_name = result_df.loc[i, 'target']
                    min_evalue = result_df.loc[i, 'evalue']
                    min_index = i
                else:
                    if result_df.loc[i, 'evalue'] < min_evalue:
                        min_evalue = result_df.loc[i, 'evalue']
                        min_index = i

            if min_index != 0:
                keep_indices += [min_index]

            templates_filtered = copy.deepcopy(result_df.iloc[keep_indices])
            templates_filtered.to_csv(pdbdir + '/final_templates_unsorted.csv')
            templates_filtered = templates_filtered.sort_values(by=['evalue'])
            templates_filtered.drop(templates_filtered.filter(regex="Unnamed"), axis=1, inplace=True)
            templates_filtered.reset_index(inplace=True, drop=True)
            templates_filtered.to_csv(pdbdir + '/final_templates.csv')

            models += [pdb]

            # option 0:
            score_0_sum = 0
            score_2_sum = 0
            total_seq_length = np.sum(np.array([len(chain_id_map[chain_id].sequence) for chain_id in chain_id_map]))
            print(templates_filtered)

            processs_list = []
            for i in range(len(templates_filtered)):
                if templates_filtered.loc[i, 'target'] == pdb:
                    print(f"filter out {pdb}")
                    continue
                alnlen = templates_filtered.loc[i, 'alnlen']
                score = alnlen / model_length
                score_0_sum += score

                print(f"comparing {pdb} and {templates_filtered.loc[i, 'target']}")
                qstart = templates_filtered.loc[i, 'qstart']
                qend = templates_filtered.loc[i, 'qend']
                extract_pdb(outputdir + '/pdb_foldseek/' + pdb, pdbdir + '/query.pdb', qstart, qend)

                target_pdb = templates_filtered.loc[i, 'target']
                tstart = templates_filtered.loc[i, 'tstart']
                tend = templates_filtered.loc[i, 'tend']
                extract_pdb(outputdir + '/pdb_foldseek/' + target_pdb, pdbdir + '/target.pdb', tstart, tend)

                processs_list.append([self.params['mmalign_program'], pdbdir + '/query.pdb', pdbdir + '/target.pdb',
                                      templates_filtered.loc[i, 'alnlen'], total_seq_length])

            pool = Pool(processes=40)
            results = pool.map(cal_mmalign, processs_list)
            pool.close()
            pool.join()

            for result in results:
                score_2_sum += result

            score_0 = score_0_sum / len(templates_filtered)
            scores_0 += [score_0]

            # option 2:
            score_2 = score_2_sum / len(templates_filtered)
            scores_2 += [score_2]

        ranking_df = pd.DataFrame({'model': models, 'score0': scores_0, 'score2': scores_2})
        ranking_df = ranking_df.sort_values(by=['score0'], ascending=False)
        ranking_df.reset_index(inplace=True, drop=True)
        ranking_df.to_csv('foldseek_qa.csv')
        os.chdir(cwd)

        return ranking_df
