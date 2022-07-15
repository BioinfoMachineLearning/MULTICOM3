import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map


def cal_mmalign(mmalign_program, inpdb, nativepdb):
    cmd = mmalign_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $2}' "
    # print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[0].rstrip('\n'))
    return tmscore


def cal_tmscore(tmscore_program, inpdb, nativepdb, tmpdir):
    cwd = os.getcwd()

    makedir_if_not_exists(tmpdir)

    os.chdir(tmpdir)

    os.system(f"cp {inpdb} inpdb.pdb")
    os.system(f"cp {nativepdb} native.pdb")

    inpdb = "inpdb.pdb"
    nativepdb = "native.pdb"

    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $3}' "
    # print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[2].rstrip('\n'))
    # cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep GDT-score | awk '{print $3}' "
    # tmscore_contents = os.popen(cmd).read().split('\n')
    # gdtscore = float(tmscore_contents[0].rstrip('\n'))

    os.chdir(cwd)

    os.system("rm -rf " + tmpdir)

    return tmscore  # , gdtscore


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_path', type=is_file, required=True)
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--workdir', type=is_dir, required=True)
    parser.add_argument('--refpdb', type=is_file, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    with open(args.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    is_homomer = False
    sequences = []
    for chain_id in chain_id_map:
        if chain_id_map[chain_id].sequence in sequences:
            continue
        sequences += [chain_id_map[chain_id].sequence]

    if len(sequences) > 1:
        print("summarizing results for hetero-multimers\n")
    else:
        print("summarizing results for homo-multimers\n")
        is_homomer = True

    monomer_qa_dir = args.workdir + '/N7_monomer_structure_evaluation'
    if not os.path.exists(monomer_qa_dir):
        monomer_qa_dir = args.workdir + '/N1_monomer_structure_evaluation'

    multimer_qa_dir = args.workdir + '/N9_multimer_structure_evaluation'
    if not os.path.exists(multimer_qa_dir):
        multimer_qa_dir = args.workdir + '/N1_quaternary_structure_evaluation'

    deep_ranking = multimer_qa_dir + '/pairwise_af_avg.ranking'
    qa_ranking = multimer_qa_dir + '/alphafold_ranking.csv'

    processed_seq = []
    first_unit = None
    if os.path.exists(monomer_qa_dir):
        for chain_id in chain_id_map:
            monomer_name = chain_id_map[chain_id].description
            monomer_seq = chain_id_map[chain_id].sequence
            if monomer_seq in processed_seq:
                continue

            print(f"#################Summarizing result for {monomer_name}#################\n")
            if first_unit is None:
                first_unit = monomer_name

            chain_qa_dir = monomer_qa_dir + '/' + monomer_name

            all_models = []
            all_alignments = []
            all_tmscores = []

            egnn_ranking = chain_qa_dir + '/egnn_selected.csv'
            if not os.path.exists(egnn_ranking):
                continue
            alignment_depth = []
            pairwise_ranking_df = pd.read_csv(egnn_ranking)
            ranked_modeles = []
            tmscores = []
            for i in range(len(pairwise_ranking_df)):
                model = pairwise_ranking_df.loc[i, 'selected_models']
                msa = chain_qa_dir + '/msa/' + model.replace('.pdb', '.a3m')
                if os.path.exists(msa):
                    alignment_depth += [len(open(msa).readlines()) / 2]
                else:
                    alignment_depth += [0]
                ranked_modeles += [model]
                tmscore = cal_mmalign(mmalign_program=params['mmalign_program'],
                                      inpdb=chain_qa_dir + '/pdb/' + model,
                                      nativepdb=args.refpdb)
                tmscores += [tmscore]

            print(f"\negnn models: {ranked_modeles}, {alignment_depth}\n")
            all_models += ranked_modeles
            all_alignments += alignment_depth
            all_tmscores += tmscores

            refine_ranking = chain_qa_dir + '/refine_selected.csv'
            refine_ranking_df = pd.read_csv(refine_ranking)
            ranked_modeles = []
            alignment_depth = []
            tmscores = []
            for i in range(5):
                model = refine_ranking_df.loc[i, 'selected_models']
                ranked_modeles += [model]
                msa = chain_qa_dir + '/msa/' + model.replace('.pdb', '.a3m')
                if os.path.exists(msa):
                    alignment_depth += [len(open(msa).readlines()) / 2]
                else:
                    alignment_depth += [0]
                tmscore = cal_mmalign(mmalign_program=params['mmalign_program'],
                                      inpdb=chain_qa_dir + '/pdb/' + model,
                                      nativepdb=args.refpdb)
                tmscores += [tmscore]

            print(f"\nrefine models: {ranked_modeles}\n")
            all_models += ranked_modeles
            all_alignments += alignment_depth
            all_tmscores += tmscores


            pairwise_ranking_df = pd.read_csv(deep_ranking)
            ranked_modeles = []
            alignment_depth = []
            tmscores = []
            for i in range(5):
                model = pairwise_ranking_df.loc[i, 'Name'] + '.pdb'
                msa = multimer_qa_dir + '/msa/' + monomer_name + '/' + model.replace('.pdb', '.monomer.a3m')
                if os.path.exists(msa):
                    alignment_depth += [len(open(msa).readlines()) / 2]
                else:
                    alignment_depth += [0]
                ranked_modeles += [model]

                tmscore = cal_mmalign(mmalign_program=params['mmalign_program'],
                                      inpdb=f"{multimer_qa_dir}/deep{i+1}/{chain_id}_top1.pdb",
                                      nativepdb=args.refpdb)

                tmscores += [tmscore]

            print(f"\ndeep models: {ranked_modeles}\n")
            all_models += ranked_modeles
            all_alignments += alignment_depth
            all_tmscores += tmscores

            qa_ranking_df = pd.read_csv(qa_ranking)
            ranked_modeles = []
            alignment_depth = []
            tmscores = []
            for i in range(5):
                model = qa_ranking_df.loc[i, 'model']
                ranked_modeles += [model]
                msa = multimer_qa_dir + '/msa/' + monomer_name + '/' + model.replace('.pdb', '.monomer.a3m')
                if os.path.exists(msa):
                    alignment_depth += [len(open(msa).readlines()) / 2]
                else:
                    alignment_depth += [0]
                tmscore = cal_mmalign(mmalign_program=params['mmalign_program'],
                                      inpdb=f"{multimer_qa_dir}/qa{i+1}/{chain_id}_top1.pdb",
                                      nativepdb=args.refpdb)

                tmscores += [tmscore]
            print(f"\nqa models: {ranked_modeles}\n")
            all_models += ranked_modeles
            all_alignments += alignment_depth
            all_tmscores += tmscores

            df = pd.DataFrame({'model': all_models, 'alignment_depth': all_alignments, 'ref_tmscore': all_tmscores})
            print(df)
            df.to_csv(f"{args.workdir}/summary_{monomer_name}_monomer.csv")

            pdb70_hit_file = f"{args.workdir}/N3_monomer_structure_generation/{monomer_name}/default/msas/pdb_hits.hhr"
            if os.path.exists(pdb70_hit_file):
                start = False
                contents = open(pdb70_hit_file).readlines()
                for i, line in enumerate(contents):
                    line = line.rstrip('\n')
                    if line == "No 1":
                        template_name = contents[i + 1].split()[0].lstrip('>')
                        evalue = contents[i + 2].split()[1].split('=')[1]
                        print(f"\nPDB70 Template: {template_name}, e-value: {evalue}")
                        break

            pdb_hit_file = f"{args.workdir}/N2_monomer_template_search/{monomer_name}/output.hhr"
            if os.path.exists(pdb_hit_file):
                start = False
                contents = open(pdb_hit_file).readlines()
                for i, line in enumerate(contents):
                    line = line.rstrip('\n')
                    if line == "No 1":
                        template_name = contents[i + 1].split()[0].lstrip('>')
                        evalue = contents[i + 2].split()[1].split('=')[1]
                        print(f"\nIn house PDB Template: {template_name}, e-value: {evalue}")
                        break

            processed_seq += [monomer_seq]

    print(f"#################Summarizing result for the multimer##################\n")

    all_models = []
    all_alignments = []
    all_tmscores = []


    alignment_depth = []
    pairwise_ranking_df = pd.read_csv(deep_ranking)
    ranked_modeles = []
    tmscores = []
    for i in range(5):
        model = pairwise_ranking_df.loc[i, 'Name'] + '.pdb'
        if first_unit is None:
            msa = ""
        else:
            if is_homomer:
                msa = multimer_qa_dir + '/msa/' + first_unit + '/' + model.replace('.pdb', '.monomer.a3m')
            else:
                msa = multimer_qa_dir + '/msa/' + first_unit + '/' + model.replace('.pdb', '.paired.a3m')
        if os.path.exists(msa):
            alignment_depth += [len(open(msa).readlines()) / 2]
        else:
            alignment_depth += [0]
        ranked_modeles += [model]

        tmscore = cal_mmalign(mmalign_program=params['mmalign_program'],
                              inpdb=multimer_qa_dir + '/pdb/' + model,
                              nativepdb=args.refpdb)

        tmscores += [tmscore]

    print(f"\ndeep models: {ranked_modeles}, {alignment_depth}\n")
    all_models += ranked_modeles
    all_alignments += alignment_depth
    all_tmscores += tmscores


    qa_ranking_df = pd.read_csv(qa_ranking)
    ranked_modeles = []
    alignment_depth = []
    tmscores = []
    for i in range(5):
        model = qa_ranking_df.loc[i, 'model']
        ranked_modeles += [model]
        if first_unit is None:
            msa = ""
        else:
            if is_homomer:
                msa = multimer_qa_dir + '/msa/' + first_unit + '/' + model.replace('.pdb', '.monomer.a3m')
            else:
                msa = multimer_qa_dir + '/msa/' + first_unit + '/' + model.replace('.pdb', '.paired.a3m')
        if os.path.exists(msa):
            alignment_depth += [len(open(msa).readlines()) / 2]
        else:
            alignment_depth += [0]
        tmscore = cal_mmalign(mmalign_program=params['mmalign_program'],
                              inpdb=multimer_qa_dir + '/pdb/' + model,
                              nativepdb=args.refpdb)

        tmscores += [tmscore]

    print(f"\nqa models: {ranked_modeles}\n")
    all_models += ranked_modeles
    all_alignments += alignment_depth
    all_tmscores += tmscores

    df = pd.DataFrame({'model': all_models, 'alignment_depth': all_alignments, 'tmscores': all_tmscores})
    print(df)
    df.to_csv(args.workdir + '/summary_multimer.csv')
