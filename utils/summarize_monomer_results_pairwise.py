import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.common.protein import read_qa_txt_as_df


def cal_tmscore(tmscore_program, inpdb, nativepdb, tmpdir):
    if nativepdb is None or not os.path.exists(nativepdb):
        return 0, 0

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
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep GDT-score | awk '{print $3}' "
    tmscore_contents = os.popen(cmd).read().split('\n')
    gdtscore = float(tmscore_contents[0].rstrip('\n'))

    os.chdir(cwd)

    os.system("rm -rf " + tmpdir)

    return tmscore, gdtscore


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--workdir', type=is_file, required=True)
    parser.add_argument('--tmpdir', type=str, required=True)
    parser.add_argument('--refpdb', type=is_file, required=False)

    args = parser.parse_args()

    makedir_if_not_exists(args.tmpdir)

    params = read_option_file(args.option_file)

    qa_dir = args.workdir + '/N4_monomer_structure_evaluation'

    if os.path.exists(args.workdir + '/N6_monomer_structure_evaluation'):
        print("Using img alignment in the pipeline")
        qa_dir = args.workdir + '/N6_monomer_structure_evaluation'

    if not os.path.exists(qa_dir):
        print("Rerunning the evaluation pipeline")
        qa_dir = args.workdir + '/N1_monomer_structure_evaluation/'
        qa_dir = args.workdir + '/N1_monomer_structure_evaluation//' + os.listdir(qa_dir)[0]

    all_models = []
    all_alignments = []
    all_tmscores = []

    egnn_ranking = qa_dir + '/egnn_selected.csv'
    alignment_depth = []
    pairwise_ranking_df = pd.read_csv(egnn_ranking)
    ranked_modeles = []
    for i in range(len(pairwise_ranking_df)):
        model = pairwise_ranking_df.loc[i, 'selected_models']
        msa = qa_dir + '/msa/' + model.replace('.pdb', '.a3m')
        if os.path.exists(msa):
            alignment_depth += [len(open(msa).readlines()) / 2]
        else:
            alignment_depth += [0]

        ranked_modeles += [model]

        tmscore, _ = cal_tmscore(tmscore_program=params['tmscore_program'],
                                 inpdb=qa_dir + '/pdb/' + model,
                                 nativepdb=args.refpdb, tmpdir=args.tmpdir)
        all_tmscores += [tmscore]

    print(f"\negnn models: {ranked_modeles}, {alignment_depth}\n")
    all_models += ranked_modeles
    all_alignments += alignment_depth

    refine_ranking = qa_dir + '/pairwise_ranking.tm'
    if os.path.exists(refine_ranking):
        refine_ranking_df = read_qa_txt_as_df(refine_ranking)
        ranked_modeles = []
        alignment_depth = []
        for i in range(5):
            model = refine_ranking_df.loc[i, 'model']
            ranked_modeles += [model]
            msa = qa_dir + '/msa/' + model.replace('.pdb', '.a3m')
            if os.path.exists(msa):
                alignment_depth += [len(open(msa).readlines()) / 2]
            else:
                alignment_depth += [0]
            tmscore, _ = cal_tmscore(tmscore_program=params['tmscore_program'],
                                     inpdb=qa_dir + '/pdb/' + model,
                                     nativepdb=args.refpdb, tmpdir=args.tmpdir)
            all_tmscores += [tmscore]

        print(f"\nrefine models: {ranked_modeles}\n")
        all_models += ranked_modeles
        all_alignments += alignment_depth

    deep_ranking = qa_dir + '/deep_selected.csv'
    alignment_depth = []
    pairwise_ranking_df = pd.read_csv(deep_ranking)
    ranked_modeles = []
    for i in range(len(pairwise_ranking_df)):
        model = pairwise_ranking_df.loc[i, 'selected_models']
        msa = qa_dir + '/msa/' + model.replace('.pdb', '.a3m')
        if os.path.exists(msa):
            alignment_depth += [len(open(msa).readlines()) / 2]
        else:
            alignment_depth += [0]
        ranked_modeles += [model]
        tmscore, _ = cal_tmscore(tmscore_program=params['tmscore_program'],
                                 inpdb=qa_dir + '/pdb/' + model,
                                 nativepdb=args.refpdb, tmpdir=args.tmpdir)
        all_tmscores += [tmscore]

    print(f"\ndeep models: {ranked_modeles}, {alignment_depth}\n")
    all_models += ranked_modeles
    all_alignments += alignment_depth

    qa_ranking = qa_dir + '/pairwise_ranking.gdt'
    if os.path.exists(refine_ranking):
        refine_ranking_df = read_qa_txt_as_df(qa_ranking)
        ranked_modeles = []
        alignment_depth = []
        for i in range(5):
            model = refine_ranking_df.loc[i, 'model']
            ranked_modeles += [model]
            msa = qa_dir + '/msa/' + model.replace('.pdb', '.a3m')
            if os.path.exists(msa):
                alignment_depth += [len(open(msa).readlines()) / 2]
            else:
                alignment_depth += [0]
            tmscore, _ = cal_tmscore(tmscore_program=params['tmscore_program'],
                                     inpdb=qa_dir + '/pdb/' + model,
                                     nativepdb=args.refpdb, tmpdir=args.tmpdir)
            all_tmscores += [tmscore]

        print(f"\nqa models: {ranked_modeles}\n")
        all_models += ranked_modeles
        all_alignments += alignment_depth

    df = pd.DataFrame({'model': all_models, 'alignment_depth': all_alignments, 'tmscore': all_tmscores})
    print(df)
    df.to_csv(args.workdir + '/summary.csv')