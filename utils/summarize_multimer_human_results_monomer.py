import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir


def cal_mmalign(mmalign_program, inpdb, nativepdb):
    cmd = mmalign_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $2}' "
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[1].rstrip('\n'))
    return tmscore


def convert_ranking_to_df(method, infile, index, pdb_dir, refpdb, tmscore_progarm,
                          ranked_field="", is_csv=True, ascending=False, index_col='index'):
    df = None
    # af=plddt_avg, enqa=score, native_scores='gdtscore
    if is_csv:
        if index_col is None:
            df = pd.read_csv(infile, index_col=index_col)
        else:
            df = pd.read_csv(infile)
        if ranked_field != 'score':
            df['score'] = df[ranked_field]
        tmscores = []
        if method == "multieva":
            df['model'] = df['Name'] + '.pdb'
        for i in range(len(df)):
            model = df.loc[i, 'model']
            tmscore = cal_mmalign(mmalign_program=tmscore_progarm,
                                     inpdb=pdb_dir + '/' + model,
                                     nativepdb=refpdb)
            tmscores += [tmscore]
        df['tmscore'] = tmscores
        df = df.sort_values(by=['score'], ascending=ascending)
        df.reset_index(inplace=True, drop=True)
        df = df.add_suffix(str(index))
        df = df[[f'model{index}', f'score{index}', f'tmscore{index}']]
    else:
        models = []
        scores = []
        for line in open(infile):
            line = line.rstrip('\n')
            contents = line.split()
            if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or \
                    contents[0] == "END":
                continue
            contents = line.split()
            model = contents[0]
            score = contents[1]
            if model.find('BML_CASP15') > 0:
                model = pathlib.Path(model).name
            if model.find('.pdb') < 0:
                model = model + '.pdb'
            models += [model]
            scores += [float(score)]
        df = pd.DataFrame({f'model{index}': models, f'score{index}': scores})
        tmscores = []
        for i in range(len(df)):
            model = df.loc[i, f'model{index}']
            if method == "SBROD":
                model = model[model.rindex('/'):]

            tmscore = cal_mmalign(mmalign_program=tmscore_progarm,
                                  inpdb=pdb_dir + '/' + model,
                                  nativepdb=refpdb)
            tmscores += [tmscore]
        df[f'tmscore{index}'] = tmscores
        df = df.sort_values(by=[f'score{index}'], ascending=False)
        df.reset_index(inplace=True)
        df = df[[f'model{index}', f'score{index}', f'tmscore{index}']]

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--workdir', type=is_file, required=True)
    parser.add_argument('--refpdb', type=str, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    qa_dir = args.workdir + '/N1_quaternary_structure_evaluation/'
    pdb_dir = qa_dir + '/pdb'

    ranking_csvs = {'af_multieva_avg': f"{qa_dir}/pairwise_af_avg.ranking",
                    'bfactor_multieva_avg': f"{qa_dir}/pairwise_bfactor_avg.ranking",
                    'af': f"{qa_dir}/alphafold_ranking.csv",
                    'bfactor': f"{qa_dir}/bfactor_ranking.csv",
                    'multieva': f"{qa_dir}/multieva.csv"}

    print(' '.join([method for method in ranking_csvs]))
    # print(ranking_csvs)
    miss_res = []
    find_all_res = True
    for method in ranking_csvs:
        if not os.path.exists(ranking_csvs[method]):
            print(f"cannot find {ranking_csvs[method]}")
            find_all_res = False
            break

    if not find_all_res:
        raise Exception("Cannot find all the results")

    all_dfs = []

    for index, method in enumerate(ranking_csvs):
        print(ranking_csvs[method])
        global_scores = None
        if method == 'af':
            global_scores = convert_ranking_to_df(method=method,
                                                  infile=ranking_csvs[method],
                                                  index=index,
                                                  pdb_dir=pdb_dir,
                                                  refpdb=args.refpdb,
                                                  tmscore_progarm=params['mmalign_program'],
                                                  ranked_field='confidence',
                                                  is_csv=True)
        elif method == 'bfactor':
            global_scores = convert_ranking_to_df(method=method,
                                                  infile=ranking_csvs[method],
                                                  index=index,
                                                  pdb_dir=pdb_dir,
                                                  refpdb=args.refpdb,
                                                  tmscore_progarm=params['mmalign_program'],
                                                  ranked_field='bfactor',
                                                  is_csv=True)
        # elif method == 'dproq_dockq' or method == 'dproq_evalue':
        #     global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
        #                                           ranked_field='score',
        #                                           is_csv=True)
        # elif method == "foldseek":
        #     global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
        #                                           ranked_field='score2',
        #                                           is_csv=True)
        elif method == "multieva":
            global_scores = convert_ranking_to_df(method=method,
                                                  infile=ranking_csvs[method],
                                                  index=index,
                                                  pdb_dir=pdb_dir,
                                                  refpdb=args.refpdb,
                                                  tmscore_progarm=params['mmalign_program'],
                                                  ranked_field='average_MMS',
                                                  is_csv=True, index_col=None)
        # elif method == 'enqa':
        #     global_scores = convert_ranking_to_df(infile=ranking_csvs[method],
        #                                           ranked_field='score',
        #                                           is_csv=True)
        elif method == "af_multieva_avg" or method == "bfactor_multieva_avg":
            global_scores = convert_ranking_to_df(method=method,
                                                  infile=ranking_csvs[method],
                                                  ranked_field='avg_score',
                                                  index=index,
                                                  pdb_dir=pdb_dir,
                                                  refpdb=args.refpdb,
                                                  tmscore_progarm=params['mmalign_program'],
                                                  is_csv=True)

        keep_indices = []
        for i in range(len(global_scores)):
            if global_scores.iloc[i, f'model{index}'].find('_A.pdb') > 0 or \
                global_scores.iloc[i, f'model{index}'].find('_B.pdb') > 0:
                continue
            keep_indices += [i]

        print(global_scores.loc[keep_indices])
        all_dfs += [global_scores.loc[keep_indices]]

    summary_df = pd.concat(all_dfs, axis=1)
    summary_df.to_csv(args.workdir + '/summary_multimer.csv')
    summary_df.head(20).to_csv(args.workdir + '/summary_multimer_20.csv')
    print(summary_df)
