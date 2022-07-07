import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir


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
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep GDT-score | awk '{print $3}' "
    tmscore_contents = os.popen(cmd).read().split('\n')
    gdtscore = float(tmscore_contents[0].rstrip('\n'))

    os.chdir(cwd)

    os.system("rm -rf " + tmpdir)

    return tmscore, gdtscore

def convert_ranking_to_df(method, infile, index,
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
            tmscores += [0]
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

            tmscores += [0]
        df[f'tmscore{index}'] = tmscores
        df = df.sort_values(by=[f'score{index}'], ascending=False)
        df.reset_index(inplace=True)
        df = df[[f'model{index}', f'score{index}', f'tmscore{index}']]

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--workdirs', type=str, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    all_dfs = {}

    for workdir in args.workdirs.split(','):

        ranking_csvs = {'af_multieva_avg': f"{workdir}/pairwise_af_avg.ranking",
                        'bfactor_multieva_avg': f"{workdir}/pairwise_bfactor_avg.ranking",
                        'af': f"{workdir}/alphafold_ranking.csv",
                        'bfactor': f"{workdir}/bfactor_ranking.csv",
                        'multieva': f"{workdir}/multieva.csv"}

        pdb_dir = workdir + '/pdb'
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

        for index, method in enumerate(ranking_csvs):
            print(ranking_csvs[method])
            global_scores = None
            if method == 'af':
                global_scores = convert_ranking_to_df(method=method,
                                                      infile=ranking_csvs[method],
                                                      index=index,
                                                      ranked_field='plddt_avg',
                                                      is_csv=True)
            elif method == 'bfactor':
                global_scores = convert_ranking_to_df(method=method,
                                                      infile=ranking_csvs[method],
                                                      index=index,
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
                                                      is_csv=True)

            if method in all_dfs:
                all_dfs[method] += [global_scores]
            else:
                all_dfs[method] = [global_scores]

    all_df_avg = []
    for index, method in enumerate(all_dfs):
        prev_df = all_dfs[method][0]
        prev_df = prev_df.add_suffix('0')
        prev_df[f'model{index}'] = prev_df[f'model{index}0']
        prev_df = prev_df.drop([f'model{index}0'], axis=1)
        print(prev_df)
        for j in range(1, len(all_dfs[method])):
            curr_df = all_dfs[method][j].add_suffix(str(j))
            curr_df[f'model{index}'] = curr_df[f'model{index}{j}']
            curr_df = curr_df.drop([f'model{index}{j}'], axis=1)
            prev_df = prev_df.merge(curr_df, on=f'model{index}')

        avg_scores = []
        for j in range(len(prev_df)):
            sum_score = 0
            for k in range(len(all_dfs[method])):
                sum_score += prev_df.loc[j, f'score{index}{k}']
            avg_scores += [sum_score/len(all_dfs[method])]

        models = list(prev_df[f'model{index}'])

        df = pd.DataFrame({f'model{index}': models, f'score{index}': avg_scores, f'tmscore{index}': [0]*len(avg_scores)})
        df = df.sort_values(by=[f'score{index}'], ascending=False)
        df.reset_index(inplace=True)
        print(method)
        print(df)
        all_df_avg += [df[[f'model{index}', f'score{index}', f'tmscore{index}']]]

    summary_df = pd.concat(all_df_avg, axis=1)
    summary_df.to_csv('summary_multimer_avg.csv')
    summary_df.head(20).to_csv('summary_multimer_avg_20.csv')
    print(summary_df)