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


def convert_ranking_to_df(method, infile, index, pdb_dir, refpdb, tmscore_progarm, tmpdir,
                          ranked_field="", is_csv=True, ascending=False):
    df = None
    # af=plddt_avg, enqa=score, native_scores='gdtscore
    if is_csv:
        df = pd.read_csv(infile)
        if ranked_field != 'score':
            df['score'] = df[ranked_field]
        tmscores = []
        for i in range(len(df)):
            model = df.loc[i, 'model']
            tmscore, _ = cal_tmscore(tmscore_program=tmscore_progarm,
                                     inpdb=pdb_dir + '/' + model,
                                     nativepdb=refpdb,
                                     tmpdir=tmpdir)
            tmscores += [tmscore]
        df['tmscore'] = tmscores
        df = df.sort_values(by=['score'], ascending=ascending)
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

            tmscore, _ = cal_tmscore(tmscore_program=tmscore_progarm,
                                     inpdb=pdb_dir + '/' + model,
                                     nativepdb=refpdb,
                                     tmpdir=tmpdir)
            tmscores += [tmscore]
        df[f'tmscore{index}'] = tmscores
        df = df.sort_values(by=[f'score{index}'], ascending=False)
        df.reset_index(inplace=True)
        df = df[[f'model{index}', f'score{index}', f'tmscore{index}']]

    return df


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--workdirs', type=str, required=True)
    parser.add_argument('--tmpdir', type=str, required=True)
    parser.add_argument('--targetname', type=str, required=True)
    parser.add_argument('--refpdb', type=str, required=True)

    args = parser.parse_args()

    makedir_if_not_exists(args.tmpdir)

    params = read_option_file(args.option_file)

    all_dfs = {}

    for workdir in args.workdirs.split(','):

        ranking_csvs = {'af_apollo_avg': f"{workdir}/pairwise_af_avg.ranking",
                        'af': f"{workdir}/alphafold_ranking.csv",
                        'apollo': f"{workdir}/pairwise_ranking.tm",
                        'deeprank': f"{workdir}/casp13/Full_TS/eva/HumanQA_gdt_prediction.txt",
                        'deeprank_dncon4': f"{workdir}/casp14/Full_TS/eva/HumanQA_gdt_prediction.txt",
                        'deeprank3_cluster': f"{workdir}/casp14/Full_TS/eva_DeepRank3/DeepRank3_Cluster.txt",
                        'deeprank3_singleqa': f"{workdir}/casp14/Full_TS/eva_DeepRank3/DeepRank3_SingleQA.txt",
                        'deeprank3_singleqa_lite': f"{workdir}/casp14/Full_TS/eva_DeepRank3/DeepRank3_SingleQA_lite.txt",
                        'enqa': f"{workdir}/enqa_ranking.csv",
                        'modfoldclust2': f"{workdir}/casp13/Full_TS/eva/ALL_scores/modfoldclust2.{target}"}

        for method in ['pcons', 'DeepQA', 'dfire', 'dncon2_long-range', 'dncon2_medium-range', 'dncon2_short-range',
                       'dope', 'OPUS', 'proq2', 'proq3', 'RWplus', 'SBROD', 'solvent',
                       'ss_sim']:
            ranking_csvs[method] = f"{workdir}/casp13/Full_TS/eva/ALL_scores/feature_{method}.{target}"

        for method in ['dist_pearson', 'dist_ssim', 'dncon4_long-range']:
            ranking_csvs[method] = f"{workdir}/casp14/Full_TS/eva_DeepRank3/ALL_scores/feature_{method}.{target}"

        # print(ranking_csvs)
        find_res = {}
        for method in ranking_csvs:
            if not os.path.exists(ranking_csvs[method]):
                print(f"cannot find {ranking_csvs[method]}")
            else:
                find_res[method] = ranking_csvs[method]

        print(' '.join([method for method in find_res]))

        for index, method in enumerate(find_res):
            print(find_res[method])
            if method == 'af':
                global_scores = convert_ranking_to_df(method=method,
                                                      infile=find_res[method],
                                                      index=index,
                                                      pdb_dir=pdb_dir,
                                                      refpdb=args.refpdb,
                                                      tmscore_progarm=params['tmscore_program'],
                                                      tmpdir=args.tmpdir,
                                                      ranked_field='plddt_avg',
                                                      is_csv=True)
            elif method == 'af_apollo_avg':
                global_scores = convert_ranking_to_df(method=method,
                                                      infile=find_res[method],
                                                      index=index,
                                                      pdb_dir=pdb_dir,
                                                      refpdb=args.refpdb,
                                                      tmscore_progarm=params['tmscore_program'],
                                                      tmpdir=args.tmpdir,
                                                      ranked_field='avg_score',
                                                      is_csv=True)
            elif method == 'af_apollo_avg_rank':
                global_scores = convert_ranking_to_df(method=method,
                                                      infile=find_res[method],
                                                      index=index,
                                                      pdb_dir=pdb_dir,
                                                      refpdb=args.refpdb,
                                                      tmscore_progarm=params['tmscore_program'],
                                                      tmpdir=args.tmpdir,
                                                      ranked_field='avg_rank',
                                                      is_csv=True, ascending=True)
            elif method == 'enqa':
                global_scores = convert_ranking_to_df(method=method,
                                                      infile=find_res[method],
                                                      index=index,
                                                      pdb_dir=pdb_dir,
                                                      refpdb=args.refpdb,
                                                      tmscore_progarm=params['tmscore_program'],
                                                      tmpdir=args.tmpdir,
                                                      ranked_field='score',
                                                      is_csv=True)
            else:
                global_scores = convert_ranking_to_df(method=method,
                                                      infile=find_res[method],
                                                      index=index,
                                                      pdb_dir=pdb_dir,
                                                      refpdb=args.refpdb,
                                                      tmscore_progarm=params['tmscore_program'],
                                                      tmpdir=args.tmpdir,
                                                      is_csv=False)

            if method in all_dfs:
                all_dfs[method] += [global_scores]
            else:
                all_dfs[method] = [global_scores]

    all_df_avg = []
    for i, method in enumerate(all_dfs):
        prev_df = all_dfs[method][0]
        for j in range(1, len(all_dfs[method])):
            prev_df = prev_df.merge(all_dfs[method][j], on=f'model{index}', suffixes=(str(j), str(j+1)))

        avg_scores = []
        for j in range(len(prev_df)):
            sum_score = 0
            for k in range(len(all_dfs[method])):
                sum_score += prev_df.loc[j, f'score{k}']
            avg_scores += [sum_score/len(all_dfs[method])]

        models = list(prev_df[f'model{index}1'])

        df = pd.DataFrame({f'model{index}': models, f'score{index}': avg_scores})
        df = df.sort_values(by=[f'score{index}'], ascending=False)
        all_df_avg += [df]

    summary_df = pd.concat(all_dfs, axis=1)
    summary_df.to_csv(f"{args.workdir}/summary_{target}.csv")
    summary_df.head(20).to_csv(f"{args.workdir}/summary_20_{target}.csv")
    print(summary_df)