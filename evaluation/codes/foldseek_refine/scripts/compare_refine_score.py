import argparse, os
import sys
import logging
import pandas as pd 
import numpy as np

def cal_mmalign(inpdb, nativepdb, outfile):
    if not os.path.exists(outfile):
        raise Exception(f"cannot find {outfile}")
    tmscore = 0
    for line in open(outfile):
        if line.find('TM-score') == 0 and line.find('Structure_1') > 0:
            tmscore = float(line.split()[1])
    return tmscore

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-n", "--nativedir", help="output folder for the pdb files", type=str, required=True)
    parser.add_argument("-w", "--workdir", help="output folder for the pdb files", type=str, required=True)
    parser.add_argument("-t", "--tmpdir", help="output folder for the pdb files", type=str, required=True)

    args = parser.parse_args()

    data_all_dict = {'targetname': [], 'ori_plddt': [], 'ori_tmscore': [], 'ref_plddt': [], 'ref_tmscore': []}
    data_inc_dict = {'targetname': [], 'ori_plddt': [], 'ori_tmscore': [], 'ref_plddt': [], 'ref_tmscore': []}
    data_all_avg_dict = {'targetname': [], 'ori_tmscore': [], 'ref_tmscore': []}
    data_inc_avg_dict = {'targetname': [], 'ori_tmscore': [], 'ref_tmscore': []}
    
    data_top1_dict = {'targetname': [], 'ori_tmscore': [], 'ref_tmscore': []}
    data_max_dict = {'targetname': [], 'ori_tmscore': [], 'ref_tmscore': []}
    
    for targetname in sorted(os.listdir(args.workdir)):
        tmpdir = args.tmpdir + '/' + targetname

        nativepdb = args.nativedir + '/' + targetname + '_filtered.pdb'
        if not os.path.exists(nativepdb):
            continue
        final_ranking = pd.read_csv(args.workdir + '/' + targetname + '/final_ranking.csv')
        model_names = set([model_name.rstrip('_ref.pdb').rstrip('_ori.pdb') for model_name in final_ranking['model'] if model_name.find('deep') < 0])

        ori_tmscores = []
        ref_tmscores = []

        ori_tmscores_inc = []
        ref_tmscores_inc = []
        for model_name in model_names:
            ori_plddt, ori_tmscore = 0, 0
            ref_plddt, ref_tmscore = 0, 0
            for i in range(len(final_ranking)):
                if final_ranking.loc[i, 'model'] == model_name + '_ori.pdb':
                    ori_plddt = final_ranking.loc[i, 'confidence']
                    ori_tmscore = cal_mmalign( inpdb=f"{args.workdir}/{targetname}/{model_name}_ori.pdb", nativepdb=nativepdb, outfile=tmpdir + '/' + model_name + '_ori_filtered_out')
                    ori_tmscores += [ori_tmscore]
                elif final_ranking.loc[i, 'model'] == model_name + '_ref.pdb':
                    ref_plddt = final_ranking.loc[i, 'confidence']
                    ref_tmscore = cal_mmalign( inpdb=f"{args.workdir}/{targetname}/{model_name}_ref.pdb", nativepdb=nativepdb, outfile=tmpdir + '/' + model_name + '_ref_filtered_out')
                    ref_tmscores += [ref_tmscore]

            data_all_dict['targetname'] += [targetname]
            data_all_dict['ori_plddt'] += [ori_plddt]
            data_all_dict['ori_tmscore'] += [ori_tmscore]
            data_all_dict['ref_plddt'] += [ref_plddt]
            data_all_dict['ref_tmscore'] += [ref_tmscore]

            if ref_plddt > ori_plddt:
                data_inc_dict['targetname'] += [targetname]
                data_inc_dict['ori_plddt'] += [ori_plddt]
                data_inc_dict['ori_tmscore'] += [ori_tmscore]
                data_inc_dict['ref_plddt'] += [ref_plddt]
                data_inc_dict['ref_tmscore'] += [ref_tmscore]

                ori_tmscores_inc += [ori_tmscore]
                ref_tmscores_inc += [ref_tmscore]

        data_top1_dict['targetname'] += [targetname]
        avg_ranking = pd.read_csv(args.workdir + '/' + targetname + '/pairwise_af_avg.ranking')
        top1_model = avg_ranking.loc[0, 'model'].replace('.pdb', '_ori.pdb')
        ori_tmscore = cal_mmalign( inpdb=f"{args.workdir}/{targetname}/{top1_model}", nativepdb=nativepdb, outfile=tmpdir + '/' + top1_model.rstrip('.pdb') + '_filtered_out')
        data_top1_dict['ori_tmscore'] += [ori_tmscore]
                
        refine_top1_model = final_ranking.loc[0, 'model']
        ref_tmscore = cal_mmalign( inpdb=f"{args.workdir}/{targetname}/{refine_top1_model}", nativepdb=nativepdb, outfile=tmpdir + '/' + refine_top1_model.rstrip('.pdb') + '_filtered_out')
        data_top1_dict['ref_tmscore'] += [ref_tmscore]

        data_max_dict['targetname'] += [targetname]
        data_max_dict['ori_tmscore'] += [np.max(np.array(ori_tmscores))]

        data_max_dict['ref_tmscore'] += [np.max(np.array(ref_tmscores))]

        data_all_avg_dict['targetname'] += [targetname]
        data_all_avg_dict['ori_tmscore'] += [np.mean(np.array(ori_tmscores))]

        data_all_avg_dict['ref_tmscore'] += [np.mean(np.array(ref_tmscores))]

        if len(ori_tmscores_inc) > 0:
            data_inc_avg_dict['targetname'] += [targetname]
            data_inc_avg_dict['ori_tmscore'] += [np.mean(np.array(ori_tmscores_inc))]
            data_inc_avg_dict['ref_tmscore'] += [np.mean(np.array(ref_tmscores_inc))]

    df_all = pd.DataFrame(data_all_dict)
    df_all.to_csv('refine_all.csv')   

    df_inc_all = pd.DataFrame(data_inc_dict)
    df_inc_all.to_csv('refine_inc_all.csv')  

    df_avg = pd.DataFrame(data_all_avg_dict)
    df_avg.to_csv('refine_all_avg.csv')   

    df_inc_avg = pd.DataFrame(data_inc_avg_dict)
    df_inc_avg.to_csv('refine_inc_avg.csv') 

    df_top1 = pd.DataFrame(data_top1_dict)
    df_top1.to_csv('refine_top1.csv')   

    df_max = pd.DataFrame(data_max_dict)
    df_max.to_csv('refine_max.csv')   
