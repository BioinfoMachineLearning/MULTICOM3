import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir


def cal_mmalign(inpdb, nativepdb, outfile):
    cmd = f"/home/bml_casp15/BML_CASP15/tools/MMalign {inpdb} {nativepdb} > {outfile}"
    os.system(cmd)
    print(cmd)
    tmscore_contents = open(outfile).readlines()
    matched_order = ""
    tmscore = ""
    src_seq = ""
    trg_seq = ""
    for i, line in enumerate(tmscore_contents):
        line = line.rstrip('\n')
        if line.startswith('Name of Chain_2:'):
            tmp = line.split()
            matched_order = tmp[len(tmp) - 1]
            matched_order = matched_order.strip()
        if line.startswith('TM-score='):
            tmscore = line.split()[1]
        if line.startswith('(":" denotes residue'):
            src_seq = tmscore_contents[i + 1]
            trg_seq = tmscore_contents[i + 3]
            break
    return matched_order, tmscore, src_seq, trg_seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--refpdb', type=str, required=True)
    args = parser.parse_args()

    mapping = ['A_A,B_B,C_C,D_D,E_E', 'F_A,G_B,H_C,I_D,J_E', 'K_A,L_B,M_C,N_D,O_E', 'P_A,Q_B,R_C,S_D,T_E']

    makedir_if_not_exists(args.outdir)

    models = []
    tmscores_C = []

    tmscores = {'tmscore1': [],
                'tmscore2': [],
                'tmscore3': [],
                'tmscore4': []}

    clashes = {'clash1': [],
               'clash2': [],
               'clash3': [],
               'clash4': [],
               'final_clash': []}

    os.system(
        f"perl /home/bml_casp15/BML_CASP15/utils/clash_check.pl {args.refpdb} &> ref.clash")
    ref_clash_count = len(open('ref.clash').readlines())

    for inpdb in os.listdir(args.indir):
        pdb_outdir = args.outdir + '/' + inpdb
        makedir_if_not_exists(pdb_outdir)
        for i in range(len(mapping)):
            os.system(f"python /home/bml_casp15/BML_CASP15/utils/reorder_chain_ids.py "
                      f"--inpdb {args.indir}/{inpdb} --outpdb {pdb_outdir}/A1B2C2_{i + 1}.pdb --mapping {mapping[i]} ")

            matched_order, tmscore, src_seq, trg_seq = cal_mmalign(f"{pdb_outdir}/A1B2C2_{i + 1}.pdb",
                                                                   args.refpdb, f"{pdb_outdir}/A1B2C2_{i + 1}.txt")
            tmscores[f'tmscore{i + 1}'].append(tmscore)
            clash_file = f"{pdb_outdir}/A1B2C2_{i + 1}.clash"
            os.system(
                f"perl /home/bml_casp15/BML_CASP15/utils/clash_check.pl {pdb_outdir}/A1B2C2_{i + 1}.pdb &> {clash_file}")
            clash_count = len(open(clash_file).readlines())
            clashes[f'clash{i + 1}'].append(clash_count)

        models += [inpdb]
        os.system(
            f"perl /home/bml_casp15/BML_CASP15/utils/clash_check.pl {args.indir}/{inpdb} &> {pdb_outdir}/total.clash")
        total_clash_count = len(open(f"{pdb_outdir}/total.clash").readlines())
        clashes['final_clash'] += [total_clash_count]

    df = pd.DataFrame({'model': models,
                       'tmscore1': tmscores['tmscore1'],
                       'clash1': clashes['clash1'],

                       'tmscore2': tmscores['tmscore2'],
                       'clash2': clashes['clash2'],

                       'tmscore3': tmscores['tmscore3'],
                       'clash3': clashes['clash3'],

                       'tmscore4': tmscores['tmscore4'],
                       'clash4': clashes['clash4'],

                       'refclash': [ref_clash_count] * len(models),
                        'final_clash': clashes['final_clash']
                       })
    print(df)
    df.to_csv('summary.csv')
