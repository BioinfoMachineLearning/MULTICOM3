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
    parser.add_argument('--targetname', type=str, required=True)
    parser.add_argument('--refpdb', type=str, required=True)
    parser.add_argument('--refpdb_A', type=str, required=True)
    parser.add_argument('--refpdb_ori', type=str, required=True)
    args = parser.parse_args()

    chain_ids = [['A', 'B', 'C', 'D'],
                 ['E', 'F', 'G', 'H'],
                 ['I', 'J', 'K', 'L']]

    makedir_if_not_exists(args.outdir)

    models = []
    tmscores_C = []
    tmscores_ori = []

    tmscores = {'tmscore1': [],
                'tmscore2': [],
                'tmscore3': []}

    matched_orders = {'matched_order1': [],
                      'matched_order2': [],
                      'matched_order3': []}

    same_seqs = {'same_seq1': [],
                 'same_seq2': [],
                 'same_seq3': []}

    for inpdb in os.listdir(args.indir):
        if not inpdb.startswith(args.targetname) or not inpdb.endswith('.pdb'):
            continue
        chain_contents = {}
        for line in open(args.indir + '/' + inpdb):
            if line.startswith('ATOM'):
                chain_id = line[21]
                if chain_id not in chain_contents:
                    chain_contents[chain_id] = [line]
                else:
                    chain_contents[chain_id] += [line]

        print(''.join([chain_id for chain_id in chain_contents]))
        pdb_outdir = args.outdir + '/' + inpdb.replace('.pdb', '')
        makedir_if_not_exists(pdb_outdir)
        matched_order, tmscore, src_seq, trg_seq = cal_mmalign(args.indir + '/' + inpdb,
                                                               args.refpdb_A, f"{pdb_outdir}/A9.txt")
        tmscores_C += [tmscore]

        matched_order, tmscore, src_seq, trg_seq = cal_mmalign(args.indir + '/' + inpdb,
                                                               args.refpdb_ori, f"{pdb_outdir}/ori.txt")
        tmscores_ori += [tmscore]

        print(tmscore)

        for i in range(len(chain_ids)):
            with open(f"{pdb_outdir}/{i + 1}.pdb", 'w') as fw:
                for chain in chain_ids[i]:
                    for line in chain_contents[chain]:
                        fw.write(line[:21] + chain + line[22:])
                    fw.write("TER\n")
            matched_order, tmscore, src_seq, trg_seq = cal_mmalign(f"{pdb_outdir}/{i + 1}.pdb",
                                                                   args.refpdb, f"{pdb_outdir}/{i + 1}.txt")
            print(matched_order)
            print(src_seq.replace('-', '').rstrip('\n'))
            print(src_seq.replace('-', '') == trg_seq.replace('-', ''))
            tmscores[f'tmscore{i + 1}'].append(tmscore)
            matched_orders[f'matched_order{i+1}'].append(matched_order)
            same_seqs[f'same_seq{i+1}'].append(src_seq.replace('-', '') == trg_seq.replace('-', ''))

        models += [inpdb]

    df = pd.DataFrame({'model': models, 'tmscore_A': tmscores_C, 'tmscore_ori': tmscores_ori,
                       'tmscore1': tmscores['tmscore1'],
                       'matched_order1': matched_orders['matched_order1'],
                       'same_seq1': same_seqs['same_seq1'],

                       'tmscore2': tmscores['tmscore2'],
                       'matched_order2': matched_orders['matched_order2'],
                       'same_seq2': same_seqs['same_seq2'],

                       'tmscore3': tmscores['tmscore3'],
                       'matched_order3': matched_orders['matched_order3'],
                       'same_seq3': same_seqs['same_seq3'],
                       })
    print(df)
    df.to_csv(f'{args.outdir}/summary.csv')
