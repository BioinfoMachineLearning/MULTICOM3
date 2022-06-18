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


def cal_tmscore(inpdb, nativepdb, outfile):
    if not os.path.exists(outfile):
        cmd = f"/home/bml_casp15/BML_CASP15/tools/TMscore {inpdb} {nativepdb} > {outfile}"
        os.system(cmd)
        print(cmd)
    tmscore_contents = open(outfile).readlines()
    common_res = 0
    tmscore = ""
    for i, line in enumerate(tmscore_contents):
        line = line.rstrip('\n')
        if line.startswith('Number of residues in common='):
            common_res = int(line.split('=')[1])
            #print(common_res)
        if line.startswith('TM-score'):
            tmscore = line.split('=')[1].split()[0]
    return common_res, tmscore


def cal_mmalign(inpdb, nativepdb, outfile):
    if not os.path.exists(outfile):
        cmd = f"/home/bml_casp15/BML_CASP15/tools/MMalign {inpdb} {nativepdb} > {outfile}"
        os.system(cmd)
        print(cmd)
    tmscore_contents = open(outfile).readlines()
    matched_order_1 = ""
    matched_order_2 = ""
    aligned_length = ""
    seqid = ""
    tmscore = ""
    for i, line in enumerate(tmscore_contents):
        line = line.rstrip('\n')
        if line.startswith('Name of Chain_1:'):
            tmp = line.split(':', 2)
            matched_order_1 = tmp[2]
            matched_order_1 = matched_order_1.split()[0]
        if line.startswith('Name of Chain_2:'):
            tmp = line.split(':', 2)
            matched_order_2 = tmp[2]
            matched_order_2 = matched_order_2.split()[0]
        if line.startswith('TM-score='):
            tmscore = line.split()[1]
        # if line.startswith('(":" denotes residue'):
        #     src_seq = tmscore_contents[i + 1]
        #     trg_seq = tmscore_contents[i + 3]
        #     break
        if line.startswith('Aligned length'):
            aligned_length = line.split(',')[0].split()[2]
            seqid = line.split(',')[2].split('=')[2]

    return matched_order_1 + ' - ' + matched_order_2, tmscore, aligned_length, seqid


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--fasta_path', type=str, required=True)
    parser.add_argument('--targetname', type=str, required=True)
    parser.add_argument('--refpdb_10', type=str, required=True)
    parser.add_argument('--refpdb_6', type=str, required=True)
    parser.add_argument('--refpdb_6_fixed', type=str, required=False)
    args = parser.parse_args()

    makedir_if_not_exists(args.outdir)

    with open(args.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    models = []
    tmscores_10 = []
    matched_orders_10 = []
    Aligned_lengths_10 = []
    seq_ids_10 = []

    tmscores_6 = []
    matched_orders_6 = []
    Aligned_lengths_6 = []
    seq_ids_6 = []

    tmscores_6_fixed = []
    matched_orders_6_fixed = []
    Aligned_lengths_6_fixed = []
    seq_ids_6_fixed = []

    clash_counts = []

    chain_tmscores = {}

    chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']

    for inpdb in os.listdir(args.indir):
        if not inpdb.startswith(args.targetname) or not inpdb.endswith('.pdb'):
            continue

        pdb_outdir = args.outdir + '/' + inpdb.replace('.pdb', '')
        makedir_if_not_exists(pdb_outdir)
        matched_order, tmscore, aligned_length, seq_id = cal_mmalign(args.indir + '/' + inpdb,
                                                                     args.refpdb_10, f"{pdb_outdir}/tenmer.tmscore")
        tmscores_10 += [tmscore]
        matched_orders_10 += [matched_order]
        Aligned_lengths_10 += [aligned_length]
        seq_ids_10 += [seq_id]
        print(matched_order)

        matched_order, tmscore, aligned_length, seq_id = cal_mmalign(args.indir + '/' + inpdb,
                                                                     args.refpdb_6, f"{pdb_outdir}/sixmer.txt")
        tmscores_6 += [tmscore]
        matched_orders_6 += [matched_order]
        Aligned_lengths_6 += [aligned_length]
        seq_ids_6 += [seq_id]

        if args.refpdb_6_fixed is not None:
            matched_order, tmscore, aligned_length, seq_id = cal_mmalign(args.indir + '/' + inpdb,
                                                                         args.refpdb_6_fixed,
                                                                         f"{pdb_outdir}/sixmer_fixed.txt")
            tmscores_6_fixed += [tmscore]
            matched_orders_6_fixed += [matched_order]
            Aligned_lengths_6_fixed += [aligned_length]
            seq_ids_6_fixed += [seq_id]
        else:
            tmscores_6_fixed += []
            matched_orders_6_fixed += []
            Aligned_lengths_6_fixed += []
            seq_ids_6_fixed += []

        pdb_outdir_10 = pdb_outdir + '/tenmer_split'
        makedir_if_not_exists(pdb_outdir_10)

        os.system(
            f"python /home/bml_casp15/BML_CASP15/utils/split_pdbs.py --inpdb {args.refpdb_10} --outdir {pdb_outdir_10}")

        pdb_outdir_6 = pdb_outdir + '/sixmer_split'
        makedir_if_not_exists(pdb_outdir_6)

        if args.refpdb_6_fixed is not None:
            os.system(
                f"python /home/bml_casp15/BML_CASP15/utils/split_pdbs.py --inpdb {args.refpdb_6_fixed} --outdir {pdb_outdir_6}")
        else:
            os.system(
                f"python /home/bml_casp15/BML_CASP15/utils/split_pdbs.py --inpdb {args.refpdb_6} --outdir {pdb_outdir_6}")

        pdb_outdir_final = pdb_outdir + '/final'
        os.system(
            f"python /home/bml_casp15/BML_CASP15/utils/split_pdbs.py --inpdb {args.indir}/{inpdb} --outdir {pdb_outdir_final}")

        refdir = pdb_outdir_6
        for chain in chains:
            if chain == 'G':
                refdir = pdb_outdir_10
            chain_pdb = pdb_outdir_final + '/' + chain + '.pdb'
            chain_pdb_reindex = chain_pdb + '.reindex'
            os.system(f"perl /home/bml_casp15/BML_CASP15/utils/reindex_pdb.pl {chain_pdb} {chain_pdb_reindex}")
            common_res, tmscore = cal_tmscore(chain_pdb_reindex, refdir + '/' + chain + '.pdb',
                                              pdb_outdir_final + '/' + chain + '.txt')
            if common_res != len(chain_id_map[chain].sequence):
                raise Exception('The sequence length does not match!')

            if chain not in chain_tmscores:
                chain_tmscores[chain] = [tmscore]
            else:
                chain_tmscores[chain] += [tmscore]

        os.system(
            f"perl /home/bml_casp15/BML_CASP15/utils/clash_check.pl {args.indir}/{inpdb} &> {pdb_outdir}/total.clash")
        total_clash_count = len(open(f"{pdb_outdir}/total.clash").readlines())
        clash_counts += [total_clash_count]
        models += [inpdb]

    df = pd.DataFrame({'model': models,
                       'clash_count': clash_counts,
                       'tmscore_10': tmscores_10,
                       'matched_order_10': matched_orders_10,
                       'aligned_length_10': Aligned_lengths_10,
                       'seqid_10': seq_ids_10,

                       'tmscore_6': tmscores_6,
                       'matched_order_6': matched_orders_6,
                       'aligned_length_6': Aligned_lengths_6,
                       'seqid_6': seq_ids_6,

                       'tmscore_6_fixed': tmscores_6_fixed,
                       'matched_order_6_fixed': matched_orders_6_fixed,
                       'aligned_length_6_fixed': Aligned_lengths_6_fixed,
                       'seqid_6_fixed': seq_ids_6_fixed,

                       'A': chain_tmscores['A'],
                       'B': chain_tmscores['B'],
                       'C': chain_tmscores['C'],
                       'D': chain_tmscores['D'],
                       'E': chain_tmscores['E'],
                       'F': chain_tmscores['F'],
                       'G': chain_tmscores['G'],
                       'H': chain_tmscores['H'],
                       'I': chain_tmscores['I'],
                       'J': chain_tmscores['J'],
                       })
    print(df)
    df.to_csv(f'{args.outdir}/summary.csv')
