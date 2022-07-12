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

def cal_tmscore(tmscore_program, inpdb, nativepdb, tmpdir):
    cwd = os.getcwd()

    makedir_if_not_exists(tmpdir)

    os.chdir(tmpdir)

    os.system(f"cp {inpdb} inpdb.pdb")
    os.system(f"cp {nativepdb} native.pdb")

    inpdb = "inpdb.pdb"
    nativepdb = "native.pdb"

    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep TM-score | awk '{print $3}' "
    print(cmd)
    tmscore_contents = os.popen(cmd).read().split('\n')
    tmscore = float(tmscore_contents[2].rstrip('\n'))
    cmd = tmscore_program + ' ' + inpdb + ' ' + nativepdb + " | grep GDT-score | awk '{print $3}' "
    tmscore_contents = os.popen(cmd).read().split('\n')
    gdtscore = float(tmscore_contents[0].rstrip('\n'))

    os.chdir(cwd)

    # os.system("rm -rf " + tmpdir)

    return tmscore, gdtscore


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--fasta_path', type=str, required=True)
    parser.add_argument('--tmscore_program', type=str, required=True)
    parser.add_argument('--tmpdir', type=str, required=True)
    args = parser.parse_args()

    models = []
    tmscores = {'tmscore1': [],
                'tmscore2': [],
                'tmscore3': []}
    same_seqs = []

    with open(args.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    pdbseq = '/bmlfast/bml_casp15/BML_CASP15//utils/pdb2seq.pl'

    for topdir in os.listdir(args.indir):

        finalpdb = args.indir + '/' + topdir + '/T1169.pdb'

        final_region_pdb_1 = args.indir + '/' + topdir + '/T1169r1.pdb'

        os.system(f'perl /home/bml_casp15/BML_CASP15/utils/extract_domain.pl {finalpdb} {final_region_pdb_1} 1 380')

        regionpdb1 = args.indir + '/' + topdir + '/region1.pdb'

        tmscore, gdtscore = cal_tmscore(args.tmscore_program, final_region_pdb_1, regionpdb1, args.tmpdir)

        tmscores['tmscore1'] += [tmscore]

        regionpdb2 = args.indir + '/' + topdir + '/region2.pdb'

        final_region_pdb_2 = args.indir + '/' + topdir + '/T1169r2.pdb'

        os.system(f'perl /home/bml_casp15/BML_CASP15/utils/extract_domain.pl {finalpdb} {final_region_pdb_2} 381 2952')

        tmscore, gdtscore = cal_tmscore(args.tmscore_program, final_region_pdb_2, regionpdb2, args.tmpdir)

        tmscores['tmscore2'] += [tmscore]

        regionpdb3 = args.indir + '/' + topdir + '/region3.pdb'

        final_region_pdb_3 = args.indir + '/' + topdir + '/T1169r3.pdb'

        os.system(f'perl /home/bml_casp15/BML_CASP15/utils/extract_domain.pl {finalpdb} {final_region_pdb_3} 2944 10000')

        tmscore, gdtscore = cal_tmscore(args.tmscore_program, final_region_pdb_3, regionpdb3, args.tmpdir)

        tmscores['tmscore3'] += [tmscore]

        seq = os.popen(f'perl {pdbseq} {finalpdb}').readlines()[0].strip()
        print(seq)
        same_seqs += [seq == input_seqs[0]]

        models += [topdir]

    df = pd.DataFrame({'model': models,
                       'tmscore1': tmscores['tmscore1'],
                       'tmscore2': tmscores['tmscore2'],
                       'tmscore3': tmscores['tmscore3'],
                       'same_seq': same_seqs
                       })
    print(df)
    df.to_csv(f'summary.csv')
