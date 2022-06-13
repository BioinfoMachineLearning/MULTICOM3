import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.monomer_structure_evaluation.pipeline_sep import extract_pkl


def split_pdb(complex_pdb, outdir):
    pre_chain = None
    i = 0
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if pre_chain is None:
            pre_chain = chain_name
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line[:21] + ' ' + line[22:])
        elif chain_name == pre_chain:
            fw.write(line[:21] + ' ' + line[22:])
        else:
            fw.close()
            i = i + 1
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            fw.write(line[:21] + ' ' + line[22:])
            pre_chain = chain_name
    fw.close()


def combine_pdb(indir, outpdb):
    with open(outpdb, 'a') as fw:
        for chain_pdb in os.listdir(indir):
            if chain_pdb.find('_reindex') > 0:
                for line in open(indir + '/' + chain_pdb):
                    chain_id = chain_pdb[0]
                    if line.startswith('ATOM'):
                        fw.write(line[:21] + chain_id + line[22:])
                fw.write("TER\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--domain_infos', type=str, required=True)
    parser.add_argument('--pklstart', type=int, required=True)
    parser.add_argument('--pklend', type=int, required=True)
    args = parser.parse_args()

    makedir_if_not_exists(args.outdir)

    pdbdir = args.indir + '/pdb'
    pkldir = args.indir + '/pkl'

    pdboutdir = args.outdir + '/pdb'
    pkloutdir = args.outdir + '/pkl'
    makedir_if_not_exists(pdboutdir)
    makedir_if_not_exists(pkloutdir)

    workdir = args.outdir + '/work'
    pdbworkdir = workdir + '/pdb'
    pklworkdir = workdir + '/pkl'
    makedir_if_not_exists(pdbworkdir)
    makedir_if_not_exists(pklworkdir)

    extract_script = '/home/bml_casp15/BML_CASP15/utils/extract_domain.pl'
    reindex_script = '/home/bml_casp15/BML_CASP15/utils/reindex_pdb.pl'

    for pdb in os.listdir(pdbdir):
        outpdbworkdir = pdbworkdir + '/' + pdb
        makedir_if_not_exists(outpdbworkdir)
        split_pdb(pdbdir + '/' + pdb, outpdbworkdir)
        for domain_info in args.domain_infos.split(','):
            chain, chain_start, chain_end = domain_info.split('_')
            os.system(f"perl {extract_script} {outpdbworkdir}/{chain}.pdb {outpdbworkdir}/{chain}.pdb_ori "
                      f"{chain_start} {chain_end}")
            os.system(f"perl {reindex_script} {outpdbworkdir}/{chain}.pdb_ori {outpdbworkdir}/{chain}.pdb_reindex")

        combine_pdb(outpdbworkdir, pdboutdir + '/' + pdb)

        extract_pkl(src_pkl=f"{pkldir}/{pdb.replace('.pdb', '.pkl')}",
                    output_pkl=f"{pkloutdir}/{pdb.replace('.pdb', '.pkl')}",
                    residue_start=args.pklstart - 1,
                    residue_end=args.pklend - 1)
