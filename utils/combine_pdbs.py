import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir


def combine_pdb(indir, outdir):
    pdbs = {}
    for pdb in os.listdir(indir):
        if pdb.find('.atom') > 0:
            if pdb[:4] not in pdbs:
                pdbs[pdb[:4]] = [pdb]
            else:
                pdbs[pdb[:4]] += [pdb]

    for pdbcode in pdbs:
        chain_pdbs = pdbs[pdbcode]
        chain_pdbs.sort()
        with open(outdir + '/' + pdbcode + '.atom', 'a') as fw:
            for chain_pdb in chain_pdbs:
                last_atom_num = -1
                last_res_num = -1
                last_res_name = ""
                final_line = ""
                for line in open(indir + '/' + chain_pdb):
                    chain_id = chain_pdb[4]
                    if line.startswith('ATOM'):
                        fw.write(line[:21] + chain_id + line[22:])
                        last_atom_num = line[6:11]
                        last_res_name = line[17:20]
                        last_res_num = line[22:26]
                        final_line = line
                fw.write("TER\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_file, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    combine_pdb(args.indir, args.outdir)
