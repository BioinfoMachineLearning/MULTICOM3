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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--start', type=int, required=True)
    parser.add_argument('--end', type=int, required=True)
    args = parser.parse_args()

    makedir_if_not_exists(args.outdir)

    pdbdir = args.indir + '/pdb_monomer'
    pkldir = args.indir + '/pkl_monomer'

    pdbdir_ori = args.outdir + '/pdb_ext'
    makedir_if_not_exists(pdbdir_ori)
    outpdbdir = args.outdir + '/pdb_monomer'
    makedir_if_not_exists(outpdbdir)
    outpkldir = args.outdir + '/pkl_monomer'
    makedir_if_not_exists(outpkldir)

    extract_script = '/home/bml_casp15/BML_CASP15/utils/extract_domain.pl'
    reindex_script = '/home/bml_casp15/BML_CASP15/utils/reindex_pdb.pl'

    for pdb in os.listdir(pdbdir):
        os.system(f"perl {extract_script} {pdbdir}/{pdb} {pdbdir_ori}/{pdb} {args.start} {args.end} ")
        os.system(f"perl {reindex_script} {pdbdir_ori}/{pdb} {outpdbdir}/{pdb}")
        extract_pkl(src_pkl=f"{pkldir}/{pdb.replace('.pdb', '.pkl')}",
                    output_pkl=f"{outpkldir}/{pdb.replace('.pdb', '.pkl')}",
                    residue_start=args.start-1,
                    residue_end=args.end-1)

