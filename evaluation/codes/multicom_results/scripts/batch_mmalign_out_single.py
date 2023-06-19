import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir


def run_command(inparams):
    mmalign_program, pdb1, pdb2, outfile = inparams
    cmd = f"{mmalign_program} {pdb1} {pdb2} -TMscore 6 -ter 1 > {outfile} "
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True)
    parser.add_argument('--nativepdb', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_file, required=True)
    parser.add_argument('--mmalign_program', type=is_file, required=True)

    args = parser.parse_args()

    process_list = []

    for native_pdb in ["1"]:

        targetname = args.nativepdb.rstrip('_filtered.pdb')

        outdir = args.outdir

        makedir_if_not_exists(outdir)

        for model in os.listdir(args.indir):
            
            os.system(f"cp {args.indir}/{model} {outdir}/{model}_filtered")

            outfile = outdir + '/' + model + '_filtered_out'

            if os.path.exists(outfile) and len(open(outfile).readlines()) >= 10:
                continue

            process_list.append([args.mmalign_program, args.nativepdb, args.indir + '/' + model, outdir + '/' + model + '_filtered_out'])

    pool = Pool(processes=80)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()
