import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.monomer_structure_evaluation.bfactor_ranking import Bfactor_qa


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_file, required=True)
    args = parser.parse_args()
    bfactorqa = Bfactor_qa()
    bfactor_ranking = bfactorqa.run(input_dir=args.indir)
    bfactor_ranking.to_csv(args.outdir + '/bfactor_ranking.csv')


