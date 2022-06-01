import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.quaternary_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_file, required=True)
    args = parser.parse_args()
    alphafold_qa = Alphafold_pkl_qa()
    alphafold_ranking = alphafold_qa.run(args.indir)
    alphafold_ranking.to_csv(args.outdir + '/alphafold_ranking.csv')

