import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir
from bml_casp15.quaternary_structure_evaluation.pairwise_mmalign import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True)
    parser.add_argument('--mmalign_program', type=is_file, required=True)

    args = parser.parse_args()

    pairwise_ranking = Pairwise_MMalign_qa(args.mmalign_program).run(args.indir)
    pairwise_ranking.to_csv('pairwise_ranking.csv')
