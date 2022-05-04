import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir
from bml_casp15.monomer_structure_refinement.util import cal_tmscore

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True)
    parser.add_argument('--native_pdb', type=is_file, required=True)
    parser.add_argument('--tmscore_program', type=is_file, required=True)
    parser.add_argument('--tmpdir', type=is_file, required=True)

    args = parser.parse_args()

    models = []
    tmscores = []
    gdtscores = []
    for model in os.listdir(args.indir):
        if model.find('.pdb') > 0:
            ref_tmscore, ref_gdtscore = cal_tmscore(args.tmscore_program, args.indir + '/' + model, args.native_pdb,
                                                    args.tmpdir)
            models += [model]
            tmscores += [ref_tmscore]
            gdtscores += [ref_gdtscore]

    df = pd.DataFrame({'model': models, 'tmscore': tmscores, 'gdtscore': gdtscores})
    df.sort_values(by=['gdtscore'], ascending=False)
    df.to_csv('result.csv')
    print(df)

