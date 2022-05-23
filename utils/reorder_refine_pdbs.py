import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.complex_templates_search import parsers
import json

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    makedir_if_not_exists(args.outdir)

    ranking_file = args.indir + '/final_ranking.csv'
    ranking_df = pd.read_csv(ranking_file)
    ranking_confidences = {}
    ranked_order = []
    model_count = 0
    for i in range(len(ranking_df)):
        model_name = ranking_df.loc[i, 'model']
        if not model_name.endswith('_ref.pdb'):
            continue

        if not os.path.exists(f"{args.indir}/{model_name}"):
            raise Exception(f"Cannot find {args.indir}/{model_name}")

        trg_model_name = f"ranked_{model_count}"

        result_output_path = f"{args.indir}/{model_name.replace('.pdb', '.pkl')}"
        with open(result_output_path, 'rb') as f:
            prediction_result = pickle.load(f)
            ranking_confidences[trg_model_name] = prediction_result['ranking_confidence']
            ranked_order.append(trg_model_name)
        os.system(f"cp {args.indir}/{model_name} {args.outdir}/{trg_model_name}.pdb")
        os.system(f"cp {result_output_path} {args.outdir}/result_{trg_model_name}.pkl")
        model_count += 1

    ranking_output_path = os.path.join(args.outdir, 'ranking_debug.json')
    with open(ranking_output_path, 'w') as f:
        label = 'iptm+ptm' if 'iptm' in prediction_result else 'plddts'
        f.write(json.dumps(
            {label: ranking_confidences, 'order': ranked_order}, indent=4))

    msadir = args.outdir + '/msas'
    makedir_if_not_exists(msadir)
    os.system(f"cp {args.indir}/{model_name.replace('.pdb', '.pkl')} {msadir}/monomer_final.a3m")




