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

def make_msa_features(msas, msa_output_dir, msa_save_path):
    """Constructs a feature dict of MSA features."""
    if not msas:
        raise ValueError('At least one MSA must be provided.')

    seen_desc = []
    seen_sequences = []
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(f'MSA {msa_index} must contain at least one sequence.')
        for sequence_index, sequence in enumerate(msa.sequences):
            if filter and sequence in seen_sequences:
                continue
            seen_sequences += [sequence]
            seen_desc += [msa.descriptions[sequence_index]]

    with open(msa_output_dir + '/' + msa_save_path, 'w') as fw:
        for (desc, seq) in zip(seen_desc, seen_sequences):
            fw.write(f'>{desc}\n{seq}\n')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    makedir_if_not_exists(args.outdir)
    ranking_confidences = {}
    ranked_order = []
    for i in range(5):
        model_name = f"ranked_{i}"
        if not os.path.exists(f"{args.indir}/{model_name}.pdb"):
            raise Exception(f"Cannot find {args.indir}/{model_name}.pdb")
        result_output_path = f"{args.indir}/{model_name}.pkl"
        with open(result_output_path, 'rb') as f:
            prediction_result = pickle.load(f)
            ranking_confidences[model_name] = prediction_result['ranking_confidence']
            ranked_order.append(model_name)
        os.system(f"cp {args.indir}/{model_name}.pdb {args.outdir}/{model_name}.pdb")
        os.system(f"cp {args.indir}/{model_name}.pkl {args.outdir}/result_{model_name}.pkl")

    ranking_output_path = os.path.join(args.outdir, 'ranking_debug.json')
    with open(ranking_output_path, 'w') as f:
        label = 'iptm+ptm' if 'iptm' in prediction_result else 'plddts'
        f.write(json.dumps(
            {label: ranking_confidences, 'order': ranked_order}, indent=4))

    msadir = args.outdir + '/msas'
    makedir_if_not_exists(msadir)

    os.system(f"cp {args.indir}/bfd_uniclust_hits.a3m {msadir}/bfd_uniclust_hits.a3m")
    os.system(f"cp {args.indir}/mgnify_hits.sto {msadir}/mgnify_hits.sto")
    os.system(f"cp {args.indir}/pdb_hits.hhr {msadir}/pdb_hits.hhr")
    os.system(f"cp {args.indir}/uniref90_hits.sto {msadir}/uniref90_hits.sto")

    uniref90_msa = parsers.parse_stockholm(open(f"{msadir}/uniref90_hits.sto").read())
    uniref90_msa = uniref90_msa.truncate(max_seqs=10000)
    bfd_msa = parsers.parse_a3m(open(f"{msadir}/bfd_uniclust_hits.a3m").read())
    mgnify_msa = parsers.parse_stockholm(open(f"{msadir}/mgnify_hits.sto").read())
    mgnify_msa = mgnify_msa.truncate(max_seqs=501)
    make_msa_features((uniref90_msa, bfd_msa, mgnify_msa), msadir, msa_save_path='monomer_final.a3m')




