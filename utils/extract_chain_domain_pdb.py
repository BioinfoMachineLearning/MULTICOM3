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


def extract_pkl(inpkl, outpkl, indices):
    with open(inpkl, 'rb') as f:
        prediction_result = pickle.load(f)
        prediction_result_new = {'plddt': prediction_result['plddt'][indices,]}
        if 'distogram' in prediction_result:
            distogram_new = prediction_result['distogram']
            distogram_new['logits'] = distogram_new['logits'][indices, indices, :]
            prediction_result_new['distogram'] = distogram_new
        with open(outpkl, 'wb') as f:
            pickle.dump(prediction_result_new, f, protocol=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('--fasta_path', type=str, required=True)
    parser.add_argument('--inpdb', type=str, required=True)
    parser.add_argument('--outpdb', type=str, required=True)
    parser.add_argument('--domain_infos', type=str, required=True)
    args = parser.parse_args()

    chain_domains = {}
    for chain_domain_info in args.domain_infos.split(','):
        chain, domain_start, domain_end = chain_domain_info.split('_')
        chain_domains[chain] = dict(start=domain_start,end=domain_end)

    for inpdb in [args.inpdb]:
        if inpdb.find('.pdb') < 0:
            continue
        chain_contents = {}
        prev_chain_name = ""
        residue_count = 0
        prev_num = 0
        chain_domain_indices = []
        for line in open(inpdb, 'r').readlines():
            if not line.startswith('ATOM'):
                continue
            residue_num = int(line[22:26].strip())
            if residue_num != prev_num:
                residue_count += 1
                prev_num = residue_num

            chain_id = line[21]
            if chain_id not in chain_domains:
                continue

            if chain_id not in chain_contents:
                if chain_id != prev_chain_name and prev_chain_name != "":
                    chain_contents[prev_chain_name]['res_end'] = residue_count - 1
                chain_contents[chain_id] = dict(content=[],
                                                sequence="",
                                                res_start=residue_count,
                                                res_end=-1)
                prev_chain_name = chain_id

            if int(chain_domains[chain_id]['start']) <= residue_num < int(chain_domains[chain_id]['end']):
                chain_contents[chain_id]['content'] += [line]
                if residue_count - 1 not in chain_domain_indices:
                    chain_domain_indices += [residue_count - 1]


        with open(args.outpdb, 'w') as fw:
            for chain_id in chain_contents:
                if len(chain_contents[chain_id]['content']) <= 1:
                    continue
                fw.writelines(chain_contents[chain_id]['content'])
                fw.write('TER\n')
            fw.write('END')

        print(chain_domain_indices)
