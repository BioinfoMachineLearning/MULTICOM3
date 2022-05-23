import os, sys, argparse, time
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import write_fasta
from bml_casp15.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from bml_casp15.monomer_structure_evaluation.pipeline_sep import get_sequence
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import pickle
import pathlib

# need to add A if using relaxation in alphafold
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'


def combine_pdb(pdbs, combine_pdb):
    with open(combine_pdb, 'w') as out:
        for chain_id, pdb in zip(PDB_CHAIN_IDS, pdbs):
            for line in open(pdb):
                if not line.startswith('ATOM'):
                    continue
                out.write(line[:21] + chain_id + line[22:])
            out.write('TER\n')


def extract_pkl(inpkl, outpkl, chain_starts, chain_ends):
    with open(inpkl, 'rb') as f:
        prediction_result = pickle.load(f)
        indices = []
        for chain_start, chain_end in zip(chain_starts, chain_ends):
            indices += list(range(chain_start, chain_end + 1))
        prediction_result_new = {'plddt': prediction_result['plddt'][indices,]}
        if 'distogram' in prediction_result:
            distogram_new = prediction_result['distogram']
            distogram_new['logits'] = distogram_new['logits'][indices, indices, :]
            prediction_result_new['distogram'] = distogram_new
        with open(outpkl, 'wb') as f:
            pickle.dump(prediction_result_new, f, protocol=4)


def extract_monomer_pdbs_reindex(complex_pdb, sequence, output_prefix):
    chain_contents = {}
    prev_chain_name = ""
    residue_count = 0
    prev_num = 0
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        residue_num = int(line[22:26].strip())
        if residue_num != prev_num:
            residue_count += 1
            prev_num = residue_num

        chain_id = line[21]
        if chain_id not in chain_contents:
            if chain_id != prev_chain_name and prev_chain_name != "":
                chain_contents[prev_chain_name]['res_end'] = residue_count - 1
                chain_contents[prev_chain_name]['sequence'] = get_sequence(chain_contents[prev_chain_name]['content'])
            chain_contents[chain_id] = dict(content=[line],
                                            sequence="",
                                            res_start=residue_count,
                                            res_end=-1)
            prev_chain_name = chain_id
        else:
            chain_contents[chain_id]['content'] += [line]

    if prev_chain_name != "":
        chain_contents[prev_chain_name]['res_end'] = residue_count
        chain_contents[prev_chain_name]['sequence'] = get_sequence(chain_contents[prev_chain_name]['content'])

    result_dict = {}
    for chain_id in chain_contents:
        if chain_contents[chain_id]['sequence'] == sequence:
            with open(f"{output_prefix}_{chain_id}.pdb_ori", 'w') as fw:
                fw.writelines(chain_contents[chain_id]['content'])
            os.system(f"perl /home/bml_casp15/BML_CASP15/utils/reindex_pdb.pl "
                      f"{output_prefix}_{chain_id}.pdb_ori {output_prefix}_{chain_id}.pdb")
            result_dict[chain_id] = dict(pdb=f"{output_prefix}_{chain_id}.pdb",
                                         pdbname=pathlib.Path(f"{output_prefix}_{chain_id}.pdb").name,
                                         chain_start=chain_contents[chain_id]['res_start'] - 1,
                                         chain_end=chain_contents[chain_id]['res_end'] - 1)
    return result_dict


def cal_contact_ratio_between_pdbs(pdb1, pdb2):
    parser = PDBParser(PERMISSIVE=1)
    structure1 = parser.get_structure('', pdb1)
    structure2 = parser.get_structure('', pdb2)
    print(pdb2)

    model1 = structure1[0]
    chain_id1 = list(model1.child_dict.keys())
    xyzPDB1 = model1[chain_id1[0]]

    model2 = structure2[0]
    chain_id2 = list(model2.child_dict.keys())
    xyzPDB2 = model2[chain_id2[0]]

    h_dist_map = np.zeros((len(xyzPDB1), len(xyzPDB2)))
    for i in range(1, len(xyzPDB1) + 1):
        for j in range(1, len(xyzPDB2) + 1):
            res_i = xyzPDB1[i]
            res_j = xyzPDB2[j]
            dist_list = []
            for atom_i in res_i:
                for atom_j in res_j:
                    if ('C' in atom_i.name or 'N' in atom_i.name or 'O' in atom_i.name or 'S' in atom_i.name) and \
                            ('C' in atom_j.name or 'N' in atom_j.name or 'O' in atom_j.name or 'S' in atom_j.name):
                        dist_list.append(atom_i - atom_j)
                    else:
                        continue
            min_dist = np.min(dist_list)
            h_dist_map[i - 1, j - 1] = min_dist

    return (h_dist_map < 6).sum() / (len(xyzPDB1) * len(xyzPDB2))


def pair_chain_ids(src_pair_dict, trg_chain_pdb_dict):
    chain_ids_1 = []
    chain_ids_2 = []
    ratios = []
    for src_pair_id in src_pair_dict:
        for trg_chain_id in trg_chain_pdb_dict:
            if trg_chain_id in src_pair_id.split('_'):
                continue
            max_ratio = 0
            pdb2 = trg_chain_pdb_dict[trg_chain_id]['pdb']
            for src_chain_idx, src_chain_id in enumerate(src_pair_id.split('_')):
                pdb1 = src_pair_dict[src_pair_id]['pdb'][src_chain_idx]
                ratio = cal_contact_ratio_between_pdbs(pdb1, pdb2)
                if ratio > max_ratio:
                    max_ratio = ratio
            chain_ids_1 += [src_pair_id]
            chain_ids_2 += [trg_chain_id]
            ratios += [max_ratio]

    df = pd.DataFrame({'chain1': chain_ids_1, 'chain2': chain_ids_2, 'ratio': ratios})
    df = df.sort_values(by=['ratio'], ascending=False)
    df.reset_index(inplace=True, drop=True)
    print(df)

    processed_ids = []
    res_pair_dict = {}
    for i in range(len(df)):
        chain1 = df.loc[i, 'chain1']
        chain2 = df.loc[i, 'chain2']
        if chain1 in processed_ids or chain2 in processed_ids:
            continue

        res_pair_dict[f"{chain1}_{chain2}"] = {'pdb': [], 'chain_start': [], 'chain_end': []}
        res_pair_dict[f"{chain1}_{chain2}"]['pdb'] += src_pair_dict[chain1]['pdb']
        res_pair_dict[f"{chain1}_{chain2}"]['pdb'] += [trg_chain_pdb_dict[chain2]['pdb']]

        res_pair_dict[f"{chain1}_{chain2}"]['chain_start'] += src_pair_dict[chain1]['chain_start']
        res_pair_dict[f"{chain1}_{chain2}"]['chain_start'] += [trg_chain_pdb_dict[chain2]['chain_start']]

        res_pair_dict[f"{chain1}_{chain2}"]['chain_end'] += src_pair_dict[chain1]['chain_end']
        res_pair_dict[f"{chain1}_{chain2}"]['chain_end'] += [trg_chain_pdb_dict[chain2]['chain_end']]

        processed_ids += [chain1]
        processed_ids += [chain2]

    return res_pair_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_path', type=str, required=True)
    parser.add_argument('--inpdb', type=str, required=True)
    parser.add_argument('--inpkl', type=str, required=False)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()

    makedir_if_not_exists(args.outdir)

    with open(args.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    # ABCDEFGHI
    # A3B3C3 -> 3 pairs of ABC

    # 1. calculate the contact ratio of each pairs: AB, AC, AD, AE .......
    # 2. merge pairs from different sequences based on the maximum contact ratio:
    # AD, BE, CF, AG,
    # DG, EH, FI,

    workdir = args.outdir + '/work'
    makedir_if_not_exists(workdir)
    multimer_pdbs = {}
    for chain_id in chain_id_map:
        chain_pdb_dict = extract_monomer_pdbs_reindex(complex_pdb=args.inpdb,
                                                      sequence=chain_id_map[chain_id].sequence,
                                                      output_prefix=f"{workdir}/{chain_id}")
        print(chain_pdb_dict)
        multimer_pdbs[chain_id] = chain_pdb_dict

    src_pair_dict = {}
    for chain_idx, chain_id in enumerate(multimer_pdbs):
        if chain_idx == 0:
            for same_chain_id in multimer_pdbs[chain_id]:
                src_pair_dict[same_chain_id] = {'pdb': [multimer_pdbs[chain_id][same_chain_id]['pdb']],
                                                'chain_start': [multimer_pdbs[chain_id][same_chain_id]['chain_start']],
                                                'chain_end': [multimer_pdbs[chain_id][same_chain_id]['chain_end']]}
        else:
            src_pair_dict = pair_chain_ids(src_pair_dict, multimer_pdbs[chain_id])

    print(src_pair_dict)

    for pair_idx, pair in enumerate(src_pair_dict):
        combine_pdb([src_pair_dict[pair]['pdb'][chain_idx] for chain_idx, chain_id in enumerate(pair.split('_'))],
                    args.outdir + '/' + pair + '.pdb')
        if args.inpkl is not None:
            extract_pkl(args.inpkl, args.outdir + '/' + pair + '.pkl',
                        src_pair_dict[pair]['chain_start'], src_pair_dict[pair]['chain_end'])
