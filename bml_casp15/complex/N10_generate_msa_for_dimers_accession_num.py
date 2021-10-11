###################################################################
# Script to generate dimers' msas using interaction calculated from RCSB_PDB
# Dependency: dimer's relation generated from RCSB_PDB
###################################################################

import os, sys, argparse, time, ftplib
from util import is_dir, is_file, read_option_file, check_dirs, check_contents, die, makedir_if_not_exists, clean_dir, \
    cal_sequence_identity_from_seq
from multiprocessing import Pool
from collections import OrderedDict
import re
import pandas as pd
from alignment import read_fasta, write_a3m, clean_a3m
from complex import write_concatenated_alignment
from string import ascii_uppercase, ascii_lowercase
import numpy as np

### ColabFold
def uni_num(ids):
    ########################################
    pa = {a: 0 for a in ascii_uppercase}
    for a in ["O", "P", "Q"]: pa[a] = 1
    ma = [[{} for k in range(6)], [{} for k in range(6)]]
    for n, t in enumerate(range(10)):
        for i in [0, 1]:
            for j in [0, 4]: ma[i][j][str(t)] = n

    for n, t in enumerate(list(ascii_uppercase) + list(range(10))):
        for i in [0, 1]:
            for j in [1, 2]: ma[i][j][str(t)] = n
        ma[1][3][str(t)] = n
    for n, t in enumerate(ascii_uppercase):
        ma[0][3][str(t)] = n
        for i in [0, 1]: ma[i][5][str(t)] = n
    ########################################
    ipdb.set_trace()
    nums = []
    for uni in ids:
        p = pa[uni[0]]
        tot, num = 1, 0
        if len(uni) == 10:
            for n, u in enumerate(reversed(uni[-4:])):
                num += ma[p][n][u] * tot
                tot *= len(ma[p][n].keys())
        for n, u in enumerate(reversed(uni[:6])):
            num += ma[p][n][u] * tot
            tot *= len(ma[p][n].keys())
        nums.append(num)
    return nums


### RoseTTAFold
def uni2idx(ids):
    '''convert uniprot ids into integers according to the structure
    of uniprot accession numbers'''
    # ipdb.set_trace()
    ids2 = [i + 'AAA0' if len(i) == 6 else i for i in ids]
    arr = np.array([list(s) for s in ids2], dtype='|S1').view(np.uint8)

    for i in [1, 5, 9]:
        arr[:, i] -= ord('0')

    arr[arr >= ord('A')] -= ord('A')
    arr[arr >= ord('0')] -= ord('0') - 26

    arr[:, 0][arr[:, 0] > ord('Q') - ord('A')] -= 3

    arr = arr.astype(np.int64)

    coef = np.array([23, 10, 26, 36, 36, 10, 26, 36, 36, 1], dtype=np.int64)
    coef = np.tile(coef[None, :], [len(ids), 1])

    c1 = [i for i, id_ in enumerate(ids) if id_[0] in 'OPQ' and len(id_) == 6]
    c2 = [i for i, id_ in enumerate(ids) if id_[0] not in 'OPQ' and len(id_) == 6]

    coef[c1] = np.array([3, 10, 36, 36, 36, 1, 1, 1, 1, 1])
    coef[c2] = np.array([23, 10, 26, 36, 36, 1, 1, 1, 1, 1])

    for i in range(1, 10):
        coef[:, -i - 1] *= coef[:, -i]

    return np.sum(arr * coef, axis=-1)


def read_a3m(a3m_file, only_id=False, max_gap_fraction=0.9):

    seqs = OrderedDict()
    contents = open(a3m_file, "r").readlines()
    current_id = None
    current_seq = ""

    regex = re.compile("^[A-Za-z0-9]+$")

    for line in contents:
        if line.startswith(">"):
            if current_id is not None:
                add = True
                if re.fullmatch(regex, current_id.split('_')[-1]) is None:
                    add = False
                if current_id.split('_')[-1].startswith("UPI"):
                    add = False
                if current_seq.count('-') / float(len(current_seq)) > max_gap_fraction:
                    add = False
                if add:
                    seqs[current_id] = current_seq

            current_id = line[1:].rstrip('\n')
            if only_id:
                current_id = current_id.split()[0]
        else:
            current_seq = line.rstrip('\n')

    if current_id is not None:
        add = True
        if re.fullmatch(regex, current_id.split('_')[-1]) is None:
            add = False
        if current_id.split('_')[-1].startswith("UPI"):
            add = False
        if current_seq.count('-') / float(len(current_seq)) > max_gap_fraction:
            add = False
        if add:
            seqs[current_id] = current_seq

    return seqs


def get_interaction_by_accession_num(alignment_1, alignment_2, chain1, chain2):
    # load the monomer alignments
    ali_1 = read_a3m(alignment_1, True, 0.9)
    ali_2 = read_a3m(alignment_2, True, 0.9)

    ids_1 = [id1.split('_')[1] for id1 in ali_1 if id1.startswith('Uni')]
    ids_2 = [id2.split('_')[1] for id2 in ali_2 if id2.startswith('Uni')]

    if len(ids_1) == 0 or len(ids_2) == 0:
        return pd.DataFrame(
        {"id_1": [chain1], "id_2": [chain2]}
    )

    hash1 = uni2idx(ids_1)
    hash2 = uni2idx(ids_2)

    idx1, idx2 = np.where(np.abs(hash1[:, None] - hash2[None, :]) < 10)

    return pd.DataFrame(
        {"id_1": ["UniRef100_" + ids_1[i] for i in idx1], "id_2": ["UniRef100_" + ids_2[j] for j in idx2]}
    )


def generate_alignments_for_heterodimers(inparams):

    chain1, chain2, monomers_msa_dir, outdir = inparams

    chain1_a3m = monomers_msa_dir + '/' + chain1 + '.a3m'
    chain2_a3m = monomers_msa_dir + '/' + chain2 + '.a3m'

    if not os.path.exists(chain1_a3m):
        #die(f"Cannot find a3m file for {chain1}: {chain1_a3m}\n")
        return

    if not os.path.exists(chain2_a3m):
        #die(f"Cannot find a3m file for {chain2}: {chain2_a3m}\n")
        return

    print(f"Processing {chain1}_{chain2}")
    outdir = f"{outdir}/{chain1}_{chain2}"
    makedir_if_not_exists(outdir)

    makedir_if_not_exists(outdir + '/' + chain1)
    clean_a3m(chain1_a3m, outdir + '/' + chain1 + '/' + chain1 + '.a3m')

    makedir_if_not_exists(outdir + '/' + chain2)
    clean_a3m(chain2_a3m, outdir + '/' + chain2 + '/' + chain2 + '.a3m')

    chain1_a3m = outdir + '/' + chain1 + '/' + chain1 + '.a3m'
    chain2_a3m = outdir + '/' + chain2 + '/' + chain2 + '.a3m'

    dimers_interaction = get_interaction_by_accession_num(chain1_a3m, chain2_a3m, chain1, chain2)

    # write concatenated alignment with distance filtering
    # TODO: save monomer alignments?
    target_header, target_seq_idx, sequences_full, sequences_monomer_1, sequences_monomer_2 = \
        write_concatenated_alignment(dimers_interaction, chain1_a3m, chain2_a3m, chain1, chain2)

    # save the alignment files
    mon_alignment_file_1 = f"{outdir}/{chain1}/{chain1}_monomer_1.a3m"
    with open(mon_alignment_file_1, "w") as of:
        write_a3m(sequences_monomer_1, of)

    mon_alignment_file_2 = f"{outdir}/{chain2}/{chain2}_monomer_2.a3m"
    with open(mon_alignment_file_2, "w") as of:
        write_a3m(sequences_monomer_2, of)

    dimers_interaction.to_csv(f"{outdir}/{chain1}_{chain2}_interact.csv", index=False)
    print(dimers_interaction)
    
    complex_ailgnment_file = f"{outdir}/{chain1}_{chain2}.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    args = parser.parse_args()

    t1 = time.time()

    params = read_option_file(args.option_file)

    abs_opt = os.path.abspath(args.option_file)

    # test(params)

    check_dirs(params, ['prosys_dir', 'complex_set_work_dir'])

    check_contents(params, ['monomers_cm_database_name', 'monomers_in_complex_database_name',
                            'monomers_not_in_complex_database_name'])

    original_work_dir = params['complex_set_work_dir'] + '/original'

    current_work_dir = params['complex_set_work_dir']

    # generate msas for heterodimers

    os.system(f"cp {original_work_dir}/heterodimers.list {current_work_dir}")

    heterodimers_list = f"{current_work_dir}/heterodimers.list"

    contents = open(heterodimers_list, "r").readlines()

    makedir_if_not_exists(current_work_dir + '/fr_dimers/accession_num')

    for content in contents:

        content = content.replace('\n', '').replace('\r', '')

        pdbcode, chain1, chain2 = content.split()

        generate_alignments_for_heterodimers([chain1, chain2,
                                             current_work_dir + '/fr',
                                             current_work_dir + '/fr_dimers/accession_num'])

    t2 = time.time()

    print("Fr library updating is done. Total time:" + (t2 - t1).__str__())
