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


def read_a3m(infile, max_gap_fraction=0.9):
    '''Read a3m MSA'''
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
               'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12,
               'O': 21, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
               'V': 18, 'W': 19, 'Y': 20, 'U': 21, 'Z': 21, 'X': 21, 'J': 21}

    parsed = []  # Save extracted msa
    species = []
    seqlen = 0
    lc = 0
    with open(infile, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith('>'):  # OX=OrganismIdentifier
                if 'OX=' in line:
                    OX = line.split('OX=')[1]
                    if len(OX) > 0:
                        species.append(int(OX.split(' ')[0]))
                    else:
                        species.append(0)
                else:
                    species.append(0)
                continue
            line = line.rstrip()
            gap_fraction = line.count('-') / float(len(line))
            if gap_fraction <= max_gap_fraction:  # Only use the lines with less than 90 % gaps
                parsed.append([mapping.get(ch, 22) for ch in line if not ch.islower()])
            else:
                if len(species) > 1:
                    species = species[:-1]  # Remove the previously stored species
                    continue
            # Check that the lengths match
            if len(parsed[-1]) != seqlen and lc >= 1:
                parsed = parsed[:-1]
                species = species[:-1]
                continue
            seqlen = len(parsed[-1])
            lc += 1

    return np.array(parsed, dtype=np.int8, order='F'), np.array(species)


def match_top_ox(ox1, ox2, msa1, msa2):
    '''Select the top ox match (first match) in each MSA and merge the sequences to a final MSA file
    in a3m format
    - number of possible combinations
    - median and std in hits per species
    '''
    # Don't remove the zeros (no OX), then the query sequences (first line)
    # will be removed
    matching_ox = np.intersect1d(ox1, ox2)

    ind1 = []  # Index to select from the individual MSAs
    ind2 = []
    ncombos = []
    # Go through all matching and select the first (top) hit
    for ox in matching_ox:
        ind1.append(min(np.argwhere(ox1 == ox)[:, 0]))
        ind2.append(min(np.argwhere(ox2 == ox)[:, 0]))

        ncombos.append(np.argwhere(ox1 == ox).shape[0] * np.argwhere(ox2 == ox).shape[0])

    # Select from MSAs and merge
    merged = np.concatenate((msa1[ind1], msa2[ind2]), axis=1)

    return merged, np.sum(ncombos), np.median(ncombos), np.std(ncombos)


def write_a3m(merged_msa, outfile):
    '''Write a3m MSA'''
    backmap = {1: 'A', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H',
               8: 'I', 9: 'K', 10: 'L', 11: 'M', 12: 'N', 13: 'P', 14: 'Q',
               15: 'R', 16: 'S', 17: 'T', 18: 'V', 19: 'W', 20: 'Y',
               21: '-'}  # Here all unusual AAs and gaps are set to the same char (same in the GaussDCA script)

    with open(outfile, 'w') as file:
        for i in range(len(merged_msa)):
            file.write('>' + str(i) + '\n')
            file.write(''.join([backmap[ch] for ch in merged_msa[i]]) + '\n')

    return None


def generate_alignments_for_heterodimers(inparams):

    chain1, chain2, monomers_msa_dir, outdir = inparams

    chain1_a3m = monomers_msa_dir + '/' + chain1 + '.a3m'
    chain2_a3m = monomers_msa_dir + '/' + chain2 + '.a3m'

    if not os.path.exists(chain1_a3m):
        # die(f"Cannot find a3m file for {chain1}: {chain1_a3m}\n")
        return

    if not os.path.exists(chain2_a3m):
        # die(f"Cannot find a3m file for {chain2}: {chain2_a3m}\n")
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

    # MSA1
    msa1, ox1 = read_a3m(chain1_a3m)
    # MSA2
    msa2, ox2 = read_a3m(chain2_a3m)

    # Match
    merged_msa, ncombos_total, median_combos, std_combos = match_top_ox(ox1, ox2, msa1, msa2)

    # save the alignment files
    mon_alignment_file_1 = f"{outdir}/{chain1}/{chain1}_monomer_1.a3m"
    write_a3m(msa1, mon_alignment_file_1)

    mon_alignment_file_2 = f"{outdir}/{chain2}/{chain2}_monomer_2.a3m"
    write_a3m(msa1, mon_alignment_file_2)

    complex_ailgnment_file = f"{outdir}/{chain1}_{chain2}.a3m"
    write_a3m(merged_msa, complex_ailgnment_file)


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

    makedir_if_not_exists(current_work_dir + '/fr_dimers/ox')

    for content in contents:
        content = content.replace('\n', '').replace('\r', '')

        pdbcode, chain1, chain2 = content.split()

        generate_alignments_for_heterodimers([chain1, chain2,
                                              current_work_dir + '/fr',
                                              current_work_dir + '/fr_dimers/ox'])

    t2 = time.time()

    print("Fr library updating is done. Total time:" + (t2 - t1).__str__())
