import argparse, os
import sys
import logging
import numpy as np
import json
import math

PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
calc_neff = '/home/bml_casp15/BML_CASP15/tools/calNf'

def read_fasta(infasta):
    targets = []
    seqs = []
    for line in open(infasta):
        line = line.rstrip('\n')
        if line[0] == ">":
            targets += [line[1:]]
        elif len(line) > 0:
            seqs += [line]
    return targets, seqs

def get_alignment_depth_v2(fasta_file, msadir):

    subunits, seqs = read_fasta(fasta_file)
    total_length = (int)(np.sum(np.array([len(seq) for seq in seqs])))

    #calculate the depth of each subunits
    all_monomer_depth = 0
    for i in range(len(seqs)):
        chain_id = PDB_CHAIN_IDS[i]
        chain_seqs = []
        multimer_a3m = msadir + '/' + chain_id + '/multimer_final.a3m'
        if os.path.exists(multimer_a3m):
            a3m_headers, a3m_seqs = read_fasta(multimer_a3m)
            chain_seqs += a3m_seqs

        monomer_a3m = msadir + '/' + chain_id + '/monomer_final.a3m'
        a3m_headers, a3m_seqs = read_fasta(monomer_a3m)
        chain_seqs += a3m_seqs
        with open(f"{chain_id}.aln", 'w') as fw:
            fw.write('\n'.join(chain_seqs))
            fw.write('\n')    

        neff_file = chain_id + '.neff'
        if not os.path.exists(neff_file):
            cmd = f"{calc_neff} {chain_id}.aln 0.8 2 > {neff_file}"
            os.system(cmd)
        Neff = float(open(neff_file).readlines()[0].rstrip('\n'))
        all_monomer_depth += len(seqs[i])/total_length * Neff

    return all_monomer_depth

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-i", "--indir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-f", "--fastadir", help="pdb name in lower case", type=str, required=True)

    args = parser.parse_args()

    for target in sorted(os.listdir(args.indir)):
        os.chdir(args.indir + '/' + target)
        alignment_depth = get_alignment_depth_v2(args.fastadir + '/' + target + '.fasta', args.indir + '/' + target) 
        print(f"{target}\t{alignment_depth}")
