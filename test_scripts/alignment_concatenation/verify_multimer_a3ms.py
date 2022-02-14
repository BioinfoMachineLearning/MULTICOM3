import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from bml_casp15.complex_alignment_generation.combine_alignment import *
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists
from bml_casp15.monomer_alignment_generation.alignment import *


def read_seqs_from_a3m(infile):
    seqs = []
    for line in open(infile):
        line = line.rstrip('\n')
        if line.startswith('>'):
            continue
        else:
            seqs += [line]
    return seqs


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--aln_dir', type=is_dir, required=True)
    parser.add_argument('-c', '--concatenate_dir', type=is_dir, required=True)
    args = parser.parse_args()
    for multimer in os.listdir(args.concatenate_dir):
        chains = multimer.split('_')
        for method in os.listdir(args.concatenate_dir + '/' + multimer):
            combined_a3m = f"{args.concatenate_dir}/{multimer}/{method}/{method}.a3m"
            if not os.path.exists(combined_a3m):
                continue
            print(f"Found {method} a3m for {multimer}")
            combined_seqs = read_seqs_from_a3m(combined_a3m)
            monomer_seqs = []
            for i in range(len(chains)):
                chain = chains[i]

                chain_seqs = read_seqs_from_a3m(f"{args.concatenate_dir}/{multimer}/{method}/{chain}_con.a3m")
                monomer_seqs += [chain_seqs]
                if len(combined_seqs) != len(chain_seqs):
                    raise Exception(
                        f"The length of combined sequence doesn't match with the monomer sequences: "
                        f"{len(combined_seqs)} and {len(chain_seqs)}")

                src_a3m = None
                aln_format = None
                if method.find('uniref_a3m') > 0:
                    src_a3m = f"{args.aln_dir}/{chain}/{chain}_uniref30.a3m"
                    aln_format = "a3m"
                elif method.find('uniprot_sto') > 0:
                    src_a3m = f"{args.aln_dir}/{chain}/{chain}_uniprot.sto"
                    aln_format = "stockholm"
                elif method.find('uniref_sto') > 0:
                    src_a3m = f"{args.aln_dir}/{chain}/{chain}_uniref90.sto"
                    aln_format = "stockholm"
                elif method.find('uniclust') >= 0:
                    src_a3m = f"{args.aln_dir}/{chain}/{chain}_uniclust30.a3m"
                    aln_format = "a3m"
                print(f"{method} : {src_a3m}")
                with open(src_a3m) as f:
                    aln = Alignment.from_file(f, format=aln_format)
                    for chain_seq in chain_seqs:
                        if chain_seq not in aln.seqs and chain_seq != aln.main_seq:
                            raise Exception(
                                f"Cannot find the sequence in the original a3m file: {chain_seq} in {src_a3m}")

            for i in range(len(combined_seqs)):
                monomer_comb_seq = "".join([monomer_seqs[j][i] for j in range(len(chains))])
                if monomer_comb_seq != combined_seqs[i]:
                    raise Exception(
                        f"The combined sequence doesn't match with the monomer sequences: {monomer_comb_seq} and {combined_seqs[i]}")
