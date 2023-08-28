import argparse, os
import sys
import logging
import numpy as np
from scipy.stats import skew
from scipy.stats import kurtosis

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-i", "--indir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-o", "--outdir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-c", "--chains", help="pdb name in lower case", type=str, required=True)


    args = parser.parse_args()
    for inpdb in os.listdir(args.indir):
        chains = set()
        for line in open(args.indir + '/' + inpdb):
            if line.startswith('ATOM'):
                chain_id = line[21]
                chains.add(chain_id)
        print(chains)
        if len(chains) == args.chains:
            os.system(f"cp {args.indir}/{inpdb} {args.outdir}")
    incomplete_pdbs = len(os.listdir(args.indir)) - len(os.listdir(args.outdir))
    print(incomplete_pdbs)