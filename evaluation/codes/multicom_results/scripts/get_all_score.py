import argparse, os
import sys
import logging
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-i", "--infile", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="pdb name in lower case", type=str, required=True)

    args = parser.parse_args()

    tmscores = []
    models = []
    for line in open(args.infile):
        contents = line.rstrip('\n').split()
        for i in range(1,len(contents)):
            models += [f"{contents[0]}_{i+1}"]
            tmscores += [float(contents[i])]

    pd.DataFrame({'model': models, 'tmscore': tmscores}).to_csv(args.outfile)

