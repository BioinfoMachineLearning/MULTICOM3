import argparse, os
import sys
import logging


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-i", "--infile", help="pdb name in lower case", type=str, required=True)

    args = parser.parse_args()

    max_score = 0
    max_model = ""
    for line in open(args.infile):
        contents = line.rstrip('\n').split()
        for i in range(1,len(contents)):
            if float(contents[i]) > max_score:
                max_score = float(contents[i])
                max_model = contents[0] + '_' + str(i)
    
    print(max_score)
    print(max_model)


