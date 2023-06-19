import argparse, os
import sys
import logging
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-i", "--indir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-m", "--mode", help="pdb name in lower case", type=str, required=True)

    args = parser.parse_args()

    for resultfile in sorted(os.listdir(args.indir)):
        if resultfile.find('T') < 0 and resultfile.find('H') < 0:
            continue
        
        print_contents = [resultfile.rstrip('.results')]

        iter_scores = []
        default_scores = []

        for line in open(args.indir + '/' + resultfile):
            line = line.rstrip('\n')
            contents = line.split()

            scores = [float(contents[i]) for i in range(1, len(contents))]
            
            if contents[0] == "default_multimer":
                default_scores += scores
            elif contents[0].find("iter") == 0:
                iter_scores += scores

        if len(iter_scores) == 0:
            continue

        # print(server_scores)
        if args.mode == "max":
            print_contents += [str(np.max(np.array(default_scores)))]
            print_contents += [str(np.max(np.array(iter_scores)))]
        else:
            print_contents += [str(np.mean(np.array(default_scores)))]
            print_contents += [str(np.mean(np.array(iter_scores)))]

        print(','.join(print_contents))

    
