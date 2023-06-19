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

        server_scores = []
        elo_scores = []
        alphafold_scores = []

        for line in open(args.indir + '/' + resultfile):
            line = line.rstrip('\n')
            contents = line.split()

            scores = [float(contents[i]) for i in range(1, len(contents))]
            
            if contents[0] == "elo":
                elo_scores += scores
            elif contents[0] == "alphafold":
                alphafold_scores += scores
            else:
                server_scores += scores

        if len(elo_scores) == 0:
            continue

        # print(server_scores)
        if args.mode == "max":
            print_contents += [str(np.max(np.array(server_scores)))]
            print_contents += [str(np.max(np.array(elo_scores)))]
            if len(alphafold_scores) == 0:
                print_contents += ["0"]
            else:
                print_contents += [str(np.max(np.array(alphafold_scores)))]
        else:
            print_contents += [str(np.mean(np.array(server_scores)))]
            print_contents += [str(np.mean(np.array(elo_scores)))]
            if len(alphafold_scores) == 0:
                print_contents += ["0"]
            else:
                print_contents += [str(np.mean(np.array(alphafold_scores)))]

        print(','.join(print_contents))

    
