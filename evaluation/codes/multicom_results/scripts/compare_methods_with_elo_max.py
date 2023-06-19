import argparse, os
import sys
import logging
import numpy as np
import json
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-i", "--indir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-o", "--outdir", help="pdb name in lower case", type=str, required=True)
    args = parser.parse_args()

    # methods = ['default', 'default_seq_temp', 'original', 'ori_seq_temp', 'colabfold', 'colab_seq_temp', 'img', 'img_seq_temp']
    methods = []
    targets = []
    for target in sorted(os.listdir(args.indir)):
        print(target)
        if  target.find('T') == 0 or target.find('H') == 0:
            targets += [target.rstrip('.results')]
        else:
            methods += [target[0:target.find('.results')]]
    
    targets = sorted(targets)
    methods = sorted(methods)

    print(targets)
    print(methods)

    data_dict = {'target': [], 'elo_max_tm': []}
    for line in open(args.indir + '/elo.results'):
        line = line.rstrip('\n')
        contents = line.split()
        if len(contents) == 1:
            continue

        tmscores = []
        for i in range(1, len(contents)):
            tmscores += [float(contents[i])]
        
        data_dict['target'] += [contents[0]]
        data_dict['elo_max_tm'] += [str(np.max(np.array(tmscores,dtype=float)))]

    for method in methods:
        data_dict[f"{method}_max_tm"] = []

    for target in data_dict['target']:
        for method in methods:
            content = os.popen(f"grep {target} {args.indir}/{method}.results").read()
            contents = content.split()
            if len(contents) == 0:
                data_dict[f"{method}_max_tm"] += [" "]
            else:
                tmscores = []
                for i in range(1, len(contents)):
                    tmscores += [float(contents[i])]
                    
                if len(tmscores) == 0:
                    data_dict[f"{method}_max_tm"] += [" "]
                else:
                    data_dict[f'{method}_max_tm'] += [str(np.max(np.array(tmscores,dtype=float)))]

    pd.DataFrame(data_dict).to_csv(args.outdir + '/method_summary.csv')



