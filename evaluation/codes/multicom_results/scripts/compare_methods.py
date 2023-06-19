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

    methods = ['default_af','default_comp','default_img','default_mul_newest','default_multimer','default_pdb','default_pdb70','default_struct','default_uniclust30','default_uniref30_22','img_multimer','lewis','pdb_iter_uniref_a3m','spec_af','spec_colab_iter','spec_comp','spec_iter_uniprot_sto','spec_iter_uniref_a3m','spec_iter_uniref_sto','spec_pdb','spec_pdb70','spec_struct','str_af','str_comp','str_iter_uniprot_sto','str_iter_uniref_a3m','str_iter_uniref_sto','str_pdb','str_pdb70','str_struct','uniclust_oxmatch_a3m','unidist_uniprot_img','unidist_uniprot_sto','unidist_uniref_a3m','unidist_uniref_sto']
    targets = []
    found_methods = []
    for target in sorted(os.listdir(args.indir)):
        print(target)
        if  target.find('T') == 0 or target.find('H') == 0:
            targets += [target.rstrip('.results')]
    
    targets = sorted(targets)
    print(targets)
    
    for method in methods:
        if os.path.exists(args.indir + '/' + method + '.results'):
            found_methods += [method]

    methods = sorted(found_methods)
    print(methods)

    for method in methods:
        data_dict = {'target': [], 'max': [], 'avg': [], 'def_max': [], 'def_avg': []}
        for line in open(args.indir + '/' + method + '.results'):
            line = line.rstrip('\n')
            contents = line.split()
            if len(contents) == 1:
                continue
            else:
                data_dict['target'] += [contents[0]]
                data_dict['max'] += [np.max(np.array(contents[1:],dtype=float))]
                data_dict['avg'] += [np.mean(np.array(contents[1:],dtype=float))]

                content = os.popen(f"grep {contents[0]} {args.indir}/default_multimer.results").read()
                contents = content.split()
                data_dict['def_max'] += [np.max(np.array(contents[1:],dtype=float))]
                data_dict['def_avg'] += [np.mean(np.array(contents[1:],dtype=float))]
        pd.DataFrame(data_dict).to_csv(args.outdir + '/' + method + '.csv')



