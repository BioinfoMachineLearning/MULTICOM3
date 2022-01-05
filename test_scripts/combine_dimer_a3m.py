import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from bml_casp15.complex_alignment_generation.combine_alignment import *
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('--dimerlist', type=str, required=True)
    parser.add_argument('-a', '--aln_dir', type=is_dir, required=True)
    parser.add_argument('-c', '--concatenate_dir', type=is_dir, required=True)
    parser.add_argument('-o', '--output_dir', type=str, required=True)
    parser.add_argument('-f', '--hhfilter', type=str, required=True)

    args = parser.parse_args()

    combine_methods = ['uniclust_oxmatch_a3m',
                       # 'pdb_interact_uniref_a3m',
                       # 'species_interact_uniref_a3m',
                       # 'uniprot_distance_uniref_a3m',
                       # 'string_interact_uniref_a3m',
                       # 'geno_dist_uniref_a3m',
                       # 'pdb_interact_uniref_sto',
                       # 'species_interact_uniref_sto',
                       # 'uniprot_distance_uniref_sto',
                       # 'string_interact_uniref_sto',
                       # 'geno_dist_uniref_sto',
                       # 'pdb_interact_uniprot_sto',
                       # 'species_interact_uniprot_sto',
                       'uniprot_distance_uniprot_sto']

    pipeline = Combine_alignment(hhfilter=args.hhfilter)

    for dimer in os.listdir(args.concatenate_dir):

        print(f"Processing {dimer}")

        chain1, chain2 = dimer.split('_')

        concatenate_inputs = []
        for combine_method in combine_methods:
            concatenate_input = {}
            if combine_method == 'uniclust_oxmatch_a3m':
                concatenate_input['chain1'] = f"{args.aln_dir}/{chain1}/{chain1}_uniclust30.a3m"
                concatenate_input['chain2'] = f"{args.aln_dir}/{chain2}/{chain2}_uniclust30.a3m"
            elif combine_method.find('_uniref_a3m') > 0:
                concatenate_input['chain1'] = f"{args.aln_dir}/{chain1}/{chain1}_uniref30.a3m"
                concatenate_input['chain2'] = f"{args.aln_dir}/{chain2}/{chain2}_uniref30.a3m"
            elif combine_method.find('_uniref_sto') > 0:
                concatenate_input['chain1'] = f"{args.aln_dir}/{chain1}/{chain1}_uniref90.sto"
                concatenate_input['chain2'] = f"{args.aln_dir}/{chain2}/{chain2}_uniref90.sto"
            elif combine_method.find('_uniprot_sto') > 0:
                concatenate_input['chain1'] = f"{args.aln_dir}/{chain1}/{chain1}_uniprot.sto"
                concatenate_input['chain2'] = f"{args.aln_dir}/{chain2}/{chain2}_uniprot.sto"
            concatenate_input[
                'interaction_csv'] = f"{args.concatenate_dir}/{chain1}_{chain2}/{combine_method}/{chain1}_{chain2}_interact.csv"
            concatenate_input[
                'concatenate_a3m'] = f"{args.concatenate_dir}/{chain1}_{chain2}/{combine_method}/{chain1}_{chain2}.a3m"
            concatenate_inputs += [concatenate_input]

        outdir = f"{args.output_dir}/{chain1}_{chain2}"
        makedir_if_not_exists(outdir)
        pipeline.combine_alignment(chain1, chain2, concatenate_inputs, outdir)



