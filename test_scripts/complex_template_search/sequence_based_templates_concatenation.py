import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.complex_templates_search.sequence_based_pipeline import *


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--dimerlist', type=str, required=True)
    parser.add_argument('--aln_dir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    makedir_if_not_exists(args.output_dir)

    pipeline = Complex_sequence_based_template_search_pipeline(params)

    for dimer in args.dimerlist.split(','):

        chain1, chain2 = dimer.split('_')

        chain1_template_a3m = f"{args.aln_dir}/{chain1}/{chain1}_uniref90.sto"

        chain2_template_a3m = f"{args.aln_dir}/{chain2}/{chain2}_uniref90.sto"

        if not os.path.exists(chain1_template_a3m):
            continue

        if not os.path.exists(chain2_template_a3m):
            continue

        monomer1_template_input = monomer_template_input(name=chain1, msa_path=chain1_template_a3m)

        monomer2_template_input = monomer_template_input(name=chain2, msa_path=chain2_template_a3m)

        pipeline.search([monomer1_template_input, monomer2_template_input], args.output_dir)