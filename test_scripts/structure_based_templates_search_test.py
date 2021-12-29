import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.complex_templates_search.structure_based_pipeline import *


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--pdb1', type=is_file, required=True)
    parser.add_argument('--pdb2', type=is_file, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    makedir_if_not_exists(args.output_dir)

    print("Start to test structure based template search pipeline")

    pipeline = Complex_structure_based_template_search_pipeline(params)

    pipeline.search([os.path.abspath(args.pdb1), os.path.abspath(args.pdb2)], args.output_dir)