import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer
from bml_casp15.common.util import is_dir, is_file, read_option_file


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--fasta', type=is_file, required=True)
    parser.add_argument('--hhblits', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_dir, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    hhblits_runner = hhblits.HHBlits(
        binary_path=args.hhblits,
        databases=[params['uniref_db_dir'] + '/' + params['uniref_db']])

    hhblits_result = hhblits_runner.query(args.fasta, args.outdir)

    print(hhblits_result)

