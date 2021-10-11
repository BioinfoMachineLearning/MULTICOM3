import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer
from bml_casp15.common.util import is_dir, is_file


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--input_fasta_path', type=is_dir, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    jackhmmer_runner = jackhmmer.Jackhmmer(
        binary_path='/data/bml_casp15/anaconda3/envs/bml_casp15/bin/jackhmmer',
        database_path=[params['uniref_db_dir'] + '/' + uniref_db + '/' + uniref_db])

    hhblits_runner = hhblits.HHBlits(
        binary_path='/data/bml_casp15/anaconda3/envs/bml_casp15/bin/hhblits',
        databases=[params['uniref_db_dir'] + '/' + uniref_db + '/' + uniref_db])

    jackhmmer_result = jackhmmer_runner.query(input_fasta_path)[0]

    hhblits_result = hhblits_runner.query(input_fasta_path)

