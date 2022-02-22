import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import check_file, check_dir, makedir_if_not_exists, check_contents, read_option_file
from bml_casp15.tertiary_structure_generation.iterative_search_pipeline import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Monomer alignment directory')
flags.DEFINE_string('inpdb_dir', None, 'Monomer alignment directory')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS

def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    check_dir(FLAGS.inpdb_dir)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    pipeline = Monomer_iterative_generation_pipeline(params)

    pipeline.search(FLAGS.fasta_path, FLAGS.inpdb_dir, FLAGS.output_dir)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'inpdb_dir',
        'output_dir'
    ])
    app.run(main)
