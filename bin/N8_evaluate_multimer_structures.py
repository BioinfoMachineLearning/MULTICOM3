import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.quaternary_structure_evaluation.pipeline import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('input_dir', None, 'Monomer model directory')
flags.DEFINE_string('output_dir', None, 'Monomer model directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    pipeline = Quaternary_structure_evaluation_pipeline(params=params)

    pipeline.process(FLAGS.input_dir, FLAGS.output_dir)

    print("The template search dimers has finished!")

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'input_dir',
        'output_dir'
    ])
    app.run(main)
