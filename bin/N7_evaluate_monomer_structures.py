import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.monomer_structure_evaluation.pipeline import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('targetname', None, 'option file')
flags.DEFINE_string('fasta_file', None, 'option file')
flags.DEFINE_string('input_monomer_dir', None, 'Monomer model directory')
flags.DEFINE_string('input_multimer_dir', None, 'Monomer model directory')
flags.DEFINE_string('chain_in_multimer', None, 'option file')
flags.DEFINE_string('output_dir', None, 'Monomer model directory')
flags.DEFINE_boolean('use_gpu', True, 'Monomer model directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    output_dir = FLAGS.output_dir + '/' + FLAGS.targetname

    makedir_if_not_exists(output_dir)

    pipeline = Monomer_structure_evaluation_pipeline(params=params, run_methods=["alphafold", "enQA"], use_gpu=FLAGS.use_gpu)

    pipeline.process(FLAGS.targetname, FLAGS.fasta_file, FLAGS.input_monomer_dir, output_dir, FLAGS.input_multimer_dir,
                     FLAGS.chain_in_multimer)

    print("The template search dimers has finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'targetname',
        'fasta_file',
        'input_monomer_dir',
        'output_dir',
    ])
    app.run(main)
