import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import check_file, check_dir, makedir_if_not_exists, check_contents, read_option_file
from bml_casp15.quaternary_structure_generation.iterative_search_pipeline_v0 import *
from absl import flags
from absl import app
import pathlib

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_list('fasta_paths', None, 'Monomer alignment directory')
flags.DEFINE_string('inpdb_dir', None, 'Monomer alignment directory')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('atomdir', None, 'Output directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    pipeline = Multimer_iterative_generation_pipeline_monomer(params)

    inpdb_dir = os.path.abspath(FLAGS.inpdb_dir)

    cwd = os.getcwd()

    all_multimer_res_all = {'targetname': [], 'tmscore': [], 'tmalign': []}
    all_multimer_res_max = {'targetname': [], 'tmscore': [], 'tmalign': []}

    for fasta_file in FLAGS.fasta_paths:
        os.chdir(cwd)
        fasta_file = os.path.abspath(fasta_file)
        output_dir = os.path.abspath(FLAGS.output_dir)
        targetname = pathlib.Path(fasta_file).stem

        if not os.path.exists(inpdb_dir):
            continue

        multimer_res_all, multimer_res_max = pipeline.search(fasta_file, inpdb_dir,
                                                             output_dir + '/' + targetname,
                                                             FLAGS.atomdir)

        all_multimer_res_all['targetname'] += multimer_res_all['targetname']
        all_multimer_res_all['tmscore'] += multimer_res_all['tmscore']
        all_multimer_res_all['tmalign'] += multimer_res_all['tmalign']

        all_multimer_res_max['targetname'] += multimer_res_max['targetname']
        all_multimer_res_max['tmscore'] += multimer_res_max['tmscore']
        all_multimer_res_max['tmalign'] += multimer_res_max['tmalign']

    cwd = os.getcwd()

    df = pd.DataFrame(all_multimer_res_all)
    df.to_csv(FLAGS.output_dir + '/all_multimer_res_all.csv')

    df = pd.DataFrame(all_multimer_res_max)
    df.to_csv(FLAGS.output_dir + '/all_multimer_res_max.csv')


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'inpdb_dir',
        'output_dir'
    ])
    app.run(main)
