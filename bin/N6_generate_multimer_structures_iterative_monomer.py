import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import check_file, check_dir, makedir_if_not_exists, check_contents, read_option_file
from bml_casp15.quaternary_structure_generation.iterative_search_pipeline_monomer import *
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

    pipeline = Multimer_iterative_generation_pipeline(params)

    all_multimer_res_all = {'targetname': [], 'model': [], 'start_lddt': [], 'end_lddt': [], 'start_tmscore': [],
                            'end_tmscore': []}
    all_multimer_res_avg = {'targetname': [], 'start_lddt': [], 'end_lddt': [], 'start_tmscore': [], 'end_tmscore': []}
    all_multimer_res_max = {'targetname': [], 'start_lddt': [], 'end_lddt': [], 'start_tmscore': [], 'end_tmscore': []}

    inpdb_dir = os.path.abspath(FLAGS.inpdb_dir)

    for fasta_path in FLAGS.fasta_paths:
        fasta_path = os.path.abspath(fasta_path)
        output_dir = os.path.abspath(FLAGS.output_dir)
        targetname = pathlib.Path(fasta_path).stem

        if not os.path.exists(inpdb_dir):
            continue

        multimer_res_all, multimer_res_avg, multimer_res_max = pipeline.search(fasta_path, inpdb_dir,
                                                                               output_dir + '/' + targetname,
                                                                               FLAGS.atomdir)

        all_multimer_res_all['targetname'] += multimer_res_all['targetname']
        all_multimer_res_all['model'] += multimer_res_all['model']
        all_multimer_res_all['start_lddt'] += multimer_res_all['start_lddt']
        all_multimer_res_all['end_lddt'] += multimer_res_all['end_lddt']
        all_multimer_res_all['start_tmscore'] += multimer_res_all['start_tmscore']
        all_multimer_res_all['end_tmscore'] += multimer_res_all['end_tmscore']

        all_multimer_res_avg['targetname'] += multimer_res_avg['targetname']
        all_multimer_res_avg['start_lddt'] += multimer_res_avg['start_lddt']
        all_multimer_res_avg['end_lddt'] += multimer_res_avg['end_lddt']
        all_multimer_res_avg['start_tmscore'] += multimer_res_avg['start_tmscore']
        all_multimer_res_avg['end_tmscore'] += multimer_res_avg['end_tmscore']

        all_multimer_res_max['targetname'] += multimer_res_max['targetname']
        all_multimer_res_max['start_lddt'] += multimer_res_max['start_lddt']
        all_multimer_res_max['end_lddt'] += multimer_res_max['end_lddt']
        all_multimer_res_max['start_tmscore'] += multimer_res_max['start_tmscore']
        all_multimer_res_max['end_tmscore'] += multimer_res_max['end_tmscore']

    df = pd.DataFrame(all_multimer_res_all)
    df.to_csv(os.path.abspath(FLAGS.output_dir) + '/all_multimer_res_all.csv')

    df = pd.DataFrame(all_multimer_res_avg)
    df.to_csv(os.path.abspath(FLAGS.output_dir) + '/all_multimer_res_avg.csv')

    df = pd.DataFrame(all_multimer_res_max)
    df.to_csv(os.path.abspath(FLAGS.output_dir) + '/all_multimer_res_max.csv')


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_paths',
        'inpdb_dir',
        'output_dir'
    ])
    app.run(main)
