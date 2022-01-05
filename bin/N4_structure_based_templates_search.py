import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import check_file, check_dir, makedir_if_not_exists, check_contents, read_option_file
from bml_casp15.complex_templates_search.structure_based_pipeline import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('model_dir', None, 'Monomer model directory')
flags.DEFINE_list('dimerlist', None, 'List of dimers: A_B,B_C')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)
    check_dir(FLAGS.model_dir)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    print("Start to test structure based template search pipeline")

    pipeline = Complex_structure_based_template_search_pipeline(params)

    for dimer in FLAGS.dimerlist:

        monomer1, monomer2 = dimer.split('_')

        monomer1_pdb = f"{FLAGS.model_dir}/{monomer1}/ranked_0.pdb"

        if not os.path.exists(monomer1_pdb):
            print(f"Cannot find teritary structure for {monomer1}: {monomer1_pdb}")
            continue

        monomer2_pdb = f"{FLAGS.model_dir}/{monomer2}/ranked_0.pdb"
        if not os.path.exists(monomer2_pdb):
            print(f"Cannot find teritary structure for {monomer2}: {monomer2_pdb}")
            continue

        outdir = FLAGS.output_dir + '/' + dimer

        makedir_if_not_exists(outdir)

        pipeline.search([monomer1_pdb, monomer2_pdb], outdir)

    print("Complex template searching has been finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'model_dir',
        'dimerlist',
        'output_dir'
    ])
    app.run(main)
