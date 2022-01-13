import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import check_file, check_dir, makedir_if_not_exists, check_contents, read_option_file
from bml_casp15.tertiary_structure_generation.pipeline import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('aln_dir', None, 'Monomer alignment directory')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def complete_result(outputdir):
    complete = True
    for i in range(0, 5):
        model = f'{outputdir}/ranked_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)
    check_dir(FLAGS.aln_dir)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    print("Start to generate monomers' structure using alphafold")

    fasta_paths = []

    for monomer in os.listdir(FLAGS.aln_dir):
        if os.path.exists(f"{FLAGS.aln_dir}/{monomer}/{monomer}.fasta"):
            if not complete_result(FLAGS.output_dir + '/' + monomer):
                fasta_paths += [f"{FLAGS.aln_dir}/{monomer}/{monomer}.fasta"]

    print(f"Total {len(fasta_paths)} monomers are generating structures")

    if len(fasta_paths) > 0:
        pipeline = Monomer_tertiary_structure_prediction_pipeline(params)
        pipeline.process(fasta_paths, FLAGS.aln_dir, FLAGS.output_dir)

    print("The prediction for monomers has finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'aln_dir',
        'output_dir'
    ])
    app.run(main)
