import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.quaternary_structure_generation.pipeline import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_list('fasta_paths', None, 'option file')
flags.DEFINE_string('aln_dir', None, 'Monomer model directory')
flags.DEFINE_string('complex_aln_dir', None, 'Monomer model directory')
flags.DEFINE_string('structure_template_dir', None, 'Monomer model directory')
flags.DEFINE_string('sequence_template_dir', None, 'Monomer model directory')
flags.DEFINE_string('monomer_model_dir', None, 'Monomer model directory')
flags.DEFINE_string('output_dir', None, 'Monomer model directory')
FLAGS = flags.FLAGS


def run_pipeline(inparams):
    params, fasta_path, aln_dir, complex_aln_dir, structure_template_dir, sequence_template_dir, monomer_model_dir, output_dir = inparams

    pipeline = Quaternary_structure_prediction_pipeline(params)

    result = None
    try:
        result = pipeline.process(fasta_path=fasta_path,
                                  aln_dir=aln_dir,
                                  complex_aln_dir=complex_aln_dir,
                                  structure_template_dir=structure_template_dir,
                                  sequence_template_dir=sequence_template_dir,
                                  monomer_model_dir=monomer_model_dir,
                                  output_dir=output_dir)
    except Exception as e:
        print(e)
        return result
    return result


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    process_list = []

    for fasta_path in FLAGS.fasta_paths:

        monomers = []
        for line in open(fasta_path):
            line = line.rstrip('\n').strip()
            if line.startswith('>'):
                monomers += [line[1:]]

        output_dir = FLAGS.output_dir + '/' + '_'.join(monomers)
        makedir_if_not_exists(output_dir)

        complex_aln_dir = FLAGS.complex_aln_dir + '/' + '_'.join(monomers)
        if not os.path.exists(complex_aln_dir):
            raise Exception(f"Cannot find complex alignment directory: {complex_aln_dir}")

        structure_template_dir = FLAGS.structure_template_dir + '/' + '_'.join(monomers)
        if not os.path.exists(structure_template_dir):
            print(f"Cannot find structure template directory: {structure_template_dir}")

        sequence_template_dir = FLAGS.sequence_template_dir + '/' + '_'.join(monomers)
        if not os.path.exists(sequence_template_dir):
            print(f"Cannot find sequence template directory: {sequence_template_dir}")

        inparams = [params, fasta_path, FLAGS.aln_dir, complex_aln_dir, \
                    structure_template_dir, sequence_template_dir, \
                    FLAGS.monomer_model_dir, output_dir]

        # process_list.append([params, [monomer1_template_input, monomer2_template_input], dimer_outdir])

        run_pipeline(inparams)

    # print(f"Total {len(process_list)} dimers to be processed")
    #
    # pool = Pool(processes=10)
    # results = pool.map(run_pipeline, process_list)
    # pool.close()
    # pool.join()

    print("The template search dimers has finished!")

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_paths',
        'aln_dir',
        'complex_aln_dir',
        'structure_template_dir',
        'sequence_template_dir',
        'monomer_model_dir',
        'output_dir'
    ])
    app.run(main)
