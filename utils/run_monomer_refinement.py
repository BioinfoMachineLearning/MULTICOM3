import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from bml_casp15.common.protein import read_qa_txt_as_df, complete_result
from bml_casp15.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline_v2, run_monomer_evaluation_pipeline, \
    run_monomer_refinement_pipeline, run_monomer_msa_pipeline_img
import pandas as pd
from absl import flags
from absl import app
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('indir', None, 'Path to multimer fastas')
flags.DEFINE_string('outdir', None, 'Output directory')
FLAGS = flags.FLAGS

def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    makedir_if_not_exists(FLAGS.outdir)

    params = read_option_file(FLAGS.option_file)

    refine_inputs = []
    for pdb_name in os.listdir(FLAGS.indir):
        if pdb_name.find('.pdb') > 0:
            refine_input = iterative_refine_pipeline.refinement_input(fasta_path=FLAGS.fasta_path,
                                                                      pdb_path=FLAGS.indir + '/' + pdb_name,
                                                                      pkl_path=FLAGS.indir + '/' + pdb_name.replace(
                                                                          '.pdb', '.pkl'),
                                                                      msa_path=FLAGS.indir + '/' + pdb_name.replace(
                                                                          '.pdb', '.a3m'))
            refine_inputs += [refine_input]

    final_dir = FLAGS.outdir + '_final'
    run_monomer_refinement_pipeline(params=params, refinement_inputs=refine_inputs, outdir=FLAGS.outdir,
                                    finaldir=final_dir, prefix="refine")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'outdir'
    ])
    app.run(main)
