import os, sys, argparse, time
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import write_fasta
from bml_casp15.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from bml_casp15.quaternary_structure_refinement import iterative_refine_pipeline_multimer
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from bml_casp15.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline_v2, run_monomer_evaluation_pipeline, run_monomer_refinement_pipeline, \
    run_concatenate_dimer_msas_pipeline, run_complex_template_search_pipeline, \
    run_quaternary_structure_generation_homo_pipeline, \
    run_quaternary_structure_generation_pipeline_foldseek, run_multimer_refinement_pipeline, \
    run_multimer_evaluation_pipeline, run_monomer_msa_pipeline_img, foldseek_iterative_monomer_input, \
    copy_same_sequence_msas, run_quaternary_structure_generation_homo_pipeline_img

from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('targetname', None, 'Path to multimer fastas')
flags.DEFINE_string('custom_msa', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('alndir', None, 'Output directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_file(FLAGS.fasta_path)

    makedir_if_not_exists(FLAGS.output_dir)

    os.chdir(params['alphafold_program_dir'])
    errormsg = ""
    uniref90_sto = FLAGS.alndir + '/' + FLAGS.targetname + '_uniref90.sto'
    if not os.path.exists(uniref90_sto):
        errormsg = errormsg + f"Cannot find uniref90 alignment for {FLAGS.targetname}: {uniref90_sto}\n"

    if len(errormsg) == 0:
        if not complete_result(FLAGS.output_dir):
            try:
                cmd = f"python {params['alphafold_program']} --fasta_path {FLAGS.fasta_path} " \
                      f"--env_dir {params['alphafold_env_dir']} " \
                      f"--database_dir {params['alphafold_database_dir']} " \
                      f"--custom_msa {FLAGS.custom_msa} " \
                      f"--uniref90_sto {uniref90_sto} " \
                      f"--output_dir {FLAGS.output_dir}"
                print(cmd)
                os.system(cmd)
            except Exception as e:
                print(e)
    else:
        print(errormsg)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)




