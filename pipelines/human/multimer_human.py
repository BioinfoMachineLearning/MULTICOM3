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
    run_multimer_evaluation_pipeline_human, run_monomer_msa_pipeline_img, foldseek_iterative_monomer_input, \
    copy_same_sequence_msas, run_quaternary_structure_generation_homo_pipeline_img

from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('stoichiometry', None, 'Output directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_file(FLAGS.fasta_path)

    makedir_if_not_exists(FLAGS.output_dir)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    monomer_qa_dir = FLAGS.output_dir + '/N1_monomer_structure_evaluation'
    multimer_model_dir = FLAGS.output_dir + '/multimer_models_ori'
    extract_multimer_model_dir = FLAGS.output_dir + '/multimer_models_extract'

    N1_outdir = FLAGS.output_dir + '/N1_quaternary_structure_evaluation'
    multimer_qa_result = run_multimer_evaluation_pipeline_human(fasta_path=FLAGS.fasta_path,
                                                          params=params, monomer_model_dir=monomer_qa_dir,
                                                          chain_id_map=chain_id_map,
                                                          indir=multimer_model_dir,
                                                          extract_dir=extract_multimer_model_dir,
                                                          outdir=N1_outdir,
                                                          stoichiometry=FLAGS.stoichiometry,
                                                          model_count=5)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir',
        'stoichiometry'
    ])
    app.run(main)
