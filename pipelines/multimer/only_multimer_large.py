import os, sys, argparse, time
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import write_fasta
from bml_casp15.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from bml_casp15.quaternary_structure_refinement import iterative_refine_pipeline_multimer
from bml_casp15.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline, run_monomer_evaluation_pipeline, run_monomer_refinement_pipeline, \
    run_concatenate_dimer_msas_pipeline, run_complex_template_search_pipeline, \
    run_quaternary_structure_generation_pipeline, \
    run_quaternary_structure_generation_pipeline_foldseek, run_multimer_refinement_pipeline, \
    run_multimer_evaluation_pipeline, run_monomer_msa_pipeline_img, foldseek_iterative_monomer_input, \
    copy_same_sequence_msas, run_quaternary_structure_generation_pipeline_default

from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('stoichiometry', None, 'stoichiometry')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(FLAGS.fasta_path)

    makedir_if_not_exists(FLAGS.output_dir)

    N1_outdir = FLAGS.output_dir + '/N1_monomer_alignments_generation'
    N1_outdir_img = FLAGS.output_dir + '/N1_monomer_alignments_generation_img'
    N2_outdir = FLAGS.output_dir + '/N2_monomer_template_search'
    N3_outdir = FLAGS.output_dir + '/N3_monomer_structure_generation'
    N4_outdir = FLAGS.output_dir + '/N4_monomer_structure_evaluation'
    img_msas = {}

    print("#################################################################################################")

    print("#################################################################################################")
    print("1-4. Start to generate monomer models")

    makedir_if_not_exists(N1_outdir)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    processed_seuqences = {}
    monomer_qas_res = {}
    for chain_id in chain_id_map:
        monomer_id = chain_id_map[chain_id].description
        monomer_sequence = chain_id_map[chain_id].sequence

        if monomer_sequence not in processed_seuqences:
            with open(f"{FLAGS.output_dir}/{monomer_id}.fasta", "w") as fw:
                write_fasta({monomer_id: monomer_sequence}, fw)
            N1_monomer_outdir = N1_outdir + '/' + monomer_id
            makedir_if_not_exists(N1_monomer_outdir)
            result = run_monomer_msa_pipeline(f"{FLAGS.output_dir}/{monomer_id}.fasta", N1_monomer_outdir, params)
            if result is None:
                raise RuntimeError(f"Program failed in step 1: monomer {monomer_id} alignment generation")

            N2_monomer_outdir = N2_outdir + '/' + monomer_id
            makedir_if_not_exists(N2_monomer_outdir)
            template_file = run_monomer_template_search_pipeline(targetname=monomer_id, sequence=monomer_id,
                                                                 a3m=f"{N1_monomer_outdir}/{monomer_id}_uniref90.sto",
                                                                 outdir=N2_monomer_outdir, params=params)
            if template_file is None:
                raise RuntimeError(f"Program failed in step 2: monomer {monomer_id} template search")

            processed_seuqences[monomer_sequence] = monomer_id

        else:
            with open(f"{FLAGS.output_dir}/{monomer_id}.fasta", "w") as fw:
                write_fasta({monomer_id: monomer_sequence}, fw)
            N1_monomer_outdir = N1_outdir + '/' + monomer_id
            makedir_if_not_exists(N1_monomer_outdir)

            copy_same_sequence_msas(srcdir=f"{N1_outdir}/{processed_seuqences[monomer_sequence]}",
                                    trgdir=N1_monomer_outdir,
                                    srcname=processed_seuqences[monomer_sequence],
                                    trgname=monomer_id)

    print("#################################################################################################")

    print("#################################################################################################")

    print("2. Start to generate complex quaternary structures")
    N2_outdir = FLAGS.output_dir + '/N2_quaternary_structure_generation'
    makedir_if_not_exists(N2_outdir)

    if not run_quaternary_structure_generation_pipeline_default(params=params,
                                                                fasta_path=FLAGS.fasta_path,
                                                                chain_id_map=chain_id_map,
                                                                aln_dir=N1_outdir,
                                                                output_dir=N2_outdir):
        print("Program failed in step 7")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
