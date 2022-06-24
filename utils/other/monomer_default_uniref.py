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

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(FLAGS.fasta_path)

    targetname = pathlib.Path(FLAGS.fasta_path).stem

    sequence = None
    for line in open(FLAGS.fasta_path):
        line = line.rstrip('\n').strip()
        if line.startswith('>'):
            targetname_in_fasta = line[1:].split()[0]
            if targetname_in_fasta != targetname:
                print("Warning: fasta file name doesn't match with fasta content!")
        else:
            sequence = line

    outdir = FLAGS.output_dir  # + '/' + targetname

    makedir_if_not_exists(outdir)

    N1_outdir = outdir + '/N1_monomer_alignments_generation'
    makedir_if_not_exists(N1_outdir)

    print("#################################################################################################")
    print(f"1. Start to generate alignments for monomers")

    result = run_monomer_msa_pipeline(fasta=FLAGS.fasta_path, outdir=N1_outdir, params=params, only_monomer=True)

    if result is None:
        raise RuntimeError('The monomer alignment generation has failed!')

    N1_outdir_img = outdir + '/N1_monomer_alignments_generation_img'
    makedir_if_not_exists(N1_outdir_img)
    img_msa = run_monomer_msa_pipeline_img(params=params, fasta=FLAGS.fasta_path, outdir=N1_outdir_img)

    print("#################################################################################################")
    print("2. Start to generate template for monomer")

    N2_outdir = outdir + '/N2_monomer_template_search'

    makedir_if_not_exists(N2_outdir)

    template_file = run_monomer_template_search_pipeline(params=params, targetname=targetname, sequence=sequence,
                                                         a3m=f"{N1_outdir}/{targetname}_uniref90.sto", outdir=N2_outdir)

    if template_file is None:
        raise RuntimeError("Program failed in step 2: monomer template search")

    print("The generation for monomer template has finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("3. Start to generate tertiary structure for monomers using alphafold")
    N3_outdir = outdir + '/N3_monomer_structure_generation'
    makedir_if_not_exists(N3_outdir)
    if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                        run_methods=['default_uniref_22'],
                                                        fasta_path=FLAGS.fasta_path,
                                                        alndir=N1_outdir, templatedir=N2_outdir, outdir=N3_outdir):
        print("Program failed in step 3: monomer structure generation")

    print("The prediction for monomers has finished!")

    print("#################################################################################################")

    print("4. Start to evaluate monomer models")

    N4_outdir = outdir + '/N4_monomer_structure_evaluation'

    makedir_if_not_exists(N4_outdir)

    result = run_monomer_evaluation_pipeline(params=params, targetname=targetname, fasta_file=FLAGS.fasta_path,
                                             input_monomer_dir=N3_outdir, outputdir=N4_outdir,
                                             generate_egnn_models=True)

    if result is None:
        raise RuntimeError("Program failed in step 4: monomer model evaluation")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
