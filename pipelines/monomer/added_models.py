import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from bml_casp15.common.protein import read_qa_txt_as_df, complete_result
from bml_casp15.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline, run_monomer_evaluation_pipeline, rerun_monomer_evaluation_pipeline, \
    run_monomer_refinement_pipeline, run_monomer_msa_pipeline_img
import pandas as pd
from absl import flags
from absl import app


flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('model_path', None, 'Path to monomer models')
flags.DEFINE_list('refine_paths', None, 'Path to monomer models')
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

    outdir = FLAGS.output_dir#  + '/' + targetname

    makedir_if_not_exists(outdir)

    print("1. Start to evaluate monomer models")

    N1_outdir = outdir + '/N1_monomer_structure_evaluation'

    os.system(f"cp -r {FLAGS.model_path} {N1_outdir}")

    makedir_if_not_exists(N1_outdir)

    result = rerun_monomer_evaluation_pipeline(params=params, targetname=targetname, fasta_file=FLAGS.fasta_path,
                                                outputdir=N1_outdir)

    if result is None:
        raise RuntimeError("Program failed in step 4: monomer model evaluation")

    print("The evaluation for monomer models has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("2. Start to refine monomer models based on the qa rankings")

    N2_outdir_avg = outdir + '/N2_monomer_structure_refinement_avg'

    makedir_if_not_exists(N2_outdir_avg)

    os.system(f"cp {result['pairwise_af_avg']} {N2_outdir_avg}")

    ref_ranking_avg = pd.read_csv(result['pairwise_af_avg'])  # apollo or average ranking or the three qas

    refined_models = {}
    for refine_path in FLAGS.refine_paths:
        files = os.listdir(refine_path)
        for file in files:
            if file != "pairwise_af_avg.ranking" and file != "alphafold_ranking.csv":
                refined_models[file] = refine_path + '/' + file

    refine_inputs = []
    for i in range(5):
        pdb_name = ref_ranking_avg.loc[i, 'model']
        if pdb_name.replace('.pdb', '') in refined_models:
            os.system(f"cp -r {refined_models[pdb_name.replace('.pdb', '')]} {N2_outdir_avg}")
        refine_input = iterative_refine_pipeline.refinement_input(fasta_path=FLAGS.fasta_path,
                                                                  pdb_path=N1_outdir + '/pdb/' + pdb_name,
                                                                  pkl_path=N1_outdir + '/pkl/' + pdb_name.replace(
                                                                      '.pdb', '.pkl'),
                                                                  msa_path=N1_outdir + '/msa/' + pdb_name.replace(
                                                                      '.pdb', '.a3m'))
        refine_inputs += [refine_input]

    final_dir = N2_outdir_avg + '_final'
    run_monomer_refinement_pipeline(params=params, refinement_inputs=refine_inputs, outdir=N2_outdir_avg,
                                    finaldir=final_dir, prefix="refine")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
