import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from multicom3.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom3.monomer_structure_refinement import iterative_refine_pipeline
from multicom3.common.protein import read_qa_txt_as_df, complete_result
from multicom3.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline_v2, run_monomer_evaluation_pipeline, \
    run_monomer_refinement_pipeline, run_monomer_msa_pipeline_img
import pandas as pd
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to monomer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_boolean('run_img', False, 'Whether to use IMG alignment to generate models')
FLAGS = flags.FLAGS


def main(argv):

    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '4.0'

    FLAGS.fasta_path = os.path.abspath(FLAGS.fasta_path)
    FLAGS.output_dir = os.path.abspath(FLAGS.output_dir)

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(FLAGS.fasta_path)

    targetname = pathlib.Path(FLAGS.fasta_path).stem
    sequence = open(FLAGS.fasta_path).readlines()[1].rstrip('\n').strip()

    outdir = FLAGS.output_dir

    makedir_if_not_exists(outdir)

    N1_outdir = os.path.join(outdir, 'N1_monomer_alignments_generation')
    makedir_if_not_exists(N1_outdir)

    print("#################################################################################################")
    print(f"1. Start to generate alignments for monomers")

    result = run_monomer_msa_pipeline(fasta=FLAGS.fasta_path, outdir=N1_outdir, params=params, only_monomer=True)

    if result is None:
        raise RuntimeError('The monomer alignment generation has failed!')

    if FLAGS.run_img:
        N1_outdir_img = os.path.join(outdir, 'N1_monomer_alignments_generation_img')
        makedir_if_not_exists(N1_outdir_img)
        img_msa = run_monomer_msa_pipeline_img(params=params, fasta=FLAGS.fasta_path, outdir=N1_outdir_img)

    print("#################################################################################################")
    print("2. Start to generate template for monomer")

    N2_outdir = os.path.join(outdir, 'N2_monomer_template_search')

    makedir_if_not_exists(N2_outdir)

    template_file = run_monomer_template_search_pipeline(params=params, targetname=targetname, sequence=sequence,
                                                         a3m=os.path.join(N1_outdir, targetname+"_uniref90.sto"), outdir=N2_outdir)

    if template_file is None:
        raise RuntimeError("Program failed in step 2: monomer template search")

    print("The generation for monomer template has finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("3. Start to generate tertiary structure for monomers using alphafold")
    N3_outdir = os.path.join(outdir, 'N3_monomer_structure_generation')
    makedir_if_not_exists(N3_outdir)
    if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                        fasta_path=FLAGS.fasta_path,
                                                        alndir=N1_outdir, templatedir=N2_outdir, outdir=N3_outdir):
        print("Program failed in step 3: monomer structure generation")


    if FLAGS.run_img:
        while not os.path.exists(img_msa):
            print(f"Waiting for img alignment: {img_msa}")
            # sleep for 5 mins
            time.sleep(300)

        print("Found img alignment, start to run monomer model generation again")

        os.system("cp " + os.path.join(N1_outdir, targetname + "_uniref90.sto") + " " + N1_outdir_img)
        if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                            run_methods=['img', 'img+seq_template'],
                                                            fasta_path=FLAGS.fasta_path,
                                                            alndir=N1_outdir_img, templatedir=N2_outdir, outdir=N3_outdir):
            print("Program failed in step 3: monomer structure generation: img")

    print("The prediction for monomers has finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("4. Start to evaluate monomer models")

    N4_outdir = os.path.join(outdir, 'N4_monomer_structure_evaluation')

    makedir_if_not_exists(N4_outdir)

    result = run_monomer_evaluation_pipeline(params=params, targetname=targetname, fasta_file=FLAGS.fasta_path,
                                             input_monomer_dir=N3_outdir, outputdir=N4_outdir,
                                             generate_egnn_models=True)

    if result is None:
        raise RuntimeError("Program failed in step 4: monomer model evaluation")

    print("The evaluation for monomer models has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("5. Start to refine monomer models based on the qa rankings")

    N5_outdir_avg = os.path.join(outdir, 'N5_monomer_structure_refinement_avg')

    makedir_if_not_exists(N5_outdir_avg)

    os.system(f"cp {result['pairwise_af_avg']} {N5_outdir_avg}")

    ref_ranking_avg = pd.read_csv(result['pairwise_af_avg'])  # apollo or average ranking or the three qas

    refine_inputs = []
    for i in range(5):
        pdb_name = ref_ranking_avg.loc[i, 'model']
        refine_input = iterative_refine_pipeline.refinement_input(fasta_path=FLAGS.fasta_path,
                                                                  pdb_path=os.path.join(N4_outdir, 'pdb', pdb_name),
                                                                  pkl_path=os.path.join(N4_outdir, 'pkl', pdb_name.replace('.pdb', '.pkl')),
                                                                  msa_path=os.path.join(N4_outdir, 'msa', pdb_name.replace('.pdb', '.a3m')))
        refine_inputs += [refine_input]

    final_dir = N5_outdir_avg + '_final'
    run_monomer_refinement_pipeline(params=params, refinement_inputs=refine_inputs, outdir=N5_outdir_avg,
                                    finaldir=final_dir, prefix="refine")

    print("The refinement for the top-ranked monomer models has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
