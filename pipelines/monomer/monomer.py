import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from bml_casp15.common.protein import read_qa_txt_as_df, complete_result
from bml_casp15.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline, run_monomer_evaluation_pipeline, \
    run_monomer_refinement_pipeline, run_monomer_msa_pipeline_img

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

    outdir = FLAGS.output_dir#  + '/' + targetname

    makedir_if_not_exists(outdir)

    N1_outdir = outdir + '/N1_monomer_alignments_generation'
    makedir_if_not_exists(N1_outdir)

    print("#################################################################################################")
    print(f"1. Start to generate alignments for monomers")

    result = run_monomer_msa_pipeline(fasta=FLAGS.fasta_path, outdir=N1_outdir, params=params)

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
    if not run_monomer_structure_generation_pipeline(params=params,
                                                     run_methods=['default', 'default+seq_template',
                                                                  'original', 'original+seq_template',
                                                                  'colabfold', 'colabfold+seq_template'],
                                                     fasta_path=FLAGS.fasta_path,
                                                     alndir=N1_outdir, templatedir=N2_outdir, outdir=N3_outdir):
        print("Program failed in step 3: monomer structure generation")

    print("The prediction for monomers has finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("4. Start to evaluate monomer models")

    N4_outdir = outdir + '/N4_monomer_structure_evaluation'

    makedir_if_not_exists(N4_outdir)

    result = run_monomer_evaluation_pipeline(params=params, targetname=targetname, fasta_file=FLAGS.fasta_path,
                                             input_monomer_dir=N3_outdir, outputdir=N4_outdir, generate_egnn_models=True)

    if result is None:
        raise RuntimeError("Program failed in step 4: monomer model evaluation")

    print("The evaluation for monomer models has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("5. Start to refine monomer models based on the qa rankings")

    N5_outdir = outdir + '/N5_monomer_structure_refinement'

    makedir_if_not_exists(N5_outdir)

    os.system(f"cp {result['apollo']} {N5_outdir}")

    ref_ranking = read_qa_txt_as_df(result['apollo'])  # apollo or average ranking or the three qas

    refine_inputs = []
    for i in range(5):
        pdb_name = ref_ranking.loc[i, 'model']
        refine_input = iterative_refine_pipeline.refinement_input(fasta_path=FLAGS.fasta_path,
                                                                  pdb_path=N4_outdir + '/pdb/' + pdb_name,
                                                                  pkl_path=N4_outdir + '/pkl/' + pdb_name.replace(
                                                                      '.pdb', '.pkl'),
                                                                  msa_path=N4_outdir + '/msa/' + pdb_name.replace(
                                                                      '.pdb', '.a3m'))
        refine_inputs += [refine_input]

    final_dir = N5_outdir + '_final'
    run_monomer_refinement_pipeline(params=params, refinement_inputs=refine_inputs,
                                    outdir=N5_outdir, finaldir=final_dir, ranking_df=ref_ranking)

    print("The refinement for the top-ranked monomer models has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("Check the alignment depth")

    default_alphafold_msa = N3_outdir + '/default/msas/monomer_final.a3m'

    if len(open(default_alphafold_msa).readlines()) < 200 * 2:
        while not os.path.exists(img_msa):
            # sleep for 5 mins
            time.sleep(300)

        print("Found img alignment, start to run monomer model generation again")

        if not run_monomer_structure_generation_pipeline(params=params,
                                                         run_methods=['img', 'img+seq_template'],
                                                         fasta_path=FLAGS.fasta_path,
                                                         alndir=N1_outdir_img, templatedir=N2_outdir, outdir=N3_outdir):
            print("Program failed in step 3: monomer structure generation")

        print("6. Start to evaluate monomer models")

        N6_outdir = outdir + '/N6_monomer_structure_evaluation'

        makedir_if_not_exists(N6_outdir)

        result = run_monomer_evaluation_pipeline(targetname=targetname, fasta_file=FLAGS.fasta_path,
                                                 inputdir=N3_outdir, outputdir=N6_outdir, generate_egnn_models=True)

        if result is None:
            raise RuntimeError("Program failed in step 6: monomer model evaluation")

        print("The evaluation for monomer models has been finished!")

        print("#################################################################################################")

        print("#################################################################################################")

        print("7. Start to refine monomer models based on the qa rankings")

        N7_outdir = outdir + '/N7_monomer_structure_refinement'

        makedir_if_not_exists(N7_outdir)

        old_ref_ranking = copy.deepcopy(ref_ranking)
        refined_models = [old_ref_ranking.loc[i, 'model'] for i in range(5)]
        ref_ranking = read_qa_txt_as_df(result['apollo'])  # apollo or average ranking or the three qas
        os.system(f"cp {result['apollo']} {N7_outdir}")

        refine_inputs = []
        for i in range(5):
            pdb_name = ref_ranking.loc[i, 'model']
            if pdb_name not in refined_models:
                refine_input = iterative_refine_pipeline.refinement_input(fasta_path=FLAGS.fasta_path,
                                                                          pdb_path=N6_outdir + '/pdb/' + pdb_name,
                                                                          pkl_path=N6_outdir + '/pkl/' + pdb_name.replace(
                                                                              '.pdb', '.pkl'),
                                                                          msa_path=N6_outdir + '/msa/' + pdb_name.replace(
                                                                              '.pdb', '.a3m'))
                refine_inputs += [refine_input]
            else:
                os.system(f"cp -r {N5_outdir}/{pdb_name} {N7_outdir}")

        final_dir = N7_outdir + '/final'
        run_monomer_refinement_pipeline(params=params, refinement_inputs=refine_inputs, outdir=N7_outdir, finaldir=final_dir)

        print("The refinement for the top-ranked monomer models has been finished!")

    else:
        print("The alphafold monomer alignment has more than 200 sequences, stop waiting for img alignment")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
