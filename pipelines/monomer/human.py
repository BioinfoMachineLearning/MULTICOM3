import os, sys, argparse, time, copy
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import read_fasta, write_fasta
from bml_casp15.monomer_alignment_generation.pipeline import *
from bml_casp15.monomer_structure_generation.pipeline import *
from bml_casp15.monomer_templates_search.sequence_based_pipeline_pdb import *
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('model_path', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def run_monomer_msa_pipeline(fasta, outdir, params):
    uniref30 = params['uniref_db_dir'] + '/' + params['uniref_db']
    uniclust30 = params['uniclust_db_dir'] + '/' + params['uniclust_db']
    uniref90_fasta = params['uniref90_fasta']
    uniprot_fasta = params['uniprot_fasta']
    smallbfd = ""  # params['smallbfd_database']
    bfd = params['bfd_database']
    mgnify = params['mgnify_database']
    hhblits_binary = params['hhblits_program']
    jackhmmer_binary = params['jackhmmer_program']
    result = None
    try:
        monomer_msa_pipeline = Monomer_alignment_generation_pipeline(jackhmmer_binary,
                                                                     hhblits_binary,
                                                                     uniref90_fasta,
                                                                     mgnify,
                                                                     smallbfd,
                                                                     bfd,
                                                                     uniref30,
                                                                     uniclust30,
                                                                     uniprot_fasta)

        result = monomer_msa_pipeline.process(fasta, outdir)
    except Exception as e:
        print(e)
    return result


def run_monomer_template_search_pipeline(targetname, sequence, a3m, outdir, params):
    template_file = None
    try:
        pipeline = monomer_sequence_based_template_search_pipeline(params)
        template_file = pipeline.search(targetname=targetname, sequence=sequence, a3m=a3m, outdir=outdir)
    except Exception as e:
        print(e)
    return template_file


def run_monomer_evaluation_pipeline(targetname, fasta_file, inputdir, outputdir):
    makedir_if_not_exists(outputdir)
    result = None
    pipeline = Monomer_structure_evaluation_pipeline(params=params,
                                                     run_methods=["apollo", "alphafold", "enQA"],
                                                     use_gpu=True)
    try:
        result = pipeline.process(targetname, fasta_file, inputdir, outputdir)
    except Exception as e:
        print(e)
    return result


def read_qa_txt(infile):
    models = []
    scores = []
    for line in open(infile):
        line = line.rstrip('\n')
        contents = line.split()
        if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or \
                contents[0] == "END":
            continue
        model, score = line.split()
        models += [model]
        scores += [float(score)]
    df = pd.DataFrame({'model': models, 'score': scores})
    df = df.sort_values(by=['score'], ascending=False)
    df.reset_index(inplace=True)
    return df


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

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(fasta_path)

    targetname = pathlib.Path(fasta_path).stem

    sequence = None
    for line in open(fasta_path):
        line = line.rstrip('\n').strip()
        if line.startswith('>'):
            targetname_in_fasta = line[1:].split()[0]
            if targetname_in_fasta != targetname:
                print("Warning: fasta file name doesn't match with fasta content!")
        else:
            sequence = line

    outdir = FLAGS.output_dir# + '/' + targetname

    makedir_if_not_exists(outdir)

    # N1_outdir = outdir + '/N1_monomer_alignments_generation'
    # makedir_if_not_exists(N1_outdir)
    #
    # print("#################################################################################################")
    # print(f"1. Start to generate alignments for monomers")
    #
    # result = run_monomer_msa_pipeline(fasta=fasta_path, outdir=N1_outdir, params=params)
    #
    # if result is None:
    #     raise RuntimeError('The monomer alignment generation has failed!')
    #
    # N1_outdir_img = outdir + '/N1_monomer_alignments_generation_img'
    # makedir_if_not_exists(N1_outdir_img)
    # monomer_msa_pipeline_img = Monomer_alignment_generation_pipeline_img(binary_path=params['deepmsa2_program'],
    #                                                                      bfd_database_path=params['bfd_database'],
    #                                                                      img_database_path=params['img_database'],
    #                                                                      metaclust_database_path=params[
    #                                                                          'metaclust_database'],
    #                                                                      mgnify_database_path=params[
    #                                                                          'mgnify_database'],
    #                                                                      uniref90_database_path=params[
    #                                                                          'uniref90_database'])
    #
    # img_msa = monomer_msa_pipeline_img.process(fasta, N1_outdir_img)
    #
    # print("#################################################################################################")
    # print("2. Start to generate template for monomer")
    #
    # N2_outdir = outdir + '/N2_monomer_template_search'
    #
    # makedir_if_not_exists(N2_outdir)
    #
    # template_file = run_monomer_template_search_pipeline(targetname=targetname, sequence=sequence,
    #                                                      a3m=a3m, outdir=N2_outdir, params=params)
    #
    # if template_file is None:
    #     raise RuntimeError("Program failed in step 2: monomer template search")
    #
    # print("The generation for monomer template has finished!")
    #
    # print("#################################################################################################")
    #
    # print("#################################################################################################")
    #
    # print("3. Start to generate tertiary structure for monomers using alphafold")
    # N3_outdir = outdir + '/N3_monomer_structure_generation'
    # makedir_if_not_exists(N3_outdir)
    # try:
    #     pipeline = Monomer_tertiary_structure_prediction_pipeline(params)
    #     pipeline.process_single(fasta_path=fasta_path,
    #                             alndir=N1_outdir,
    #                             template_dir=N2_outdir,
    #                             outdir=N3_outdir,
    #                             run_methods=['default', 'original', 'colabfold', 'seq_temp_pdb'])
    # except Exception as e:
    #     print(e)
    #     print("Program failed in step 3: monomer structure generation")
    #
    # print("The prediction for monomers has finished!")
    #
    # print("#################################################################################################")

    N3_outdir = FLAGS.model_dir

    print("#################################################################################################")

    print("4. Start to evaluate monomer models")

    N4_outdir = outdir + '/N4_monomer_structure_evaluation'

    makedir_if_not_exists(N4_outdir)

    result = run_monomer_evaluation_pipeline(targetname=targetname, fasta_file=fasta_path,
                                             inputdir=N3_outdir, outputdir=N4_outdir)

    if result is None:
        raise RuntimeError("Program failed in step 4: monomer model evaluation")

    print("The evaluation for monomer models has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("5. Start to refine monomer models based on the qa rankings")

    N5_outdir = outdir + '/N5_monomer_structure_refinement'

    makedir_if_not_exists(N5_outdir)

    ref_ranking = read_qa_txt(result['apollo'])  # apollo or average ranking or the three qas

    pipeline = iterative_refine_pipeline.Monomer_iterative_refinement_pipeline_server(params=params)
    refine_inputs = []
    for i in range(5):
        pdb_name = ref_ranking.loc[i, 'model']
        refine_input = iterative_refine_pipeline.refinement_input(fasta_path=fasta_path,
                                                                  pdb_path=N4_outdir + '/pdb/' + pdb_name,
                                                                  pkl_path=N4_outdir + '/pkl/' + pdb_name.replace(
                                                                      '.pdb', '.pkl'),
                                                                  msa_path=N4_outdir + '/msa/' + pdb_name.replace(
                                                                      '.pdb', '.a3m'))
        refine_inputs += [refine_input]
    pipeline.search(refinement_inputs=refine_inputs, outdir=N5_outdir)

    final_dir = N5_outdir + '/final'
    makedir_if_not_exists(final_dir)

    pipeline = Monomer_refinement_model_selection()
    pipeline.select_v1(indir=N5_outdir, outdir=final_dir)

    print("The refinement for the top-ranked monomer models has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("Check the alignment depth")

    default_alphafold_msa = N3_outdir + '/default/msas/monomer_final.a3m'

    if len(open(default_alphafold_msa).readlines()) < 200:
        while not os.path.exists(img_msa):
            # sleep for 5 mins
            time.sleep(300)

        print("Found img alignment, start to run monomer model generation again")

        pipeline = Monomer_tertiary_structure_prediction_pipeline(params)
        pipeline.process_single(fasta_path=fasta_path,
                                alndir=N1_outdir_img,
                                outdir=N3_outdir,
                                run_methods=['img'])

        print("6. Start to evaluate monomer models")

        N6_outdir = outdir + '/N6_monomer_structure_evaluation'

        makedir_if_not_exists(N6_outdir)

        result = run_monomer_evaluation_pipeline(targetname=targetname, fasta_file=fasta_path,
                                                 inputdir=N3_outdir, outputdir=N6_outdir)

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
        ref_ranking = read_qa_txt(result['apollo'])  # apollo or average ranking or the three qas

        pipeline = iterative_refine_pipeline.Monomer_iterative_refinement_pipeline_server(params=params)
        refine_inputs = []
        for i in range(5):
            pdb_name = ref_ranking.loc[i, 'model']
            if pdb_name not in refined_models:
                refine_input = iterative_refine_pipeline.refinement_input(fasta_path=fasta_path,
                                                                          pdb_path=N6_outdir + '/pdb/' + pdb_name,
                                                                          pkl_path=N6_outdir + '/pkl/' + pdb_name.replace(
                                                                              '.pdb', '.pkl'),
                                                                          msa_path=N6_outdir + '/msa/' + pdb_name.replace(
                                                                              '.pdb', '.a3m'))
                refine_inputs += [refine_input]
            else:
                os.system(f"cp -r {N5_outdir}/{pdb_name} {N7_outdir}")

        pipeline.search(refinement_inputs=refine_inputs, outdir=N7_outdir)

        print("The refinement for the top-ranked monomer models has been finished!")

        final_dir = N7_outdir + '/final'
        makedir_if_not_exists(final_dir)

        pipeline = Monomer_refinement_model_selection()
        pipeline.select_v1(indir=N7_outdir, outdir=final_dir)

    else:
        print("The alphafold monomer alignment has more than 200 sequences, stop waiting for img alignment")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
