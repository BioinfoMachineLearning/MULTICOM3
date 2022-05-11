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
    run_quaternary_structure_generation_homo_pipeline_v2, \
    run_quaternary_structure_generation_pipeline_foldseek, run_multimer_refinement_pipeline, \
    run_multimer_evaluation_pipeline, run_monomer_msa_pipeline_img, foldseek_iterative_monomer_input, \
    copy_same_sequence_msas

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
    N5_outdir = FLAGS.output_dir + '/N5_monomer_structure_refinement'
    N6_outdir = FLAGS.output_dir + '/N6_monomer_structure_evaluation'
    N7_outdir = FLAGS.output_dir + '/N7_monomer_structure_refinement'
    img_msas = {}

    print("#################################################################################################")

    print("#################################################################################################")
    print("1-3. Start to generate monomer models")

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

            N1_monomer_outdir_img = N1_outdir_img + '/' + monomer_id
            makedir_if_not_exists(N1_monomer_outdir_img)
            img_msas[chain_id] = run_monomer_msa_pipeline_img(params=params,
                                                              fasta=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                              outdir=N1_monomer_outdir_img)

            N2_monomer_outdir = N2_outdir + '/' + monomer_id
            makedir_if_not_exists(N2_monomer_outdir)
            template_file = run_monomer_template_search_pipeline(targetname=monomer_id, sequence=monomer_id,
                                                                 a3m=f"{N1_monomer_outdir}/{monomer_id}_uniref90.sto",
                                                                 outdir=N2_monomer_outdir, params=params)
            if template_file is None:
                raise RuntimeError(f"Program failed in step 2: monomer {monomer_id} template search")

            N3_monomer_outdir = N3_outdir + '/' + monomer_id
            makedir_if_not_exists(N3_monomer_outdir)
            if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                                run_methods=['default', 'default+seq_template',
                                                                             'default_uniclust30',
                                                                             'original', 'original+seq_template',
                                                                             'colabfold', 'colabfold+seq_template'],
                                                                fasta_path=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                                alndir=N1_monomer_outdir,
                                                                templatedir=N2_monomer_outdir,
                                                                outdir=N3_monomer_outdir):
                print(f"Program failed in step 3: monomer {monomer_id} structure generation")

            N4_monomer_outdir = N4_outdir + '/' + monomer_id
            makedir_if_not_exists(N4_monomer_outdir)
            result = run_monomer_evaluation_pipeline(params=params,
                                                     targetname=monomer_id,
                                                     fasta_file=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                     input_monomer_dir=N3_outdir + '/' + monomer_id,
                                                     outputdir=N4_monomer_outdir, generate_egnn_models=True)
            if result is None:
                raise RuntimeError(f"Program failed in step 4: monomer {monomer_id} model evaluation")
            monomer_qas_res[monomer_id] = result

            N5_monomer_outdir = N5_outdir + '/' + monomer_id
            makedir_if_not_exists(N5_monomer_outdir)
            final_dir = N5_monomer_outdir + '_final'

            os.system(f"cp {result['pairwise_af_avg']} {N5_monomer_outdir}")
            ref_ranking = pd.read_csv(result['pairwise_af_avg'])  # apollo or average ranking or the three qas
            refine_inputs = []
            for i in range(5):
                pdb_name = ref_ranking.loc[i, 'model']
                refine_input = iterative_refine_pipeline.refinement_input(
                    fasta_path=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                    pdb_path=f"{N4_monomer_outdir}/pdb/{pdb_name}",
                    pkl_path=f"{N4_monomer_outdir}/pkl/{pdb_name.replace('.pdb', '.pkl')}",
                    msa_path=f"{N4_monomer_outdir}/msa/{pdb_name.replace('.pdb', '.a3m')}")
                refine_inputs += [refine_input]

            run_monomer_refinement_pipeline(params=params, refinement_inputs=refine_inputs,
                                            outdir=N5_monomer_outdir, finaldir=final_dir, prefix="refine")

            print("The refinement for the top-ranked monomer models has been finished!")

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

            N1_monomer_outdir_img = N1_outdir_img + '/' + monomer_id
            makedir_if_not_exists(N1_monomer_outdir_img)
            img_msas[chain_id] = run_monomer_msa_pipeline_img(params=params, fasta=FLAGS.fasta_path,
                                                              outdir=N1_monomer_outdir_img)

            N2_monomer_outdir = N2_outdir + '/' + monomer_id
            makedir_if_not_exists(N2_monomer_outdir)

            # make a copy
            N3_monomer_outdir = N3_outdir + '/' + monomer_id
            N4_monomer_outdir = N4_outdir + '/' + monomer_id
            N5_monomer_outdir = N5_outdir + '/' + monomer_id

            if not os.path.exists(N3_monomer_outdir):
                os.system(f"cp -r {N3_outdir}/{processed_seuqences[monomer_sequence]} {N3_monomer_outdir}")
            os.system(f"cp -r {N4_outdir}/{processed_seuqences[monomer_sequence]} {N4_monomer_outdir}")

            if not os.path.exists(N4_monomer_outdir):
                for msa in os.listdir(N4_monomer_outdir + '/msa'):
                    os.system(
                        f"sed -i 's/>{processed_seuqences[monomer_sequence]}/>{monomer_id}/g' {N4_monomer_outdir}/msa/{msa}")

            if not os.path.exists(N5_monomer_outdir):
                os.system(f"cp -r {N5_outdir}/{processed_seuqences[monomer_sequence]} {N5_monomer_outdir}")

    print("#################################################################################################")

    print("#################################################################################################")

    print("#################################################################################################")
    print("8. Start to generate complex alignments")

    N8_outdir = FLAGS.output_dir + '/N8_complex_alignments_concatenation'
    makedir_if_not_exists(N8_outdir)

    try:
        concat_methods = ['pdb_interact', 'species_interact', 'uniclust_oxmatch']
        run_concatenate_dimer_msas_pipeline(
            multimer=','.join([chain_id_map[chain_id].description for chain_id in chain_id_map]),
            run_methods=concat_methods,
            monomer_aln_dir=N1_outdir, outputdir=N8_outdir, params=params, is_homomers=True)

    except Exception as e:
        print(e)
        print("Program failed in step 8")

    print("#################################################################################################")

    print("#################################################################################################")

    print("9. Start to search complex templates based on monomer structures")

    N9_outdir = FLAGS.output_dir + '/N9_complex_templates_search'

    run_complex_template_search_pipeline(multimers=[chain_id_map[chain_id].description for chain_id in chain_id_map],
                                         monomer_aln_dir=N1_outdir,
                                         monomer_model_dir=N3_outdir,
                                         outdir=N9_outdir, params=params)

    print("#################################################################################################")

    print("#################################################################################################")

    print("10. Start to generate complex quaternary structures")
    N10_outdir = FLAGS.output_dir + '/N10_quaternary_structure_generation'
    makedir_if_not_exists(N10_outdir)

    if not run_quaternary_structure_generation_homo_pipeline_v2(params=params,
                                                                fasta_path=FLAGS.fasta_path,
                                                                chain_id_map=chain_id_map,
                                                                aln_dir=N1_outdir,
                                                                complex_aln_dir=N8_outdir,
                                                                template_dir=N9_outdir,
                                                                monomer_model_dir=N3_outdir,
                                                                output_dir=N10_outdir):
        print("Program failed in step 10")

    print("Complex quaternary structure generation has been finished!")

    print("#################################################################################################")

    processed_seuqences = []

    for chain_id in chain_id_map:
        monomer_id = chain_id_map[chain_id].description
        monomer_sequence = chain_id_map[chain_id].sequence

        if monomer_sequence not in processed_seuqences:

            default_alphafold_msa = N3_monomer_outdir + '/default/msas/monomer_final.a3m'

            if len(open(default_alphafold_msa).readlines()) < 500 * 2:
                while not os.path.exists(img_msas[chain_id]):
                    # sleep for 5 mins
                    time.sleep(300)

                print("Found img alignment, start to run monomer model generation again")

                if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                                    run_methods=['img', 'img+seq_template'],
                                                                    fasta_path=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                                    alndir=N1_monomer_outdir_img,
                                                                    templatedir=N2_monomer_outdir,
                                                                    outdir=N3_monomer_outdir):
                    print("Program failed in step 3: monomer structure generation")

                print("6. Start to evaluate monomer models")

                N6_monomer_outdir = N6_outdir + '/' + monomer_id

                makedir_if_not_exists(N6_monomer_outdir)

                result = run_monomer_evaluation_pipeline(params=params,
                                                         targetname=monomer_id,
                                                         fasta_file=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                         input_monomer_dir=N3_monomer_outdir,
                                                         outputdir=N6_monomer_outdir, generate_egnn_models=True)

                if result is None:
                    raise RuntimeError("Program failed in step 6: monomer model evaluation")
                monomer_qas_res[monomer_id] = result

                print("The evaluation for monomer models has been finished!")

                print(
                    "#################################################################################################")

                print(
                    "#################################################################################################")

                print("7. Start to refine monomer models based on the qa rankings")

                N7_monomer_outdir = N7_outdir + '/' + monomer_id

                makedir_if_not_exists(N7_monomer_outdir)

                old_ref_ranking = copy.deepcopy(ref_ranking)
                refined_models = [old_ref_ranking.loc[i, 'model'] for i in range(5)]
                ref_ranking = pd.read_csv(result['pairwise_af_avg'])  # apollo or average ranking or the three qas
                os.system(f"cp {result['pairwise_af_avg']} {N7_monomer_outdir}")

                refine_inputs = []
                for i in range(5):
                    pdb_name = ref_ranking.loc[i, 'model']
                    if pdb_name not in refined_models:
                        refine_input = iterative_refine_pipeline.refinement_input(
                            fasta_path=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                            pdb_path=N6_monomer_outdir + '/pdb/' + pdb_name,
                            pkl_path=N6_monomer_outdir + '/pkl/' + pdb_name.replace(
                                '.pdb', '.pkl'),
                            msa_path=N6_monomer_outdir + '/msa/' + pdb_name.replace(
                                '.pdb', '.a3m'))
                        refine_inputs += [refine_input]
                    else:
                        os.system(f"cp -r {N5_monomer_outdir}/{pdb_name} {N7_monomer_outdir}")

                final_dir = N7_monomer_outdir + '_final'
                run_monomer_refinement_pipeline(params=params, refinement_inputs=refine_inputs,
                                                outdir=N7_monomer_outdir,
                                                finaldir=final_dir)

                print("The refinement for the top-ranked monomer models has been finished!")

            else:
                print("The alphafold monomer alignment has more than 200 sequences, stop waiting for img alignment")
                break

            processed_seuqences[monomer_sequence] = monomer_id

        else:
            # make a copy
            N3_monomer_outdir = N3_outdir + '/' + monomer_id

            if not os.path.exists(N3_monomer_outdir):
                os.system(f"cp -r {N3_outdir}/{processed_seuqences[monomer_sequence]}/img* {N3_monomer_outdir}")

            N6_monomer_outdir = N6_outdir + '/' + monomer_id
            if not os.path.exists(N6_monomer_outdir):
                if os.path.exists(N6_outdir + '/' + processed_seuqences[monomer_sequence]):
                    os.system(f"cp -r {N6_outdir}/{processed_seuqences[monomer_sequence]} {N6_monomer_outdir}")
                    for msa in os.listdir(N6_monomer_outdir + '/msa'):
                        os.system(
                            f"sed -i 's/>{processed_seuqences[monomer_sequence]}/>{monomer_id}/g' {N6_monomer_outdir}/msa/{msa}")

            N7_monomer_outdir = N7_outdir + '/' + monomer_id
            if not os.path.exists(N7_monomer_outdir):
                if os.path.exists(N7_outdir + '/' + processed_seuqences[monomer_sequence]):
                    os.system(f"cp -r {N7_outdir}/{processed_seuqences[monomer_sequence]} {N7_monomer_outdir}")

    print("#################################################################################################")

    print("11. Start to run multimer iterative generation pipeline using top-ranked monomer models")

    qa_result_dir = N4_outdir
    if os.path.exists(N6_outdir):
        qa_result_dir = N6_outdir

    pipeline_inputs = []
    for i in range(2):
        monomer_pdb_dirs = {}
        monomer_alphafold_a3ms = {}
        pdb_name = None
        for chain_id in chain_id_map:
            monomer_id = chain_id_map[chain_id].description
            monomer_ranking = read_qa_txt_as_df(monomer_qas_res[monomer_id]['apollo'])
            pdb_name = monomer_ranking.loc[i, 'model']
            break

        for chain_id in chain_id_map:
            monomer_id = chain_id_map[chain_id].description
            monomer_pdb_dirs[chain_id] = f"{qa_result_dir}/{monomer_id}/pdb/{pdb_name}"
            monomer_alphafold_a3ms[chain_id] = f"{qa_result_dir}/{monomer_id}/msa/{pdb_name.replace('.pdb', '.a3m')}"

        print(monomer_alphafold_a3ms)
        pipeline_inputs += [foldseek_iterative_monomer_input(monomer_pdb_dirs=monomer_pdb_dirs,
                                                             monomer_alphafold_a3ms=monomer_alphafold_a3ms)]

    if os.path.exists(N6_outdir):
        for chain_id in chain_id_map:
            monomer_id = chain_id_map[chain_id].description
            monomer_pdb_dirs[chain_id] = f"{qa_result_dir}/{monomer_id}/pdb/img_1.pdb"
            monomer_alphafold_a3ms[chain_id] = f"{qa_result_dir}/{monomer_id}/msa/img_1.a3m"

        pipeline_inputs += [foldseek_iterative_monomer_input(monomer_pdb_dirs=monomer_pdb_dirs,
                                                             monomer_alphafold_a3ms=monomer_alphafold_a3ms)]

    if not run_quaternary_structure_generation_pipeline_foldseek(params=params, fasta_path=FLAGS.fasta_path,
                                                                 chain_id_map=chain_id_map,
                                                                 pipeline_inputs=pipeline_inputs, outdir=N10_outdir,
                                                                 is_homomers=True):
        print("Program failed in step 6 iterative")

    print("Complex quaternary structure generation has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("8. Start to evaluate multimer models")

    N11_outdir = FLAGS.output_dir + '/N11_multimer_structure_evaluation'
    multimer_qa_result = run_multimer_evaluation_pipeline(fasta_path=FLAGS.fasta_path,
                                                          params=params, monomer_model_dir=qa_result_dir,
                                                          chain_id_map=chain_id_map,
                                                          indir=N10_outdir, outdir=N11_outdir,
                                                          stoichiometry=FLAGS.stoichiometry)

    print("#################################################################################################")

    print("#################################################################################################")

    print("9. Start to refine multimer models based on the qa rankings")

    N12_outdir = FLAGS.output_dir + '/N12_multimer_structure_refinement'

    makedir_if_not_exists(N12_outdir)
    ref_ranking = pd.read_csv(multimer_qa_result['pairwise_af_avg'])  # apollo or average ranking or the three qas

    refine_inputs = []
    for i in range(5):
        pdb_name = ref_ranking.loc[i, 'model']
        msa_paths = {}
        for chain_id in chain_id_map:
            msa_paths[chain_id] = dict(
                paired_msa=f"{N11_outdir}/msa/{chain_id_map[chain_id].description}/{pdb_name.replace('.pdb', '')}.paired.a3m",
                monomer_msa=f"{N11_outdir}/msa/{chain_id_map[chain_id].description}/{pdb_name.replace('.pdb', '')}.monomer.a3m")
        print(msa_paths)
        refine_input = iterative_refine_pipeline_multimer.refinement_input_multimer(chain_id_map=chain_id_map,
                                                                                    fasta_path=FLAGS.fasta_path,
                                                                                    pdb_path=N11_outdir + '/pdb/' + pdb_name,
                                                                                    pkl_path=N11_outdir + '/pkl/' + pdb_name.replace(
                                                                                        '.pdb', '.pkl'),
                                                                                    msa_paths=msa_paths)
        refine_inputs += [refine_input]

    final_dir = N12_outdir + '_final'
    run_multimer_refinement_pipeline(params=params, refinement_inputs=refine_inputs, outdir=N12_outdir,
                                     finaldir=final_dir, stoichiometry="homomer")

    print("The refinement for the top-ranked multimer models has been finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
