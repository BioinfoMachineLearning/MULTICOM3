import os, sys, argparse, time
from multiprocessing import Pool
from multicom3.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom3.monomer_alignment_generation.alignment import write_fasta
from multicom3.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from multicom3.multimer_structure_refinement import iterative_refine_pipeline_multimer
from multicom3.monomer_structure_refinement import iterative_refine_pipeline
from multicom3.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline_v2, run_monomer_evaluation_pipeline, run_monomer_refinement_pipeline, \
    run_monomer_msas_concatenation_pipeline, run_monomer_templates_concatenation_pipeline, \
    run_multimer_structure_generation_pipeline_v2, \
    run_multimer_structure_generation_pipeline_foldseek, run_multimer_refinement_pipeline, \
    run_multimer_evaluation_pipeline, run_monomer_msa_pipeline_img, foldseek_iterative_monomer_input, \
    copy_same_sequence_msas

from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fasta')
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

    makedir_if_not_exists(FLAGS.output_dir)

    N1_outdir = FLAGS.output_dir + '/N1_monomer_alignments_generation'
    N1_outdir_img = FLAGS.output_dir + '/N1_monomer_alignments_generation_img'
    N2_outdir = FLAGS.output_dir + '/N2_monomer_template_search'
    N3_outdir = FLAGS.output_dir + '/N3_monomer_structure_generation'
    img_msas = {}

    print("#################################################################################################")

    print("#################################################################################################")
    print("1-3. Start to generate monomer models")

    makedir_if_not_exists(N1_outdir)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(fasta_string=input_fasta_str, sort_sequence=True)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    processed_seuqences = {}
    for chain_id in chain_id_map:
        monomer_id = chain_id
        monomer_sequence = chain_id_map[chain_id].sequence

        if monomer_sequence not in processed_seuqences:

            with open(f"{FLAGS.output_dir}/{monomer_id}.fasta", "w") as fw:
                write_fasta({chain_id: monomer_sequence}, fw)
            N1_monomer_outdir = N1_outdir + '/' + monomer_id
            makedir_if_not_exists(N1_monomer_outdir)
            result = run_monomer_msa_pipeline(f"{FLAGS.output_dir}/{monomer_id}.fasta", N1_monomer_outdir, params)
            if result is None:
                raise RuntimeError(f"Program failed in step 1: monomer {monomer_id} alignment generation")

            if FLAGS.run_img:
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
                                                                fasta_path=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                                alndir=N1_monomer_outdir,
                                                                templatedir=N2_monomer_outdir,
                                                                outdir=N3_monomer_outdir):
                print(f"Program failed in step 3: monomer {monomer_id} structure generation")

            processed_seuqences[monomer_sequence] = monomer_id
        else:
            with open(f"{FLAGS.output_dir}/{monomer_id}.fasta", "w") as fw:
                write_fasta({chain_id: monomer_sequence}, fw)
            N1_monomer_outdir = N1_outdir + '/' + monomer_id
            makedir_if_not_exists(N1_monomer_outdir)

            copy_same_sequence_msas(srcdir=f"{N1_outdir}/{processed_seuqences[monomer_sequence]}",
                                    trgdir=N1_monomer_outdir,
                                    srcname=processed_seuqences[monomer_sequence],
                                    trgname=monomer_id)

            N1_monomer_outdir_img = N1_outdir_img + '/' + monomer_id
            makedir_if_not_exists(N1_monomer_outdir_img)

            N2_monomer_outdir = N2_outdir + '/' + monomer_id
            makedir_if_not_exists(N2_monomer_outdir)

            # make a copy
            N3_monomer_outdir = N3_outdir + '/' + monomer_id

            if not os.path.exists(N3_monomer_outdir):
                os.system(f"cp -r {N3_outdir}/{processed_seuqences[monomer_sequence]} {N3_monomer_outdir}")

    print("#################################################################################################")

    print("#################################################################################################")
    print("4. Start to generate complex alignments")

    N4_outdir = FLAGS.output_dir + '/N4_monomer_alignments_concatenation'
    makedir_if_not_exists(N4_outdir)

    try:
        concat_methods = ['pdb_interact', 'species_interact', 'uniclust_oxmatch',
                          'string_interact', 'uniprot_distance']
        run_monomer_msas_concatenation_pipeline(
            multimer=','.join([chain_id for chain_id in chain_id_map]),
            run_methods=concat_methods,
            monomer_aln_dir=N1_outdir, outputdir=N4_outdir, params=params)
    except Exception as e:
        print(e)
        print("Program failed in step 5")

    print("#################################################################################################")

    print("#################################################################################################")

    print("5. Start to search complex templates based on monomer structures")

    N5_outdir = FLAGS.output_dir + '/N5_monomer_templates_concatenation'

    run_monomer_templates_concatenation_pipeline(multimers=[chain_id for chain_id in chain_id_map],
                                         monomer_aln_dir=N1_outdir,
                                         monomer_model_dir=N3_outdir,
                                         outdir=N5_outdir, params=params)

    print("#################################################################################################")

    print("#################################################################################################")

    print("6. Start to generate complex multimer structures")
    N6_outdir = FLAGS.output_dir + '/N6_multimer_structure_generation'
    makedir_if_not_exists(N6_outdir)

    run_methods = ['default',
                    'default+sequence_based_template_pdb70',
                    'default+structure_based_template',
                    'default+sequence_based_template_pdb',
                    'default+sequence_based_template_complex_pdb',
                    'default+alphafold_model_templates',
                    'uniclust_oxmatch_a3m',
                    'pdb_interact_uniref_a3m',
                    'species_interact_uniref_a3m',
                    'species_interact_uniref_a3m+sequence_based_template_pdb70',
                    'species_interact_uniref_a3m+structure_based_template',
                    'species_interact_uniref_a3m+sequence_based_template_pdb',
                    'species_interact_uniref_a3m+sequence_based_template_complex_pdb',
                    'species_interact_uniref_a3m+alphafold_model_templates',
                    'uniprot_distance_uniref_a3m',
                    'string_interact_uniref_a3m',
                    'pdb_interact_uniref_sto',
                    'species_interact_uniref_sto',
                    'uniprot_distance_uniref_sto',
                    'string_interact_uniref_sto',
                    'string_interact_uniref_sto+sequence_based_template_pdb70',
                    'string_interact_uniref_sto+structure_based_template',
                    'string_interact_uniref_sto+sequence_based_template_pdb',
                    'string_interact_uniref_sto+sequence_based_template_complex_pdb',
                    'string_interact_uniref_sto+alphafold_model_templates',
                    'pdb_interact_uniprot_sto',
                    'species_interact_uniprot_sto',
                    'uniprot_distance_uniprot_sto',
                    'string_interact_uniprot_sto']

    if not run_multimer_structure_generation_pipeline_v2(params=params,
                                                           fasta_path=FLAGS.fasta_path,
                                                           chain_id_map=chain_id_map,
                                                           aln_dir=N1_outdir,
                                                           complex_aln_dir=N4_outdir,
                                                           template_dir=N5_outdir,
                                                           monomer_model_dir=N3_outdir,
                                                           output_dir=N6_outdir,
                                                           run_methods=run_methods):
        print("Program failed in step 7")

    print("Multimer structure generation has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    if FLAGS.run_img:
        img_wait_list = [chain_id for chain_id in chain_id_map]

        img_processed_list = []
        processed_seuqences = {}
        while len(img_processed_list) != len(img_wait_list):
            for chain_id in img_wait_list:
                if chain_id in img_processed_list:
                    continue
                monomer_sequence = chain_id_map[chain_id].sequence
                monomer_id = chain_id
                if monomer_sequence not in processed_seuqences:
                    if os.path.exists(img_msas[chain_id]):
                        print("Found img alignment, start to run monomer model generation again")

                        N1_monomer_outdir_img = N1_outdir_img + '/' + monomer_id
                        os.system(f"cp {N1_outdir}/{monomer_id}/{monomer_id}_uniref90.sto {N1_monomer_outdir_img}")

                        if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                                            run_methods=['img', 'img+seq_template'],
                                                                            fasta_path=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                                            alndir=N1_monomer_outdir_img,
                                                                            templatedir=N2_outdir + '/' + monomer_id,
                                                                            outdir=N3_outdir + '/' + monomer_id):
                            print("Program failed in step 3: monomer structure generation")
                        processed_seuqences[monomer_sequence] = monomer_id
                        img_processed_list += [chain_id]
                else:
                    copy_same_sequence_msas(srcdir=f"{N1_outdir}_img/{processed_seuqences[monomer_sequence]}",
                                            trgdir=f"{N1_outdir}_img/{monomer_id}",
                                            srcname=processed_seuqences[monomer_sequence],
                                            trgname=monomer_id)

                    # make a copy
                    N3_monomer_outdir = N3_outdir + '/' + monomer_id

                    if os.path.exists(N3_monomer_outdir):
                        os.system(f"rm -rf {N3_monomer_outdir}")

                    makedir_if_not_exists(N3_monomer_outdir)
                    os.system(f"cp -r {N3_outdir}/{processed_seuqences[monomer_sequence]} {N3_monomer_outdir}")

            # sleep for 5 mins
            time.sleep(300)

    print("7. Start to evaluate monomer models")

    N7_outdir = FLAGS.output_dir + '/N7_monomer_structure_evaluation'
    monomer_qas_res = {}
    processed_seuqences = {}
    for chain_id_idx, chain_id in enumerate(chain_id_map):
        monomer_id = chain_id
        monomer_sequence = chain_id_map[chain_id].sequence
        if monomer_sequence not in processed_seuqences:
            N7_monomer_outdir = N7_outdir + '/' + monomer_id
            makedir_if_not_exists(N7_monomer_outdir)
            result = run_monomer_evaluation_pipeline(params=params,
                                                     targetname=monomer_id,
                                                     fasta_file=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                     input_monomer_dir=N3_outdir + '/' + monomer_id,
                                                     input_multimer_dir=N6_outdir,
                                                     outputdir=N7_monomer_outdir, generate_egnn_models=True)
            if result is None:
                raise RuntimeError(f"Program failed in step 7: monomer {monomer_id} model evaluation")
            monomer_qas_res[monomer_id] = result
            processed_seuqences[monomer_sequence] = monomer_id
        else:
            # make a copy
            N7_monomer_outdir = N7_outdir + '/' + monomer_id
            os.system(f"cp -r {N7_outdir}/{processed_seuqences[monomer_sequence]} {N7_monomer_outdir}")
            for msa in os.listdir(N7_monomer_outdir + '/msa'):
                os.system(
                    f"sed -i 's/>{processed_seuqences[monomer_sequence]}/>{monomer_id}/g' {N7_monomer_outdir}/msa/{msa}")

    print("#################################################################################################")

    print("#################################################################################################")

    print("8. Start to run multimer iterative generation pipeline using top-ranked monomer models")

    pipeline_inputs = []
    for i in range(2):
        monomer_pdb_dirs = {}
        monomer_alphafold_a3ms = {}
        for chain_id in chain_id_map:
            monomer_id = chain_id
            monomer_ranking = pd.read_csv(monomer_qas_res[monomer_id]['apollo_monomer'])
            pdb_name = monomer_ranking.loc[i, 'model']
            monomer_pdb_dirs[chain_id] = f"{N7_outdir}/{monomer_id}/pdb/{pdb_name}"
            monomer_alphafold_a3ms[chain_id] = f"{N7_outdir}/{monomer_id}/msa/{pdb_name.replace('.pdb', '.a3m')}"
        pipeline_inputs += [foldseek_iterative_monomer_input(monomer_pdb_dirs=monomer_pdb_dirs,
                                                             monomer_alphafold_a3ms=monomer_alphafold_a3ms)]

    if not run_multimer_structure_generation_pipeline_foldseek(params=params, fasta_path=FLAGS.fasta_path,
                                                                 chain_id_map=chain_id_map,
                                                                 pipeline_inputs=pipeline_inputs, outdir=N6_outdir):
        print("Program failed in step 6 iterative")

    print("Multimer structure generation has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("9. Start to evaluate multimer models")

    N8_outdir = FLAGS.output_dir + '/N8_multimer_structure_evaluation'
    multimer_qa_result = run_multimer_evaluation_pipeline(fasta_path=FLAGS.fasta_path,
                                                          params=params,
                                                          chain_id_map=chain_id_map,
                                                          indir=N6_outdir, outdir=N8_outdir)

    print("#################################################################################################")

    print("#################################################################################################")

    print("10. Start to refine multimer models based on the qa rankings")

    N9_outdir = FLAGS.output_dir + '/N9_multimer_structure_refinement'

    makedir_if_not_exists(N9_outdir)
    ref_ranking = pd.read_csv(multimer_qa_result['pairwise_af_avg'])  # apollo or average ranking or the three qas

    refine_inputs = []
    for i in range(5):
        pdb_name = ref_ranking.loc[i, 'model']
        msa_paths = {}
        for chain_id in chain_id_map:
            msa_paths[chain_id] = dict(
                paired_msa=f"{N8_outdir}/msa/{chain_id}/{pdb_name.replace('.pdb', '')}.paired.a3m",
                monomer_msa=f"{N8_outdir}/msa/{chain_id}/{pdb_name.replace('.pdb', '')}.monomer.a3m")

        refine_input = iterative_refine_pipeline_multimer.refinement_input_multimer(chain_id_map=chain_id_map,
                                                                                    fasta_path=FLAGS.fasta_path,
                                                                                    pdb_path=N8_outdir + '/pdb/' + pdb_name,
                                                                                    pkl_path=N8_outdir + '/pkl/' + pdb_name.replace(
                                                                                        '.pdb', '.pkl'),
                                                                                    msa_paths=msa_paths)
        refine_inputs += [refine_input]

    final_dir = N9_outdir + '_final'
    stoichiometry = "heteromer"

    run_multimer_refinement_pipeline(chain_id_map=chain_id_map,
                                     params=params, refinement_inputs=refine_inputs, outdir=N9_outdir,
                                     finaldir=final_dir, stoichiometry=stoichiometry)

    print("The refinement for the top-ranked multimer models has been finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
