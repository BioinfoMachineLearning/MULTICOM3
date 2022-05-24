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
        print(processed_seuqences)
        monomer_id = chain_id_map[chain_id].description
        monomer_sequence = chain_id_map[chain_id].sequence

        with open(f"{FLAGS.output_dir}/{monomer_id}.fasta", "w") as fw:
            write_fasta({monomer_id: monomer_sequence}, fw)
        N1_monomer_outdir = N1_outdir + '/' + monomer_id
        makedir_if_not_exists(N1_monomer_outdir)

        result = run_monomer_msa_pipeline(f"{FLAGS.output_dir}/{monomer_id}.fasta", N1_monomer_outdir, params)
        if result is None:
            raise RuntimeError(f"Program failed in step 1: monomer {monomer_id} alignment generation")

        N1_monomer_outdir_img = N1_outdir_img + '/' + monomer_id
        makedir_if_not_exists(N1_monomer_outdir_img)
        img_msas[chain_id] = run_monomer_msa_pipeline_img(params=params, fasta=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                          outdir=N1_monomer_outdir_img)

        N2_monomer_outdir = N2_outdir + '/' + monomer_id
        makedir_if_not_exists(N2_monomer_outdir)
        template_file = run_monomer_template_search_pipeline(targetname=monomer_id, sequence=monomer_id,
                                                             a3m=f"{N1_monomer_outdir}/{monomer_id}_uniref90.sto",
                                                             outdir=N2_monomer_outdir, params=params)
        if template_file is None:
            raise RuntimeError(f"Program failed in step 2: monomer {monomer_id} template search")

        if monomer_sequence not in processed_seuqences:
            N3_monomer_outdir = N3_outdir + '/' + monomer_id
            makedir_if_not_exists(N3_monomer_outdir)
            if not run_monomer_structure_generation_pipeline(params=params,
                                                             run_methods=['default'],
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
                                                     input_monomer_dir=N3_monomer_outdir,
                                                     outputdir=N4_monomer_outdir)
            if result is None:
                raise RuntimeError(f"Program failed in step 4: monomer {monomer_id} model evaluation")

            monomer_qas_res[chain_id] = result
            processed_seuqences[monomer_sequence] = monomer_id
        else:
            N3_monomer_outdir = N3_outdir + '/' + monomer_id
            N4_monomer_outdir = N4_outdir + '/' + monomer_id

            if not os.path.exists(N3_monomer_outdir):
                os.system(f"cp -r {N3_outdir}/{processed_seuqences[monomer_sequence]} {N3_monomer_outdir}")

            if not os.path.exists(N4_monomer_outdir):
                os.system(f"cp -r {N4_outdir}/{processed_seuqences[monomer_sequence]} {N4_monomer_outdir}")

            result = run_monomer_evaluation_pipeline(params=params,
                                                     targetname=monomer_id,
                                                     fasta_file=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                     input_monomer_dir=N3_monomer_outdir,
                                                     outputdir=N4_monomer_outdir)
            monomer_qas_res[chain_id] = result

            for msa in os.listdir(N4_monomer_outdir + '/msa'):
                os.system(
                    f"sed -i 's/>{processed_seuqences[monomer_sequence]}/>{monomer_id}/g' {N4_monomer_outdir}/msa/{msa}")

    print("#################################################################################################")

    print("#################################################################################################")
    print("5. Start to generate complex alignments")

    N5_outdir = FLAGS.output_dir + '/N5_complex_alignments_concatenation'
    makedir_if_not_exists(N5_outdir)

    try:
        concat_methods = ['pdb_interact', 'species_interact', 'uniclust_oxmatch',
                          'string_interact', 'uniprot_distance', 'species_colabfold_interact']
        run_concatenate_dimer_msas_pipeline(
            multimer=','.join([chain_id_map[chain_id].description for chain_id in chain_id_map]),
            run_methods=concat_methods,
            monomer_aln_dir=N1_outdir, outputdir=N5_outdir, params=params)
    except Exception as e:
        print(e)
        print("Program failed in step 5")

    print("#################################################################################################")

    print("#################################################################################################")

    print("6. Start to search complex templates based on monomer structures")

    N6_outdir = FLAGS.output_dir + '/N6_complex_templates_search'

    run_complex_template_search_pipeline(multimers=[chain_id_map[chain_id].description for chain_id in chain_id_map],
                                         monomer_aln_dir=N1_outdir,
                                         monomer_model_dir=N3_outdir,
                                         outdir=N6_outdir, params=params)

    print("#################################################################################################")

    print("#################################################################################################")

    print("7. Start to generate complex quaternary structures")
    N7_outdir = FLAGS.output_dir + '/N7_quaternary_structure_generation'
    makedir_if_not_exists(N7_outdir)

    if not run_quaternary_structure_generation_pipeline(params=params,
                                                        fasta_path=FLAGS.fasta_path,
                                                        chain_id_map=chain_id_map,
                                                        aln_dir=N1_outdir,
                                                        complex_aln_dir=N5_outdir,
                                                        template_dir=N6_outdir,
                                                        monomer_model_dir=N3_outdir,
                                                        output_dir=N7_outdir):
        print("Program failed in step 7")

    pipeline_inputs = []
    for i in range(2):
        monomer_pdb_dirs = {}
        monomer_alphafold_a3ms = {}
        for chain_id in chain_id_map:
            monomer_id = chain_id_map[chain_id].description
            monomer_ranking = read_qa_txt_as_df(monomer_qas_res[chain_id]['apollo'])
            pdb_name = monomer_ranking.loc[i, 'model']
            monomer_pdb_dirs[chain_id] = f"{N4_outdir}/{monomer_id}/pdb/{pdb_name}"
            monomer_alphafold_a3ms[chain_id] = f"{N4_outdir}/{monomer_id}/msa/{pdb_name.replace('.pdb', '.a3m')}"
        pipeline_inputs += [foldseek_iterative_monomer_input(monomer_pdb_dirs=monomer_pdb_dirs,
                                                             monomer_alphafold_a3ms=monomer_alphafold_a3ms)]

    if not run_quaternary_structure_generation_pipeline_foldseek(params=params, fasta_path=FLAGS.fasta_path,
                                                                 chain_id_map=chain_id_map,
                                                                 pipeline_inputs=pipeline_inputs, outdir=N7_outdir):
        print("Program failed in step 7 iterative")

    print("Complex quaternary structure generation has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("8. Start to evaluate multimer models")

    N8_outdir = FLAGS.output_dir + '/N8_multimer_structure_evaluation'
    multimer_qa_result = run_multimer_evaluation_pipeline(params=params, fasta_path=FLAGS.fasta_path,
                                                          chain_id_map=chain_id_map, monomer_model_dir=N3_outdir,
                                                          indir=N7_outdir,
                                                          outdir=N8_outdir, stoichiometry=FLAGS.stoichiometry)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
