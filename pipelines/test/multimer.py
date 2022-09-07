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
    run_quaternary_structure_generation_pipeline_v2, \
    run_quaternary_structure_generation_pipeline_foldseek, run_multimer_refinement_pipeline, \
    run_multimer_evaluation_pipeline, run_monomer_msa_pipeline_img, foldseek_iterative_monomer_input, \
    copy_same_sequence_msas

from absl import flags
from absl import app
import copy
import pandas as pd
import pathlib

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fastadir', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    makedir_if_not_exists(FLAGS.output_dir)

    for fasta_file in os.listdir(FLAGS.fastadir):
        outputdir = FLAGS.output_dir + pathlib.Path(FLAGS.fastadir + '/' + fasta_file).stem
        makedir_if_not_exists(outputdir)

        N1_outdir = outputdir + '/N1_monomer_alignments_generation'

        with open(FLAGS.fastadir + '/' + fasta_file) as f:
            input_fasta_str = f.read()
        input_seqs, input_descs = parse_fasta(input_fasta_str)
        chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                           descriptions=input_descs)

        processed_seuqences = {}
        for chain_id in chain_id_map:
            monomer_id = chain_id_map[chain_id].description
            monomer_sequence = chain_id_map[chain_id].sequence

            if monomer_sequence not in processed_seuqences:

                with open(f"{outputdir}/{monomer_id}.fasta", "w") as fw:
                    write_fasta({monomer_id: monomer_sequence}, fw)
                N1_monomer_outdir = N1_outdir + '/' + monomer_id
                makedir_if_not_exists(N1_monomer_outdir)
                result = run_monomer_msa_pipeline(f"{outputdir}/{monomer_id}.fasta", N1_monomer_outdir, params)
                if result is None:
                    raise RuntimeError(f"Program failed in step 1: monomer {monomer_id} alignment generation")

                processed_seuqences[monomer_sequence] = monomer_id
            else:
                with open(f"{outputdir}/{monomer_id}.fasta", "w") as fw:
                    write_fasta({monomer_id: monomer_sequence}, fw)
                N1_monomer_outdir = N1_outdir + '/' + monomer_id
                makedir_if_not_exists(N1_monomer_outdir)

                copy_same_sequence_msas(srcdir=f"{N1_outdir}/{processed_seuqences[monomer_sequence]}",
                                        trgdir=N1_monomer_outdir,
                                        srcname=processed_seuqences[monomer_sequence],
                                        trgname=monomer_id)

        N4_outdir = outputdir + '/N4_complex_alignments_concatenation'
        makedir_if_not_exists(N4_outdir)

        try:
            concat_methods = ['pdb_interact', 'species_interact', 'uniclust_oxmatch',
                              'string_interact', 'uniprot_distance', 'species_colabfold_interact']
            run_concatenate_dimer_msas_pipeline(
                multimer=','.join([chain_id_map[chain_id].description for chain_id in chain_id_map]),
                run_methods=concat_methods,
                monomer_aln_dir=N1_outdir, outputdir=N4_outdir, params=params)
        except Exception as e:
            print(e)
            print("Program failed in step 5")

        print("6. Start to generate complex quaternary structures")
        N6_outdir = outputdir + '/N6_quaternary_structure_generation'
        makedir_if_not_exists(N6_outdir)

        run_methods_part1 = ['default_uniclust30',
                             'uniclust_oxmatch_a3m',
                             'pdb_interact_uniref_a3m',
                             'species_interact_uniref_a3m',
                             'uniprot_distance_uniref_a3m',
                             'string_interact_uniref_a3m',
                             # 'geno_dist_uniref_a3m',
                             'pdb_interact_uniref_sto',
                             'species_interact_uniref_sto',
                             'uniprot_distance_uniref_sto',
                             'string_interact_uniref_sto',
                             # 'geno_dist_uniref_sto',
                             'pdb_interact_uniprot_sto',
                             'species_interact_uniprot_sto',
                             'uniprot_distance_uniprot_sto',
                             'string_interact_uniprot_sto']

        if not run_quaternary_structure_generation_pipeline_v2(params=params,
                                                               fasta_path=FLAGS.fastadir + '/' + fasta_file,
                                                               chain_id_map=chain_id_map,
                                                               aln_dir=N1_outdir,
                                                               complex_aln_dir=N4_outdir,
                                                               template_dir="",
                                                               monomer_model_dir="",
                                                               output_dir=N6_outdir,
                                                               run_methods=run_methods_part1):
            print("Program failed in step 7")

        print("Complex quaternary structure generation has been finished!")

        print("#################################################################################################")

        print("#################################################################################################")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fastadir',
        'output_dir'
    ])
    app.run(main)