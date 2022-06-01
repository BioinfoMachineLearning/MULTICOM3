import os, sys, argparse, time
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import write_fasta
from bml_casp15.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from bml_casp15.quaternary_structure_refinement import iterative_refine_pipeline_multimer
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from bml_casp15.common.pipeline import run_monomer_evaluation_pipeline_human

from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('chain', 'A', 'Output directory')
flags.DEFINE_boolean('refine', False, 'Output directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_file(FLAGS.fasta_path)

    makedir_if_not_exists(FLAGS.output_dir)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    monomer_model_dir = FLAGS.output_dir + '/monomer_models'
    multimer_model_dir = FLAGS.output_dir + '/multimer_models'

    N1_outdir = FLAGS.output_dir + '/N1'
    N2_outdir = FLAGS.output_dir + '/N2'
    makedir_if_not_exists(N1_outdir)
    processed_seuqences = {}
    for chain_id in chain_id_map:
        if chain_id != FLAGS.chain:
            continue
        monomer_id = chain_id_map[chain_id].description
        monomer_sequence = chain_id_map[chain_id].sequence
        if monomer_sequence not in processed_seuqences:
            with open(f"{FLAGS.output_dir}/{monomer_id}.fasta", "w") as fw:
                write_fasta({monomer_id: monomer_sequence}, fw)
            N1_monomer_outdir = N1_outdir + '/' + monomer_id
            makedir_if_not_exists(N1_monomer_outdir)
            result = run_monomer_evaluation_pipeline_human(params=params,
                                                           targetname=monomer_id,
                                                           fasta_file=f"{FLAGS.output_dir}/{monomer_id}.fasta",
                                                           input_monomer_dir=monomer_model_dir + '/' + monomer_id,
                                                           input_multimer_dir=multimer_model_dir,
                                                           outputdir=N1_monomer_outdir, generate_egnn_models=True,
                                                           model_count=5)
            if result is None:
                raise RuntimeError(f"Program failed in step 7: monomer {monomer_id} model evaluation")

            if FLAGS.refine:
                makedir_if_not_exists(N2_outdir)

                print("2. Start to refine monomer models based on the qa rankings")

                N2_monomer_outdir = N2_outdir + '/' + monomer_id

                makedir_if_not_exists(N2_monomer_outdir)

                os.system(f"cp {result['casp13']} {N2_monomer_outdir}")

                ref_ranking = read_qa_txt_as_df(result['casp13'])  # apollo or average ranking or the three qas

                refine_inputs = []
                for i in range(5):
                    pdb_name = ref_ranking.loc[i, 'model'] + '.pdb'
                    refine_input = iterative_refine_pipeline.refinement_input(fasta_path=FLAGS.fasta_path,
                                                                              pdb_path=N1_monomer_outdir + '/pdb/' + pdb_name,
                                                                              pkl_path=N1_monomer_outdir + '/pkl/' + pdb_name.replace(
                                                                                  '.pdb', '.pkl'),
                                                                              msa_path=N1_monomer_outdir + '/msa/' + pdb_name.replace(
                                                                                  '.pdb', '.a3m'))
                    refine_inputs += [refine_input]

                final_dir = N2_monomer_outdir + '_final'
                run_monomer_refinement_pipeline_human(params=params, ranking_df=ref_ranking,
                                                      refinement_inputs=refine_inputs, outdir=N2_monomer_outdir,
                                                      finaldir=final_dir, prefix="refine",
                                                      refpdb=N1_monomer_outdir+'/pdb/Elo_0.pdb')


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
