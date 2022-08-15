import os, sys, argparse, time, copy, pathlib
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from bml_casp15.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from bml_casp15.common.pipeline import run_multimer_refinement_pipeline
import pandas as pd
from absl import flags
from absl import app
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, clean_dir
from bml_casp15.quaternary_structure_refinement import iterative_refine_pipeline_multimer

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('indir', None, 'Path to multimer fastas')
flags.DEFINE_string('outdir', None, 'Output directory')
FLAGS = flags.FLAGS

def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    makedir_if_not_exists(FLAGS.outdir)

    params = read_option_file(FLAGS.option_file)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    print("1111111111111111111111111111111111111")
    refine_inputs = []
    for pdb_name in os.listdir(FLAGS.indir):
        if pdb_name.find('.pdb') < 0:
            continue
        print(pdb_name)
        msa_paths = {}
        for chain_id in chain_id_map:
            msa_paths[chain_id] = dict(
                paired_msa=f"{FLAGS.indir}/{chain_id_map[chain_id].description}/{pdb_name.replace('.pdb', '')}.paired.a3m",
                monomer_msa=f"{FLAGS.indir}/{chain_id_map[chain_id].description}/{pdb_name.replace('.pdb', '')}.monomer.a3m")
        print(msa_paths)
        refine_input = iterative_refine_pipeline_multimer.refinement_input_multimer(chain_id_map=chain_id_map,
                                                                                    fasta_path=FLAGS.fasta_path,
                                                                                    pdb_path=FLAGS.indir + '/' + pdb_name,
                                                                                    pkl_path=FLAGS.indir + '/' + pdb_name.replace(
                                                                                        '.pdb', '.pkl'),
                                                                                    msa_paths=msa_paths)
        refine_inputs += [refine_input]

    final_dir = FLAGS.outdir + '_final'
    run_multimer_refinement_pipeline(chain_id_map=chain_id_map,
                                     params=params, refinement_inputs=refine_inputs, outdir=FLAGS.outdir,
                                     finaldir=final_dir, stoichiometry="homomer")

    print("The refinement for the top-ranked multimer models has been finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'outdir'
    ])
    app.run(main)
