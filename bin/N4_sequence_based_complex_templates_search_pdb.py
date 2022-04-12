import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.complex_templates_search.sequence_based_pipeline_pdb import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('aln_dir', None, 'Monomer model directory')
flags.DEFINE_list('dimerlist', None, 'List of dimers: A_B,B_C')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def run_pipeline(inparams):
    
    params, monomer_inputs, outdir = inparams

    pipeline = Complex_sequence_based_template_search_pipeline(params)

    pipeline.search(monomer_inputs, outdir)


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    pipeline = Complex_sequence_based_template_search_pipeline(params)

    process_list = []

    for dimer in FLAGS.dimerlist:
        
        chain1, chain2 = dimer.split('_')

        chain1_template_a3m = f"{FLAGS.aln_dir}/{chain1}/{chain1}_uniref90.sto"

        chain2_template_a3m = f"{FLAGS.aln_dir}/{chain2}/{chain2}_uniref90.sto"

        if not os.path.exists(chain1_template_a3m):
            continue

        if not os.path.exists(chain2_template_a3m):
            continue

        seq1 = open(f"{FLAGS.aln_dir}/{chain1}/{chain1}.fasta").readlines()[1].rstrip('\n')
        seq2 = open(f"{FLAGS.aln_dir}/{chain2}/{chain2}.fasta").readlines()[1].rstrip('\n')

        monomer1_template_input = monomer_template_input(name=chain1, msa_path=chain1_template_a3m, hmm_path="", seq=seq1)

        monomer2_template_input = monomer_template_input(name=chain2, msa_path=chain2_template_a3m, hmm_path="", seq=seq2)

        dimer_outdir = f"{FLAGS.output_dir}/{chain1}_{chain2}"

        if not os.path.exists(dimer_outdir + '/concatenated_template_idx.csv'):
            process_list.append([params, [monomer1_template_input, monomer2_template_input], dimer_outdir])

        # run_pipeline([params, [monomer1_template_input, monomer2_template_input], dimer_outdir])

    print(f"Total {len(process_list)} dimers to be processed")

    pool = Pool(processes=15)
    results = pool.map(run_pipeline, process_list)
    pool.close()
    pool.join()

    print("The template search dimers has finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'aln_dir',
        'dimerlist',
        'output_dir'
    ])
    app.run(main)
