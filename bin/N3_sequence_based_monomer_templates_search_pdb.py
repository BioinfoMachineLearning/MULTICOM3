import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.monomer_templates_search.sequence_based_pipeline_pdb import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('aln_dir', None, 'Monomer model directory')
flags.DEFINE_list('fasta_paths', None, 'List of dimers: A,B')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def run_pipeline(inparams):
    
    params, targetname, sequence, a3m, outdir = inparams

    pipeline = monomer_sequence_based_template_search_pipeline(params)

    pipeline.search(targetname=targetname, sequence=sequence, a3m=a3m, outdir=outdir)


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    process_list = []

    for fasta_path in FLAGS.fasta_paths:
        targetname = None
        sequence = None
        for line in open(fasta_path):
            line = line.rstrip('\n').strip()
            if line.startswith('>'):
                targetname = line[1:]
            else:
                sequence = line

        template_a3m = f"{FLAGS.aln_dir}/{targetname}/{targetname}_uniref90.sto"

        if not os.path.exists(template_a3m):
            continue

        outdir = f"{FLAGS.output_dir}/{targetname}"

        if not os.path.exists(f"{outdir}/sequence_templates.csv"):
            process_list.append([params, targetname, sequence, template_a3m, outdir])

            # run_pipeline([params, targetname, sequence, template_a3m, outdir])

    print(f"Total {len(process_list)} dimers to be processed")

    pool = Pool(processes=20)
    results = pool.map(run_pipeline, process_list)
    pool.close()
    pool.join()

    print("The template search dimers has finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'aln_dir',
        'fasta_paths',
        'output_dir'
    ])
    app.run(main)
