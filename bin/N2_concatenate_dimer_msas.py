import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import check_file, check_dir, makedir_if_not_exists, check_contents, read_option_file
from bml_casp15.complex_alignment_generation.pipeline import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('aln_dir', None, 'Monomer alignment directory')
flags.DEFINE_list('dimerlist', None, 'List of dimers: A_B,B_C')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)
    check_dir(FLAGS.aln_dir)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    print("Start to generate alignments for targets")
    alignments = []

    for dimer in FLAGS.dimerlist:
        chain1, chain2 = dimer.split('_')
        chain1_aln_dir = FLAGS.aln_dir + '/' + chain1
        if os.path.exists(chain1_aln_dir):
            chain1_a3ms = {'name': chain1,
                           'uniref30_a3m': f"{chain1_aln_dir}/{chain1}_uniref30.a3m",
                           'uniref90_sto': f"{chain1_aln_dir}/{chain1}_uniref90.sto",
                           'uniprot_sto': f"{chain1_aln_dir}/{chain1}_uniprot.sto",
                           'uniclust30_a3m': f"{chain1_aln_dir}/{chain1}_uniclust30.a3m"}
                           #'mgnify_sto': f"{chain1_aln_dir}/{chain1}_mgnify.sto",
                           #'bfd_a3m': f"{chain1_aln_dir}/{chain1}_bfd.a3m"}
        else:
            chain1_a3ms = {'name': chain1}

        chain2_aln_dir = FLAGS.aln_dir + '/' + chain2
        if os.path.exists(chain2_aln_dir):
            chain2_a3ms = {'name': chain2,
                           'uniref30_a3m': f"{chain2_aln_dir}/{chain2}_uniref30.a3m",
                           'uniref90_sto': f"{chain2_aln_dir}/{chain2}_uniref90.sto",
                           'uniprot_sto': f"{chain2_aln_dir}/{chain2}_uniprot.sto",
                           'uniclust30_a3m': f"{chain2_aln_dir}/{chain2}_uniclust30.a3m"}
                           #'mgnify_sto': f"{chain2_aln_dir}/{chain2}_mgnify.sto",
                           #'bfd_a3m': f"{chain2_aln_dir}/{chain2}_bfd.a3m"}
        else:
            chain2_a3ms = {'name': chain2}

        alignment = {'chain1': chain1_a3ms, 'chain2': chain2_a3ms, 'outdir': f"{FLAGS.output_dir}/{chain1}_{chain2}"}

        add = True
        for a3ms in [chain1_a3ms, chain2_a3ms]:
            if len(a3ms) == 1:
                add = False
            for key in a3ms:
                if key.find('uni') >= 0:
                    if not os.path.exists(a3ms[key]):
                        add = False
                    else:
                        contents = open(a3ms[key]).readlines()
                        if len(contents) == 0:
                            print(f"File: {a3ms[key]} is empty!")
                            add = False
                            os.system(f"rm {a3ms[key]}")

        if add:
            alignments += [alignment]

    print(f"Total {len(alignments)} pairs can be concatenated")

    if len(alignments) > 0:
        print("Start to concatenate alignments for dimers")
        pipeline = Complex_alignment_concatenation_pipeline(params, multiprocess=True, process_num=5)
        alignments = pipeline.concatenate(alignments, params['hhfilter_program'])
    print("The alignment concatenation for dimers has finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'aln_dir',
        'dimerlist',
        'output_dir'
    ])
    app.run(main)
