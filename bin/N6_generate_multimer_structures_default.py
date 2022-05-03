import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.quaternary_structure_generation.pipeline_default import *
from absl import flags
from absl import app
from bml_casp15.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map


flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'option file')
flags.DEFINE_string('aln_dir', None, 'Monomer model directory')
flags.DEFINE_string('output_dir', None, 'Monomer model directory')
FLAGS = flags.FLAGS


def run_pipeline(inparams):
    params, fasta_path, aln_dir, output_dir = inparams

    with open(fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    pipeline = Quaternary_structure_prediction_pipeline(params)

    result = None
    try:
        result = pipeline.process(fasta_path=fasta_path,
                                  chain_id_map=chain_id_map,
                                  aln_dir=aln_dir,
                                  output_dir=output_dir)
    except Exception as e:
        print(e)
        return result
    return result


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    for fasta_path in os.listdir(FLAGS.fasta_path):
        fasta_path = FLAGS.fasta_path + '/' + fasta_path
        monomers = []
        for line in open(fasta_path):
            line = line.rstrip('\n').strip()
            if line.startswith('>'):
                monomers += [line[1:]]

        output_dir = FLAGS.output_dir + '/' + '_'.join(monomers)
        makedir_if_not_exists(output_dir)

        inparams = [params, fasta_path, FLAGS.aln_dir, output_dir]

        # process_list.append(inparams)

        run_pipeline(inparams)

    print("The template search dimers has finished!")

if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'aln_dir',
        'output_dir'
    ])
    app.run(main)
