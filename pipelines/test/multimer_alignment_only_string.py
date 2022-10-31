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
from bml_casp15.complex_alignment_generation.pipeline_v3 import *

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

    alignments = []

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

        chains = [chain_id_map[chain_id].description for chain_id in chain_id_map]
        # alignment = {'outdir': f"{outputdir}/{'_'.join(chains)}"}
        alignment = {'outdir': N4_outdir}
        for i in range(len(chains)):
            chain = chains[i]
            chain_aln_dir = N1_outdir + '/' + chain
            if os.path.exists(chain_aln_dir):
                chain_a3ms = {'name': chain,
                            'colabfold_a3m': f"{chain_aln_dir}/{chain}_colabfold.a3m",
                            'uniref30_a3m': f"{chain_aln_dir}/{chain}_uniref30.a3m",
                            'uniref90_sto': f"{chain_aln_dir}/{chain}_uniref90.sto",
                            'uniprot_sto': f"{chain_aln_dir}/{chain}_uniprot.sto",
                            'uniclust30_a3m': f"{chain_aln_dir}/{chain}_uniclust30.a3m"}
            else:
                chain_a3ms = {'name': chain}
            alignment[f"chain{i + 1}"] = chain_a3ms

        complete = True
        for name in alignment:
            if name == 'outdir':
                continue
            a3ms = alignment[name]
            if len(a3ms) == 1:
                complete = False
            for key in a3ms:
                if key.find('uni') >= 0:
                    if not os.path.exists(a3ms[key]):
                        complete = False
                    else:
                        contents = open(a3ms[key]).readlines()
                        if len(contents) == 0:
                            print(f"File: {a3ms[key]} is empty!")
                            complete = False
                            os.system(f"rm {a3ms[key]}")

        if complete:
            alignments += [alignment]
        else:
            print("The a3ms for dimers are not complete!")

    print(f"Total {len(alignments)} pairs can be concatenated")

    print("Start to concatenate alignments for dimers")

    concat_methods = ['string_interact']

    complex_alignment_concatenation_pipeline = Complex_alignment_concatenation_pipeline(params=params, run_methods=concat_methods)
    alignments = complex_alignment_concatenation_pipeline.concatenate(alignments, params['hhfilter_program'], is_homomers=False)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fastadir',
        'output_dir'
    ])
    app.run(main)
