import os, sys, argparse, time
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import read_fasta, write_fasta
from bml_casp15.monomer_alignment_generation.pipeline import Monomer_alignment_generation_pipeline
from bml_casp15.complex_alignment_generation.pipeline import *
from bml_casp15.tertiary_structure_generation.pipeline import *
from bml_casp15.complex_templates_search.structure_based_pipeline import *
from bml_casp15.complex_templates_search.sequence_based_pipeline import *
from bml_casp15.quaternary_structure_generation.pipeline import *
from absl import flags
from absl import app

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_list('fasta_paths', None, 'Path to multimer fastas')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS


def run_monomer_msa_pipeline(inparams):
    fasta, outdir, params = inparams

    uniref30 = params['uniref_db_dir'] + '/' + params['uniref_db']
    uniclust30 = params['uniclust_db_dir'] + '/' + params['uniclust_db']
    uniref90_fasta = params['uniref90_fasta']
    uniprot_fasta = params['uniprot_fasta']
    smallbfd = ""  # params['smallbfd_database']
    bfd = params['bfd_database']
    mgnify = params['mgnify_database']

    hhblits_binary = params['hhblits_program']
    jackhmmer_binary = params['jackhmmer_program']

    result = None
    try:
        monomer_msa_pipeline = Monomer_alignment_generation_pipeline(jackhmmer_binary,
                                                                     hhblits_binary,
                                                                     uniref90_fasta,
                                                                     mgnify,
                                                                     smallbfd,
                                                                     bfd,
                                                                     uniref30,
                                                                     uniclust30,
                                                                     uniprot_fasta)

        result = monomer_msa_pipeline.process(fasta, outdir)
    except Exception as e:
        return result
    return result


def run_concatenate_dimer_msas_pipeline(dimer, monomer_aln_dir, outputdir, params):
    chain1, chain2 = dimer.split('_')
    chain1_aln_dir = monomer_aln_dir + '/' + chain1
    if os.path.exists(chain1_aln_dir):
        chain1_a3ms = {'name': chain1,
                       'uniref30_a3m': f"{chain1_aln_dir}/{chain1}_uniref30.a3m",
                       'uniref90_sto': f"{chain1_aln_dir}/{chain1}_uniref90.sto",
                       'uniprot_sto': f"{chain1_aln_dir}/{chain1}_uniprot.sto",
                       'uniclust30_a3m': f"{chain1_aln_dir}/{chain1}_uniclust30.a3m"}
    else:
        chain1_a3ms = {'name': chain1}

    chain2_aln_dir = monomer_aln_dir + '/' + chain2
    if os.path.exists(chain2_aln_dir):
        chain2_a3ms = {'name': chain2,
                       'uniref30_a3m': f"{chain2_aln_dir}/{chain2}_uniref30.a3m",
                       'uniref90_sto': f"{chain2_aln_dir}/{chain2}_uniref90.sto",
                       'uniprot_sto': f"{chain2_aln_dir}/{chain2}_uniprot.sto",
                       'uniclust30_a3m': f"{chain2_aln_dir}/{chain2}_uniclust30.a3m"}
    else:
        chain2_a3ms = {'name': chain2}

    alignment = {'chain1': chain1_a3ms, 'chain2': chain2_a3ms, 'outdir': outputdir}

    complete = True
    for a3ms in [chain1_a3ms, chain2_a3ms]:
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
                        # os.system(f"rm {a3ms[key]}")

    if complete:

        alignments = [alignment]

        print(f"Total {len(alignments)} pairs can be concatenated")

        print("Start to concatenate alignments for dimers")

        complex_alignment_concatenation_pipeline = Complex_alignment_concatenation_pipeline(params)

        alignments = complex_alignment_concatenation_pipeline.concatenate(alignments, params['hhfilter_program'])

    else:

        print("The a3ms for dimers are not complete!")


def complete_result(outputdir):
    complete = True
    for i in range(0, 5):
        model = f'{outputdir}/ranked_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    for fasta_path in FLAGS.fasta_paths:

        check_file(fasta_path)

        print(fasta_path)

        dimer_name = os.path.basename(fasta_path[0:fasta_path.index('.fasta')])

        dimer_outdir = FLAGS.output_dir + '/' + dimer_name

        makedir_if_not_exists(dimer_outdir)

        N1_outdir = dimer_outdir + '/N1_monomer_alignments_generation'
        makedir_if_not_exists(N1_outdir)
        fasta_paths = []
        monomer_list = []
        with open(fasta_path) as fileobj:
            for seq_id, seq in read_fasta(fileobj):
                with open(f"{dimer_outdir}/{seq_id}.fasta", "w") as fw:
                    write_fasta({seq_id: seq}, fw)

                outdir = N1_outdir + '/' + seq_id
                makedir_if_not_exists(outdir)
                monomer_list.append([f"{dimer_outdir}/{seq_id}.fasta", outdir, params])
                fasta_paths += [f"{dimer_outdir}/{seq_id}.fasta"]

        print("#################################################################################################")
        print(f"1. Start to generate alignments for monomers, total {len(monomer_list)} monomers to be processed")
        pool = Pool(processes=5)
        results = pool.map(run_monomer_msa_pipeline, monomer_list)
        pool.close()
        pool.join()

        for result in results:
            if result is None:
                raise RuntimeError('Some of the monomers alignment generation has failed!')

        print("#################################################################################################")
        print("2. Start to generate alignments for targets")
        N2_outdir = dimer_outdir + '/N2_complex_alignments_concatenation'
        makedir_if_not_exists(N2_outdir)

        try:
            run_concatenate_dimer_msas_pipeline(dimer_name, N1_outdir, N2_outdir, params)
        except Exception as e:
            print(e)
            print("Program failed in step 2")

        print("#################################################################################################")

        print("#################################################################################################")
        print("3. Start to generate tertiary structure for monomers using alphafold")
        N3_outdir = dimer_outdir + '/N3_monomer_tertiary_structure_generation'
        makedir_if_not_exists(N3_outdir)

        to_be_processed_fasta = []
        if len(fasta_paths) > 0:
            for fasta_path in fasta_paths:
                outdir = N3_outdir + '/' + pathlib.Path(fasta_path).stem
                if not complete_result(outdir):
                    to_be_processed_fasta += [fasta_path]

        if len(to_be_processed_fasta) > 0:
            print(f"Total {len(to_be_processed_fasta)} monomers are generating structures")
            try:
                pipeline = Monomer_tertiary_structure_prediction_pipeline(params)
                pipeline.process(to_be_processed_fasta, N1_outdir, N3_outdir)
            except Exception as e:
                print(e)
                print("Program failed in step 3")

        print("The prediction for monomers has finished!")

        print("#################################################################################################")

        print("#################################################################################################")
        print("4. Start to search complex templates based on monomer structures")
        N4_outdir = dimer_outdir + '/N4_structure_based_templates_search'
        makedir_if_not_exists(N4_outdir)

        monomer1, monomer2 = dimer_name.split('_')

        monomer1_pdb = f"{N3_outdir}/{monomer1}/ranked_0.pdb"

        if not os.path.exists(monomer1_pdb):
            print(f"Cannot find teritary structure for {monomer1}: {monomer1_pdb}")
            continue

        monomer2_pdb = f"{N3_outdir}/{monomer2}/ranked_0.pdb"
        if not os.path.exists(monomer2_pdb):
            print(f"Cannot find teritary structure for {monomer2}: {monomer2_pdb}")
            continue

        os.system(f"cp {monomer1_pdb} {N4_outdir}/{monomer1}.pdb")
        monomer1_pdb = f"{N4_outdir}/{monomer1}.pdb"

        os.system(f"cp {monomer2_pdb} {N4_outdir}/{monomer2}.pdb")
        monomer2_pdb = f"{N4_outdir}/{monomer2}.pdb"

        try:
            pipeline = Complex_structure_based_template_search_pipeline(params)

            pipeline.search([monomer1_pdb, monomer2_pdb], N4_outdir)
        except Exception as e:
            print(e)
            print("Program failed in step 4")

        print("Complex structure based template searching has been finished!")

        print("#################################################################################################")

        print("#################################################################################################")
        print("5. Start to search complex templates based on sequence similarity")
        N5_outdir = dimer_outdir + '/N5_sequence_based_templates_search'
        makedir_if_not_exists(N5_outdir)

        monomer1_template_a3m = f"{N1_outdir}/{monomer1}/{monomer1}_uniref90.sto"

        monomer2_template_a3m = f"{N1_outdir}/{monomer2}/{monomer2}_uniref90.sto"

        if not os.path.exists(monomer1_template_a3m):
            continue

        if not os.path.exists(monomer2_template_a3m):
            continue

        seq1 = open(f"{N1_outdir}/{chain1}/{chain1}.fasta").readlines()[1].rstrip('\n')
        seq2 = open(f"{N1_outdir}/{chain2}/{chain2}.fasta").readlines()[1].rstrip('\n')

        monomer1_template_input = monomer_template_input(name=monomer1, msa_path=monomer1_template_a3m, seq=seq1)

        monomer2_template_input = monomer_template_input(name=monomer2, msa_path=monomer2_template_a3m, seq=seq2)

        try:
            pipeline = Complex_sequence_based_template_search_pipeline(params)
            pipeline.search([monomer1_template_input, monomer2_template_input], N5_outdir)
        except Exception as e:
            print(e)
            print("Program failed in step 5")

        print("Complex sequence based template searching has been finished!")

        print("#################################################################################################")

        print("#################################################################################################")
        print("6. Start to generate complex quaternary structures")
        N6_outdir = dimer_outdir + '/N6_quaternary_structure_generation'
        makedir_if_not_exists(N6_outdir)

        try:
            pipeline = Quaternary_structure_prediction_pipeline(params)
            result = pipeline.process(fasta_path=fasta_path,
                                      aln_dir=N1_outdir,
                                      complex_aln_dir=N2_outdir,
                                      structure_template_dir=N4_outdir,
                                      sequence_template_dir=N5_outdir,
                                      monomer_model_dir=N3_outdir,
                                      output_dir=N6_outdir)
        except Exception as e:
            print(e)
            print("Program failed in step 6")

        print("Complex quaternary structure generation has been finished!")

        print("#################################################################################################")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_paths',
        'output_dir'
    ])
    app.run(main)
