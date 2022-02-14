import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
#from bml_casp15.complex_alignment_generation.pipeline.concatenate import Concatenate
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer
from bml_casp15.monomer_alignment_generation.pipeline import Monomer_alignment_generation_pipeline


def generate_a3ms_for_single_seq(inparams):

    fasta, outdir, params = inparams

    uniref30 = params['uniref_db_dir'] + '/' + params['uniref_db']
    uniclust30 = params['uniclust_db_dir'] + '/' + params['uniclust_db']
    uniref90_fasta = params['uniref90_fasta']
    uniprot_fasta = params['uniprot_fasta']
    smallbfd = params['smallbfd_database']
    bfd = params['bfd_database']
    mgnify = params['mgnify_database']

    hhblits_binary = params['hhblits_program']
    jackhmmer_binary = params['jackhmmer_program']

    pipeline = Monomer_alignment_generation_pipeline(jackhmmer_binary,
                                                     hhblits_binary,
                                                     uniref90_fasta,
                                                     mgnify,
                                                     smallbfd,
                                                     bfd,
                                                     uniref30,
                                                     uniclust30,
                                                     uniprot_fasta)

    return pipeline.process(fasta, outdir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--fastadir', type=is_dir, required=True)
    parser.add_argument('--dimerlist', type=is_file, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--finish_list', type=str, required=False)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    makedir_if_not_exists(args.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    process_list = []
    print("Start to generate alignments for monomers")
    monomer_list = []
    with open(args.dimerlist) as f:
        for line in f:
            chain1, chain2 = line.rstrip('\n').split()
            if chain1 not in monomer_list:
                monomer_list += [chain1]
            if chain2 not in monomer_list:
                monomer_list += [chain2]

    finish_list = []
    if args.finish_list is not None:
        finish_list = [line.rstrip('\n') for line in open(args.finish_list).readlines()]

    for monomer in monomer_list:
        if not os.path.exists(f"{args.fastadir}/{monomer}.fasta"):
            print(f"Cannot find fasta file for {monomer}!")
            continue
        if monomer in finish_list:
            continue

        outdir = args.output_dir + f"/{monomer}"
        makedir_if_not_exists(outdir)
        os.system(f"cp {args.fastadir}/{monomer}.fasta {outdir}")
        process_list.append([f"{outdir}/{monomer}.fasta", outdir, params])

    print(f"Total {len(process_list)} monomers to be processed")
    pool = Pool(processes=8)
    results = pool.map(generate_a3ms_for_single_seq, process_list)
    pool.close()
    pool.join()
