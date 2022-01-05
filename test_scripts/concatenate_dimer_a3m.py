import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.complex_alignment_generation.pipeline import *
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--dimerlist', type=is_file, required=True)
    parser.add_argument('--aln_dir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    makedir_if_not_exists(args.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    process_list = []
    print("Start to generate alignments for targets")
    alignments = []
    with open(args.dimerlist) as f:
        for line in f:
            chain1, chain2 = line.rstrip('\n').split()

            chain1_aln_dir = args.aln_dir + '/' + chain1
            if os.path.exists(chain1_aln_dir):
                chain1_a3ms = {'name': chain1,
                               'uniref30_a3m': f"{chain1_aln_dir}/{chain1}_uniref30.a3m",
                               'uniref90_sto': f"{chain1_aln_dir}/{chain1}_uniref90.sto",
                               'uniprot_sto': f"{chain1_aln_dir}/{chain1}_uniprot.sto",
                               'uniclust30_a3m': f"{chain1_aln_dir}/{chain1}_uniclust30.a3m",
                               'mgnify_sto': f"{chain1_aln_dir}/{chain1}_mgnify.sto",
                               'bfd_a3m': f"{chain1_aln_dir}/{chain1}_bfd.a3m"}
            else:
                chain1_a3ms = {'name': chain1}

            chain2_aln_dir = args.aln_dir + '/' + chain2
            if os.path.exists(chain2_aln_dir):
                chain2_a3ms = {'name': chain2,
                               'uniref30_a3m': f"{chain2_aln_dir}/{chain2}_uniref30.a3m",
                               'uniref90_sto': f"{chain2_aln_dir}/{chain2}_uniref90.sto",
                               'uniprot_sto': f"{chain2_aln_dir}/{chain2}_uniprot.sto",
                               'uniclust30_a3m': f"{chain2_aln_dir}/{chain2}_uniclust30.a3m",
                               'mgnify_sto': f"{chain2_aln_dir}/{chain2}_mgnify.sto",
                               'bfd_a3m': f"{chain2_aln_dir}/{chain2}_bfd.a3m"}
            else:
                chain2_a3ms = {'name': chain2}

            alignment = {'chain1': chain1_a3ms, 'chain2': chain2_a3ms, 'outdir': f"{args.output_dir}/{chain1}_{chain2}"}

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
                                os.system(f"rm {a3ms[key]}")

            if add:
                alignments += [alignment]

    print(f"Total {len(alignments)} pairs can be concatenated")

    print("Start to concatenate alignments for dimers")

    fast_methods = ["fused", "pdb_interact", "species_interact", "uniclust_oxmatch", "uniprot_distance"]
    pipeline = Complex_alignment_concatenation_pipeline(params, concatenate_methods=fast_methods, multiprocess=True, process_num=5)
    alignments = pipeline.concatenate(alignments, params['hhfilter_program'])

    # slow_methods = ["string_interact"]
    # pipeline = Complex_alignment_concatenation_pipeline(params, concatenate_methods=slow_methods, multiprocess=True,
    #                                                     process_num=5)
    # alignments = pipeline.concatenate(alignments, params['hhfilter_program'])



