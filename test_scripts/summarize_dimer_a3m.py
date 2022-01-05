import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.complex_alignment_generation.pipeline import *
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer
from bml_casp15.monomer_alignment_generation.alignment import *

def cal_msa_depth(msa_file):
    #print(msa_file)
    format = "a3m"
    if msa_file.find('.sto') > 0:
        format = "stockholm"
    with open(msa_file) as f:
        aln = Alignment.from_file(f, format=format)
    depth = len(aln.seqs)
    return depth


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
    #print("Start to generate alignments for targets")
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

            alignment = {'chain1': chain1_a3ms,
                         'chain2': chain2_a3ms,
                         'outdir': f"{args.output_dir}/{chain1}_{chain2}",
                         'fused': f"{args.output_dir}/{chain1}_{chain2}/fused.a3m",
                         'uniclust_oxmatch_a3m': f"{args.output_dir}/{chain1}_{chain2}/uniclust_oxmatch_a3m/{chain1}_{chain2}.a3m",
                         'pdb_interact_uniref_a3m': f"{args.output_dir}/{chain1}_{chain2}/pdb_interact_uniref_a3m/{chain1}_{chain2}.a3m",
                         'species_interact_uniref_a3m': f"{args.output_dir}/{chain1}_{chain2}/species_interact_uniref_a3m/{chain1}_{chain2}.a3m",
                         'uniprot_distance_uniref_a3m': f"{args.output_dir}/{chain1}_{chain2}/uniprot_distance_uniref_a3m/{chain1}_{chain2}.a3m",
                         'string_interact_uniref_a3m': f"{args.output_dir}/{chain1}_{chain2}/string_interact_uniref_a3m/{chain1}_{chain2}.a3m",
                         'geno_dist_uniref_a3m': f"{args.output_dir}/{chain1}_{chain2}/geno_dist_uniref_a3m/{chain1}_{chain2}.a3m",
                         'pdb_interact_uniref_sto': f"{args.output_dir}/{chain1}_{chain2}/pdb_interact_uniref_sto/{chain1}_{chain2}.a3m",
                         'species_interact_uniref_sto': f"{args.output_dir}/{chain1}_{chain2}/species_interact_uniref_sto/{chain1}_{chain2}.a3m",
                         'uniprot_distance_uniref_sto': f"{args.output_dir}/{chain1}_{chain2}/uniprot_distance_uniref_sto/{chain1}_{chain2}.a3m",
                         'string_interact_uniref_sto': f"{args.output_dir}/{chain1}_{chain2}/string_interact_uniref_sto/{chain1}_{chain2}.a3m",
                         'geno_dist_uniref_sto': f"{args.output_dir}/{chain1}_{chain2}/geno_dist_uniref_sto/{chain1}_{chain2}.a3m",
                         'pdb_interact_uniprot_sto': f"{args.output_dir}/{chain1}_{chain2}/pdb_interact_uniprot_sto/{chain1}_{chain2}.a3m",
                         'species_interact_uniprot_sto': f"{args.output_dir}/{chain1}_{chain2}/species_interact_uniprot_sto/{chain1}_{chain2}.a3m",
                         'uniprot_distance_uniprot_sto': f"{args.output_dir}/{chain1}_{chain2}/uniprot_distance_uniprot_sto/{chain1}_{chain2}.a3m",
                         'string_interact_uniprot_sto': f"{args.output_dir}/{chain1}_{chain2}/string_interact_uniprot_sto/{chain1}_{chain2}.a3m",
                         'geno_dist_uniprot_sto': f"{args.output_dir}/{chain1}_{chain2}/geno_dist_uniprot_sto/{chain1}_{chain2}.a3m",
                         }

            add = True
            for key in alignment:
                if key.find('uni') >= 0:
                    if not os.path.exists(alignment[key]):
                        add = False

            if add:
                alignments += [alignment]
                depths = []
                output_uniref_a3m = False
                output_uniref_sto = False
                output_uniprot_sto = False
                for key in alignment:
                    #print(key)
                    if key == "fused":
                        contents = open(alignment[key]).readlines()
                        depths += [str(int(len(contents)/2))]

                    if key.find('uniref_a3m') > 0 and not output_uniref_a3m:
                        depth1 = cal_msa_depth(chain1_a3ms['uniref30_a3m'])
                        depth2 = cal_msa_depth(chain2_a3ms['uniref30_a3m'])
                        depths += [f"{depth1}:{depth2}"]
                        output_uniref_a3m = True

                    if key.find('uniref_sto') > 0 and not output_uniref_sto:
                        depth1 = cal_msa_depth(chain1_a3ms['uniref90_sto'])
                        depth2 = cal_msa_depth(chain2_a3ms['uniref90_sto'])
                        depths += [f"{depth1}:{depth2}"]
                        output_uniref_sto = True


                    if key.find('uniprot_sto') > 0 and not output_uniprot_sto:
                        depth1 = cal_msa_depth(chain1_a3ms['uniprot_sto'])
                        depth2 = cal_msa_depth(chain2_a3ms['uniprot_sto'])
                        depths += [f"{depth1}:{depth2}"]
                        output_uniprot_sto = True


                    if key.find('uni') >= 0:
                        depths += [str(cal_msa_depth(alignment[key]))]

                depths_str = "\t".join(depths)
                print(f"{chain1}_{chain2}\t{depths_str}")



    #print(f"Total {len(alignments)} pairs can be concatenated")

    #print("Start to do summary for dimers")





