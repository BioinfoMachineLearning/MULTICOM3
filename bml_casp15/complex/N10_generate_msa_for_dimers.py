###################################################################
#Script to generate dimers' msas using interaction calculated by species from the uniprot information
#Dependency: HHblits_db (UniRef)
###################################################################


import os, sys, argparse, time, ftplib
from util import is_dir, is_file, read_option_file, check_dirs, check_contents, die, makedir_if_not_exists, clean_dir, cal_sequence_identity_from_seq
from multiprocessing import Pool
from collections import OrderedDict
import re
import pandas as pd
from complex import extract_header_annotation, read_species_annotation_table, most_similar_by_organism, write_concatenated_alignment
from alignment import read_a3m, write_a3m, clean_a3m, cal_identities_to_target


def load_monomer_info(name, a3m_file):

    annotation = extract_header_annotation(a3m_file)

    annotation_table = read_species_annotation_table(annotation)

    if annotation_table is None:
        return None, None

    seqs = read_a3m(a3m_file)

    similarities = cal_identities_to_target(seqs, seqs[name])

    most_similar_in_species = most_similar_by_organism(similarities, annotation_table)

    return annotation_table, most_similar_in_species


def generate_alignments_for_heterodimers(inparams):

    params, chain1, chain2, monomers_msa_dir, outdir = inparams

    chain1_a3m = monomers_msa_dir + '/' + chain1 + '.a3m'
    chain2_a3m = monomers_msa_dir + '/' + chain2 + '.a3m'

    if not os.path.exists(chain1_a3m):
        #die(f"Cannot find a3m file for {chain1}: {chain1_a3m}\n")
        return

    if not os.path.exists(chain2_a3m):
        #die(f"Cannot find a3m file for {chain2}: {chain2_a3m}\n")
        return

    print(f"Processing {chain1}_{chain2}")
    outdir = f"{outdir}/{chain1}_{chain2}"
    makedir_if_not_exists(outdir)

    makedir_if_not_exists(outdir + '/' + chain1)
    clean_a3m(chain1_a3m, outdir + '/' + chain1 + '/' + chain1 + '.a3m')

    makedir_if_not_exists(outdir + '/' + chain2)
    clean_a3m(chain2_a3m, outdir + '/' + chain2 + '/' + chain2 + '.a3m')

    chain1_a3m = outdir + '/' + chain1 + '/' + chain1 + '.a3m'
    chain2_a3m = outdir + '/' + chain2 + '/' + chain2 + '.a3m'

    annotation_table_1, most_similar_in_species_1 = load_monomer_info(chain1, chain1_a3m)
    print(annotation_table_1)
    print(most_similar_in_species_1)
    annotation_table_2, most_similar_in_species_2 = load_monomer_info(chain2, chain2_a3m)
    print(annotation_table_2)
    print(most_similar_in_species_2)
    if annotation_table_1 is None or annotation_table_2 is None:
        species_intersection = pd.DataFrame({'id_1': [chain1], 'id_2': [chain2]})
    else:
        species_intersection = most_similar_in_species_1.merge(
            most_similar_in_species_2,
            how="inner",  # takes the intersection
            on="species",  # merges on species identifiers
            suffixes=("_1", "_2")
        )

    # write concatenated alignment with distance filtering
    # TODO: save monomer alignments?
    target_header, target_seq_idx, sequences_full, sequences_monomer_1, sequences_monomer_2 = \
        write_concatenated_alignment(species_intersection, chain1_a3m, chain2_a3m, chain1, chain2)

    # save the alignment files
    mon_alignment_file_1 = f"{outdir}/{chain1}/{chain1}_monomer_1.a3m"
    with open(mon_alignment_file_1, "w") as of:
        write_a3m(sequences_monomer_1, of)

    if annotation_table_1 is not None:
        annotation_table_1.to_csv(f"{outdir}/{chain1}/{chain1}_annotation.table", index=False)
        most_similar_in_species_1.to_csv(f"{outdir}/{chain1}/{chain1}_species.csv", index=False)

    mon_alignment_file_2 = f"{outdir}/{chain2}/{chain2}_monomer_2.a3m"
    with open(mon_alignment_file_2, "w") as of:
        write_a3m(sequences_monomer_2, of)

    if annotation_table_2 is not None:
        annotation_table_2.to_csv(f"{outdir}/{chain2}/{chain2}_annotation.table", index=False)
        most_similar_in_species_2.to_csv(f"{outdir}/{chain2}/{chain2}_species.csv", index=False)

    species_intersection.to_csv(f"{outdir}/{chain1}_{chain2}_interact.csv")
    print(species_intersection)
    complex_ailgnment_file = f"{outdir}/{chain1}_{chain2}.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)


def test(params):

    generate_alignments_for_heterodimers([params, '1A5XA', '1A62A', '/home/jl4mc/multicom4s/multicom4s/tool/update_complex_dbs/a3m'])

    die("111111111")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    args = parser.parse_args()

    t1 = time.time()

    params = read_option_file(args.option_file)

    abs_opt = os.path.abspath(args.option_file)

    # test(params)

    check_dirs(params, ['prosys_dir', 'complex_set_work_dir'])

    check_contents(params, ['monomers_cm_database_name', 'monomers_in_complex_database_name', 'monomers_not_in_complex_database_name'])

    original_work_dir = params['complex_set_work_dir'] + '/original'

    current_work_dir = params['complex_set_work_dir']

    # generate msas for heterodimers

    os.system(f"cp {original_work_dir}/heterodimers.list {current_work_dir}")

    heterodimers_list = f"{current_work_dir}/heterodimers.list"

    contents = open(heterodimers_list, "r").readlines()

    makedir_if_not_exists(current_work_dir + '/fr_dimers/best_hit')

    for content in contents:

        content = content.replace('\n', '').replace('\r', '')

        pdbcode, chain1, chain2 = content.split()

        generate_alignments_for_heterodimers([params, chain1, chain2, current_work_dir + '/fr',
                                             current_work_dir + '/fr_dimers/best_hit'])

    t2 = time.time()

    print("Fr library updating is done. Total time:" + (t2 - t1).__str__())




