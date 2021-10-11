###################################################################
#Script to generate dimers' msas using interaction calculated by genomic distance
#Dependency: ena_genome_location_table and uniprot_to_embl_table
###################################################################

import os, sys, argparse, time, ftplib, re
from util import is_dir, is_file, read_option_file, check_dirs, check_contents, die, makedir_if_not_exists, clean_dir, cal_sequence_identity_from_seq
from multiprocessing import Pool
from collections import defaultdict, OrderedDict
import pandas as pd
from complex import extract_cds_ids, extract_embl_annotation, add_full_header, find_possible_partners, best_reciprocal_matching, write_concatenated_alignment
from alignment import read_fasta, read_a3m, write_a3m, clean_a3m, cal_identities_to_target


def load_monomer_info(name, a3m_file, uniprot_to_embl_table, ena_genome_location_table):

    cds_ids = extract_cds_ids(a3m_file, uniprot_to_embl_table)

    # print(cds_ids)
    # extract genome location information from ENA

    genome_location_table = extract_embl_annotation(
        cds_ids,
        ena_genome_location_table)

    # print(genome_location_table)

    genome_location_table = add_full_header(
        genome_location_table, a3m_file
    )

    # print(genome_location_table)
    return genome_location_table


def generate_alignments_for_heterodimers(inparams):

    params, chain1, chain2, monomers_msa_dir, outdir = inparams

    uniprot_to_embl_table = params['uniprot_to_embl_table']
    ena_genome_location_table = params['ena_genome_location_table']

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

    gene_location_table_1 = load_monomer_info(chain1, chain1_a3m, uniprot_to_embl_table, ena_genome_location_table)

    gene_location_table_2 = load_monomer_info(chain2, chain2_a3m, uniprot_to_embl_table, ena_genome_location_table)

    print(gene_location_table_1)
    print(gene_location_table_2)

    possible_partners = None
    if gene_location_table_1.empty or gene_location_table_2.empty:
        id_pairing = pd.DataFrame({'id_1': [chain1], 'id_2': [chain2]})
    else:
        # find all possible matches
        possible_partners = find_possible_partners(
            gene_location_table_1, gene_location_table_2
        )

        # find the best reciprocal matches
        id_pairing_unfiltered = best_reciprocal_matching(possible_partners)

        # filter best reciprocal matches by genome distance threshold
        if params["genome_distance_threshold"]:
            distance_threshold = int(params["genome_distance_threshold"])
            id_pairing = id_pairing_unfiltered[id_pairing_unfiltered.distance.astype('int64') < distance_threshold]
        else:
            id_pairing = id_pairing_unfiltered

        id_pairing.loc[:, "id_1"] = id_pairing.loc[:, "uniprot_id_1"]
        id_pairing.loc[:, "id_2"] = id_pairing.loc[:, "uniprot_id_2"]

    print(id_pairing)
    # write concatenated alignment with distance filtering
    # TODO: save monomer alignments?
    target_header, target_seq_idx, sequences_full, sequences_monomer_1, sequences_monomer_2 = \
        write_concatenated_alignment(id_pairing, chain1_a3m, chain2_a3m, chain1, chain2)

    # save the alignment files
    mon_alignment_file_1 = f"{outdir}/{chain1}/{chain1}_monomer_1.a3m"
    with open(mon_alignment_file_1, "w") as of:
        write_a3m(sequences_monomer_1, of)

    gene_location_table_1.to_csv(f"{outdir}/{chain1}/{chain1}_gene_loc.table", index=False)

    mon_alignment_file_2 = f"{outdir}/{chain2}/{chain2}_monomer_2.a3m"
    with open(mon_alignment_file_2, "w") as of:
        write_a3m(sequences_monomer_2, of)

    gene_location_table_2.to_csv(f"{outdir}/{chain2}/{chain2}_gene_loc.table", index=False)

    if possible_partners is not None:
        possible_partners.to_csv(f"{outdir}/{chain1}_{chain2}_partners.csv", index=False)

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

    check_contents(params, ['uniprot_to_embl_table', 'ena_genome_location_table'])

    original_work_dir = params['complex_set_work_dir'] + '/original'

    current_work_dir = params['complex_set_work_dir']

    # generate msas for heterodimers

    os.system(f"cp {original_work_dir}/heterodimers.list {current_work_dir}")

    heterodimers_list = f"{current_work_dir}/heterodimers.list"

    contents = open(heterodimers_list, "r").readlines()

    makedir_if_not_exists(current_work_dir + '/fr_dimers/geno_dist')

    for content in contents:
        content = content.replace('\n', '').replace('\r', '')

        pdbcode, chain1, chain2 = content.split()

        generate_alignments_for_heterodimers([params, chain1, chain2, current_work_dir + '/fr',
                                              current_work_dir + '/fr_dimers/geno_dist'])

    t2 = time.time()

    print("Fr library updating is done. Total time:" + (t2 - t1).__str__())




