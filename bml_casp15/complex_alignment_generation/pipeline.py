import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_dir, is_file, read_option_file, makedir_if_not_exists
from bml_casp15.complex_alignment_generation.geno_dist import Geno_interact
from bml_casp15.complex_alignment_generation.pdb_interact import PDB_interact
from bml_casp15.complex_alignment_generation.species_interact import Species_interact
from bml_casp15.complex_alignment_generation.string_interact import STRING_interact
from bml_casp15.complex_alignment_generation.uniclust_oxmatch import UNICLUST_oxmatch
from bml_casp15.complex_alignment_generation.uniprot_distance import UNIPROT_distance
from bml_casp15.monomer_alignment_generation.alignment import *
from bml_casp15.common.util import makedir_if_not_exists

mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
           'G': 6, 'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12,
           'O': 21, 'P': 13, 'Q': 14, 'R': 15, 'S': 16, 'T': 17,
           'V': 18, 'W': 19, 'Y': 20, 'U': 21, 'Z': 21, 'X': 21, 'J': 21}

backmap = {1: 'A', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H',
           8: 'I', 9: 'K', 10: 'L', 11: 'M', 12: 'N', 13: 'P', 14: 'Q',
           15: 'R', 16: 'S', 17: 'T', 18: 'V', 19: 'W', 20: 'Y', 21: '-'}


def fused_msa(sequences1, sequences2, fuse_msa_name):
    parsed_sequences1 = []
    for sequence in sequences1:
        gap_fraction = sequence.count('-') / float(len(sequence))
        if gap_fraction <= 0.9:  # Only use the lines with less than 90 % gaps
            parsed_sequences1.append([mapping.get(ch, 22) for ch in sequence if not ch.islower()])
    parsed_sequences1 = np.array(parsed_sequences1, dtype=np.int8, order='F')

    parsed_sequences2 = []
    for sequence in sequences2:
        gap_fraction = sequence.count('-') / float(len(sequence))
        if gap_fraction <= 0.9:  # Only use the lines with less than 90 % gaps
            parsed_sequences2.append([mapping.get(ch, 22) for ch in sequence if not ch.islower()])
    parsed_sequences2 = np.array(parsed_sequences2, dtype=np.int8, order='F')

    # Construct entire a3m matrix
    fused = np.zeros((parsed_sequences1.shape[0] + parsed_sequences2.shape[0],
                      parsed_sequences1.shape[1] + parsed_sequences2.shape[1]))
    fused[:] = 21  # Assign gaps
    # Assign a3m1
    fused[:parsed_sequences1.shape[0], :parsed_sequences1.shape[1]] = parsed_sequences1
    # Assign a3m2
    fused[parsed_sequences1.shape[0]:, parsed_sequences1.shape[1]:] = parsed_sequences2
    # Write the fused MSA
    write_fused_a3m(fused, fuse_msa_name)


def write_fused_a3m(fused, outfile):
    '''Write a3m MSA'''
    with open(outfile, 'w') as file:
        for i in range(len(fused)):
            file.write('>' + str(i) + '\n')
            file.write(''.join([backmap[ch] for ch in fused[i]]) + '\n')
    return None


def parse_header(header):
    # discard anything in header that might come after the
    # first whitespace (by convention this is typically annotation)
    header = header.split()[0]

    # try to extract region from sequence header
    m = re.search("(.+)/(\d+)-(\d+)", header)
    if m:
        id_, start_str, end_str = m.groups()
        region_start, region_end = int(start_str), int(end_str)
        return id_, region_start, region_end
    else:
        # cannot find region, so just give back sequence iD
        return header, None, None


def write_concatenated_alignment(id_pairing, alignment_1, alignment_2):
    def _prepare_header(id1, id2):

        return "{}_{}".format(id1, id2)

    sequences_to_write = []  # list of (header,seq1,seq2) tuples

    target_header = _prepare_header(alignment_1.main_id, alignment_2.main_id)

    sequences_to_write.append((target_header, alignment_1.main_seq, alignment_2.main_seq))

    # create other headers and sequences
    seen_seqs = []
    filter_pair_ids = {'id_1': [], 'id_2': [], 'index_1': [], 'index_2': []}
    for id1, id2 in zip(id_pairing.id_1, id_pairing.id_2):

        # prepare the concatenated header
        # print(alignment_1.headers)
        # print(id1)
        combine_seq = alignment_1[id1] + alignment_2[id2]
        if combine_seq in seen_seqs:
            continue

        seen_seqs += [combine_seq]
        header1, _, _ = parse_header(alignment_1.headers[id1][0])
        header2, _, _ = parse_header(alignment_2.headers[id2][0])
        concatenated_header = _prepare_header(header1, header2)

        # save the information
        sequences_to_write.append(
            (
                concatenated_header,
                alignment_1[id1],
                alignment_2[id2]
            )
        )

        filter_pair_ids['id_1'] += [id1]
        filter_pair_ids['id_2'] += [id2]
        filter_pair_ids['index_1'] += [alignment_1.id_to_index[id1] + 1]
        filter_pair_ids['index_2'] += [alignment_2.id_to_index[id2] + 1]

    # concatenate strings
    sequences_full = OrderedDict([
        (header, seq1 + seq2) for header, seq1, seq2 in sequences_to_write
    ])

    sequences_monomer_1 = OrderedDict([
        (header, seq1) for header, seq1, seq2 in sequences_to_write
    ])

    sequences_monomer_2 = OrderedDict([
        (header, seq2) for header, seq1, seq2 in sequences_to_write
    ])

    return target_header, sequences_full, sequences_monomer_1, sequences_monomer_2, pd.DataFrame(filter_pair_ids)


def write_dimer_a3ms(pair_ids, aln_1, aln_2, hhfilter, outdir):
    makedir_if_not_exists(outdir)

    target_header, sequences_full, sequences_monomer_1, sequences_monomer_2, pair_ids = write_concatenated_alignment(
        pair_ids,
        aln_1,
        aln_2)

    # save the alignment files
    mon_alignment_file_1 = f"{outdir}/{aln_1.main_id}_monomer_1.a3m"
    with open(mon_alignment_file_1, "w") as of:
        write_a3m(sequences_monomer_1, of)

    mon_alignment_file_2 = f"{outdir}/{aln_2.main_id}_monomer_2.a3m"
    with open(mon_alignment_file_2, "w") as of:
        write_a3m(sequences_monomer_2, of)

    pair_ids.to_csv(f"{outdir}/{aln_1.main_id}_{aln_2.main_id}_interact.csv", index=False)

    complex_ailgnment_file = f"{outdir}/{aln_1.main_id}_{aln_2.main_id}.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)

    return {'aln_file': complex_ailgnment_file, 'pair_ids': pair_ids}


def concatenate_alignments(inparams):
    runners, alignment, methods, hhfilter = inparams

    chain1_a3ms = alignment["chain1"]
    chain2_a3ms = alignment["chain2"]
    outdir = alignment['outdir']

    print(f"Concatenating {chain1_a3ms['name']} and {chain2_a3ms['name']}")

    try:

        makedir_if_not_exists(outdir)

        uniref_sto_aln1, uniref_sto_aln2 = None, None
        with open(chain1_a3ms["uniref90_sto"]) as f:
            uniref_sto_aln1 = Alignment.from_file(f, format="stockholm")
        with open(chain2_a3ms["uniref90_sto"]) as f:
            uniref_sto_aln2 = Alignment.from_file(f, format="stockholm")

        uniclust_a3m_aln1, uniclust_a3m_aln2 = None, None
        with open(chain1_a3ms["uniclust30_a3m"]) as f:
            uniclust_a3m_aln1 = Alignment.from_file(f, format="a3m")
        with open(chain2_a3ms["uniclust30_a3m"]) as f:
            uniclust_a3m_aln2 = Alignment.from_file(f, format="a3m")

        uniref_a3m_aln1, uniref_a3m_aln2 = None, None
        if chain1_a3ms["uniref30_a3m"] is not None:
            with open(chain1_a3ms["uniref30_a3m"]) as f:
                uniref_a3m_aln1 = Alignment.from_file(f, format="a3m")
        if chain2_a3ms["uniref30_a3m"] is not None:
            with open(chain2_a3ms["uniref30_a3m"]) as f:
                uniref_a3m_aln2 = Alignment.from_file(f, format="a3m")

        uniprot_sto_aln1, uniprot_sto_aln2 = None, None
        if chain1_a3ms["uniprot_sto"] is not None:
            with open(chain1_a3ms["uniprot_sto"]) as f:
                uniprot_sto_aln1 = Alignment.from_file(f, format="stockholm")
        if chain2_a3ms["uniprot_sto"] is not None:
            with open(chain2_a3ms["uniprot_sto"]) as f:
                uniprot_sto_aln2 = Alignment.from_file(f, format="stockholm")

        for method in methods:
            if method == "fused":
                sequences1 = []
                sequences2 = []

                sequences1 += uniref_sto_aln1.seqs
                sequences1 += uniclust_a3m_aln1.seqs
                with open(chain1_a3ms["mgnify_sto"]) as f:
                    ali = next(read_stockholm(f))
                    sequences1 += list(ali.seqs.values())

                with open(chain1_a3ms["bfd_a3m"]) as f:
                    seqs = read_a3m(f)
                    sequences1 += list(seqs.values())

                sequences2 += uniref_sto_aln2.seqs
                sequences2 += uniclust_a3m_aln2.seqs
                with open(chain2_a3ms["mgnify_sto"]) as f:
                    ali = next(read_stockholm(f))
                    sequences2 += list(ali.seqs.values())

                with open(chain2_a3ms["bfd_a3m"]) as f:
                    seqs = read_a3m(f)
                    sequences2 += list(seqs.values())

                fused_msa(sequences1, sequences2, outdir + '/fused.a3m')
                alignment['fused_msa'] = outdir + '/fused.a3m'

            elif method == "geno_dist":
                if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                    # print(uniref_a3m_aln1.headers)

                    pair_ids = runners['Geno_interact'].get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)
                    alignment["geno_dist_uniref_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1, uniref_a3m_aln2,
                                                                         hhfilter,
                                                                         outdir + '/geno_dist_uniref_a3m')
                    print(f"geno_dist_uniref_a3m: {len(pair_ids)} pairs")

                if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                    pair_ids = runners['Geno_interact'].get_interactions(uniref_sto_aln1, uniref_sto_aln2)
                    alignment["geno_dist_uniref_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1, uniref_sto_aln2,
                                                                         hhfilter,
                                                                         outdir + '/geno_dist_uniref_sto')
                    print(f"geno_dist_uniref_sto: {len(pair_ids)} pairs")

                if uniprot_sto_aln1 is not None and uniprot_sto_aln2 is not None:
                    pair_ids = runners['Geno_interact'].get_interactions(uniprot_sto_aln1, uniprot_sto_aln2)
                    alignment["geno_dist_uniprot_sto"] = write_dimer_a3ms(pair_ids, uniprot_sto_aln1,
                                                                          uniprot_sto_aln2, hhfilter,
                                                                          outdir + '/geno_dist_uniprot_sto')
                    print(f"geno_dist_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "pdb_interact":
                if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                    pair_ids = runners['pdb_interact'].get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)
                    alignment["pdb_interact_uniref_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1,
                                                                            uniref_a3m_aln2, hhfilter,
                                                                            outdir + '/pdb_interact_uniref_a3m')
                    print(f"pdb_interact_uniref_a3m: {len(pair_ids)} pairs")

                if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                    pair_ids = runners['pdb_interact'].get_interactions(uniref_sto_aln1, uniref_sto_aln2)
                    alignment["pdb_interact_uniref_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1,
                                                                            uniref_sto_aln2, hhfilter,
                                                                            outdir + '/pdb_interact_uniref_sto')
                    print(f"pdb_interact_uniref_sto: {len(pair_ids)} pairs")

                if uniprot_sto_aln1 is not None and uniprot_sto_aln2 is not None:
                    pair_ids = runners['pdb_interact'].get_interactions(uniprot_sto_aln1, uniprot_sto_aln2)
                    alignment["pdb_interact_uniprot_sto"] = write_dimer_a3ms(pair_ids, uniprot_sto_aln1,
                                                                             uniprot_sto_aln2, hhfilter,
                                                                             outdir + '/pdb_interact_uniprot_sto')
                    print(f"pdb_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "species_interact":
                if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                    pair_ids = Species_interact.get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)
                    alignment["species_interact_uniref_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1,
                                                                                uniref_a3m_aln2, hhfilter,
                                                                                outdir + '/species_interact_uniref_a3m')
                    print(f"species_interact_uniref_a3m: {len(pair_ids)} pairs")

                if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                    pair_ids = Species_interact.get_interactions(uniref_sto_aln1, uniref_sto_aln2)
                    alignment["species_interact_uniref_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1,
                                                                                uniref_sto_aln2, hhfilter,
                                                                                outdir + '/species_interact_uniref_sto')
                    print(f"species_interact_uniref_sto: {len(pair_ids)} pairs")

                if uniprot_sto_aln1 is not None and uniprot_sto_aln2 is not None:
                    pair_ids = Species_interact.get_interactions(uniprot_sto_aln1, uniprot_sto_aln2)
                    alignment["species_interact_uniprot_sto"] = write_dimer_a3ms(pair_ids, uniprot_sto_aln1,
                                                                                 uniprot_sto_aln2, hhfilter,
                                                                                 outdir + '/species_interact_uniprot_sto')
                    print(f"species_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "string_interact":
                if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                    pair_ids = runners['string_interact'].get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)
                    alignment["string_interact_uniref_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1,
                                                                               uniref_a3m_aln2, hhfilter,
                                                                               outdir + '/string_interact_uniref_a3m')
                    print(f"string_interact_uniref_a3m: {len(pair_ids)} pairs")

                if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                    pair_ids = runners['string_interact'].get_interactions(uniref_sto_aln1, uniref_sto_aln2)
                    alignment["string_interact_uniref_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1,
                                                                               uniref_sto_aln2, hhfilter,
                                                                               outdir + '/string_interact_uniref_sto')
                    print(f"string_interact_uniref_sto: {len(pair_ids)} pairs")

                if uniprot_sto_aln1 is not None and uniprot_sto_aln2 is not None:
                    pair_ids = runners['string_interact'].get_interactions(uniprot_sto_aln1, uniprot_sto_aln2)
                    alignment["string_interact_uniprot_sto"] = write_dimer_a3ms(pair_ids, uniprot_sto_aln1,
                                                                                uniprot_sto_aln2, hhfilter,
                                                                                outdir + '/string_interact_uniprot_sto')
                    print(f"string_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "uniclust_oxmatch":
                if uniclust_a3m_aln1 is not None and uniclust_a3m_aln2 is not None:
                    pair_ids = UNICLUST_oxmatch.get_interactions(uniclust_a3m_aln1, uniclust_a3m_aln2)
                    alignment["uniclust_oxmatch_a3m"] = write_dimer_a3ms(pair_ids, uniclust_a3m_aln1,
                                                                         uniclust_a3m_aln2, hhfilter,
                                                                         outdir + '/uniclust_oxmatch_a3m')
                    print(f"uniclust_oxmatch_a3m: {len(pair_ids)} pairs")

            elif method == "uniprot_distance":
                if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                    pair_ids = UNIPROT_distance.get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)
                    alignment["uniprot_distance_uniref_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1,
                                                                                uniref_a3m_aln2, hhfilter,
                                                                                outdir + '/uniprot_distance_uniref_a3m')
                    print(f"uniprot_distance_uniref_a3m: {len(pair_ids)} pairs")

                if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                    pair_ids = UNIPROT_distance.get_interactions(uniref_sto_aln1, uniref_sto_aln2)
                    alignment["uniprot_distance_uniref_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1,
                                                                                uniref_sto_aln2, hhfilter,
                                                                                outdir + '/uniprot_distance_uniref_sto')
                    print(f"uniprot_distance_uniref_sto: {len(pair_ids)} pairs")

                if uniprot_sto_aln1 is not None and uniprot_sto_aln2 is not None:
                    pair_ids = UNIPROT_distance.get_interactions(uniprot_sto_aln1, uniprot_sto_aln2)
                    alignment["uniprot_distance_uniprot_sto"] = write_dimer_a3ms(pair_ids, uniprot_sto_aln1,
                                                                                 uniprot_sto_aln2, hhfilter,
                                                                                 outdir + '/uniprot_distance_uniprot_sto')
                    print(f"uniprot_distance_uniprot_sto: {len(pair_ids)} pairs")

    except Exception as e:
        print(e)
        return None


class Complex_alignment_concatenation_pipeline:

    def __init__(self, params, multiprocess=False, process_num=1):

        self.params = params

        self.methods = params['concatenate_methods'].split(',')

        print("Using methods:")
        print(self.methods)

        self.runners = {}

        if "geno_dist" in self.methods:
            self.Geno_interact_runner = Geno_interact(self.params['uniprot_to_embl_table'],
                                                      self.params['ena_genome_location_table'],
                                                      self.params['genome_distance_threshold'])
            self.runners['Geno_interact'] = self.Geno_interact_runner

        if "pdb_interact" in self.methods:
            self.pdb_interact_runner = PDB_interact(self.params['uniprot2pdb_mapping_file'], self.params['dimers_list'])
            self.pdb_interact_runner.load_data()
            self.runners['pdb_interact'] = self.pdb_interact_runner

        if "string_interact" in self.methods:
            self.string_interact_runner = STRING_interact(self.params['string2uniprot_map'])
            self.string_interact_runner.load_data(750)
            self.runners['string_interact'] = self.string_interact_runner

        self.multiprocess = multiprocess
        self.process_num = process_num

    def concatenate(self, alignments, hhfilter):
        res_alignments = []
        if self.multiprocess:
            concatenate_list = []
            for alignment in alignments:
                concatenate_list.append([self.runners, alignment, self.methods, hhfilter])
            pool = Pool(processes=self.process_num)
            res_alignments = pool.map(concatenate_alignments, concatenate_list)
            pool.close()
            pool.join()
        else:
            for alignment in alignments:
                res_alignments += [concatenate_alignments([self.runners, alignment, self.methods, hhfilter])]
        return res_alignments
