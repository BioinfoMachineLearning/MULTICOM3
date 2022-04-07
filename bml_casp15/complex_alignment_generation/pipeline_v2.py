import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_dir, is_file, read_option_file, makedir_if_not_exists
from bml_casp15.complex_alignment_generation.geno_dist import Geno_interact
from bml_casp15.complex_alignment_generation.pdb_interact_v2 import PDB_interact_v2
from bml_casp15.complex_alignment_generation.species_interact_v2 import Species_interact_v2
from bml_casp15.complex_alignment_generation.string_interact_v2 import STRING_interact_v2
from bml_casp15.complex_alignment_generation.uniclust_oxmatch_v2 import UNICLUST_oxmatch_v2
from bml_casp15.complex_alignment_generation.uniprot_distance_v2 import UNIPROT_distance_v2
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


def write_concatenated_alignment(id_pairing, alignments):
    def _prepare_header(ids):
        return "_____".join(ids)

    sequences_to_write = {}  # list of (header,seq1,seq2) tuples

    target_header = _prepare_header([alignment.main_id for alignment in alignments])

    sequences_to_write[target_header] = [alignment.main_seq for alignment in alignments]
    # print(sequences_to_write)

    # create other headers and sequences
    seen_seqs = {'full': []}
    filter_pair_ids = {}
    for i in range(len(alignments)):
        filter_pair_ids[f"id_{i + 1}"] = [alignments[i].main_id]
        filter_pair_ids[f"index_{i + 1}"] = [0]
        seen_seqs[alignments[i].main_id] = []

    pair_id_count = 0
    for i in range(len(id_pairing)):
        seqs = []
        headers = []

        seen_monomer_seq = False
        for j in range(len(alignments)):
            id = id_pairing.loc[i, f"id_{j + 1}"]
            seq = alignments[j][id]
            if seq in seen_seqs[alignments[j].main_id]:
                seen_monomer_seq = True
                break

        if seen_monomer_seq:
            continue

        for j in range(len(alignments)):
            id = id_pairing.loc[i, f"id_{j + 1}"]
            seq = alignments[j][id]
            seen_seqs[alignments[j].main_id] += [seq]
            seqs += [seq]
            header, _, _ = parse_header(alignments[j].headers[id][0])
            headers += [header]

        combine_seq = "".join(seqs)
        if combine_seq in seen_seqs['full']:
            continue

        seen_seqs['full'] += [combine_seq]

        concatenated_header = _prepare_header(headers)

        # save the information
        sequences_to_write[concatenated_header] = seqs

        for j in range(len(alignments)):
            filter_pair_ids[f'id_{j + 1}'] += [id_pairing.loc[i, f"id_{j + 1}"]]
            filter_pair_ids[f'index_{j + 1}'] += [pair_id_count + 1]

        pair_id_count += 1

    sequences_full = OrderedDict()
    sequences_monomers = {}
    for i in range(len(alignments)):
        sequences_monomers[alignments[i].main_id] = OrderedDict()

    for header_full in sequences_to_write:
        seq_full = "".join(sequences_to_write[header_full])
        sequences_full[header_full] = seq_full

        headers = header_full.split('_____')
        for i in range(len(alignments)):
            sequences_monomers[alignments[i].main_id][headers[i]] = sequences_to_write[header_full][i]

    return sequences_full, sequences_monomers, pd.DataFrame(filter_pair_ids)


def write_multimer_a3ms(pair_ids, alignments, outdir, method):
    outdir = outdir + '/' + method

    makedir_if_not_exists(outdir)

    sequences_full, sequences_monomers, pair_ids = write_concatenated_alignment(pair_ids, alignments)

    pair_ids.to_csv(f"{outdir}/{method}_interact.csv", index=False)

    complex_alignment_file = f"{outdir}/{method}.a3m"
    with open(complex_alignment_file, "w") as of:
        write_a3m(sequences_full, of)

    # save the alignment files
    for monomer_id in sequences_monomers:
        mon_alignment_file = f"{outdir}/{monomer_id}_con.a3m"
        with open(mon_alignment_file, "w") as of:
            write_a3m(sequences_monomers[monomer_id], of)

    return {'aln_file': complex_alignment_file, 'pair_ids': pair_ids}


def concatenate_alignments(inparams):
    runners, alignment, methods, hhfilter = inparams

    outdir = alignment['outdir']
    print(f"Concatenating {' and '.join([alignment[chain]['name'] for chain in alignment if chain != 'outdir'])}")
    try:

        makedir_if_not_exists(outdir)

        uniref_sto_alignments = []
        uniclust_a3m_alignments = []
        uniref_a3m_alignments = []
        uniprot_sto_alignments = []

        for chain in alignment:
            if chain == "outdir":
                continue
            with open(alignment[chain]["uniref90_sto"]) as f:
                uniref_sto_alignments += [Alignment.from_file(f, format="stockholm")]
            with open(alignment[chain]["uniclust30_a3m"]) as f:
                uniclust_a3m_alignments += [Alignment.from_file(f, format="a3m")]
            with open(alignment[chain]["uniref30_a3m"]) as f:
                uniref_a3m_alignments += [Alignment.from_file(f, format="a3m")]
            with open(alignment[chain]["uniprot_sto"]) as f:
                uniprot_sto_alignments += [Alignment.from_file(f, format="stockholm")]

        for method in methods:
            if method == "pdb_interact":
                if len(uniref_a3m_alignments) > 0:
                    pair_ids = runners['pdb_interact'].get_interactions(uniref_a3m_alignments)
                    alignment["pdb_interact_uniref_a3m"] = write_multimer_a3ms(pair_ids, uniref_a3m_alignments,
                                                                               outdir, 'pdb_interact_uniref_a3m')
                    print(f"pdb_interact_uniref_a3m: {len(pair_ids)} pairs")

                if len(uniref_sto_alignments) > 0:
                    pair_ids = runners['pdb_interact'].get_interactions(uniref_sto_alignments)
                    alignment["pdb_interact_uniref_sto"] = write_multimer_a3ms(pair_ids, uniref_sto_alignments,
                                                                               outdir, 'pdb_interact_uniref_sto')
                    print(f"pdb_interact_uniref_sto: {len(pair_ids)} pairs")

                if len(uniprot_sto_alignments) > 0:
                    pair_ids = runners['pdb_interact'].get_interactions(uniprot_sto_alignments)
                    alignment["pdb_interact_uniprot_sto"] = write_multimer_a3ms(pair_ids, uniprot_sto_alignments,
                                                                                outdir, 'pdb_interact_uniprot_sto')
                    print(f"pdb_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "species_interact":
                if len(uniref_a3m_alignments) > 0:
                    pair_ids = Species_interact_v2.get_interactions(uniref_a3m_alignments)
                    alignment["species_interact_uniref_a3m"] = write_multimer_a3ms(pair_ids, uniref_a3m_alignments,
                                                                                   outdir, 'species_interact_uniref_a3m')
                    print(f"species_interact_uniref_a3m: {len(pair_ids)} pairs")

                if len(uniref_sto_alignments) > 0:
                    pair_ids = Species_interact_v2.get_interactions(uniref_sto_alignments)
                    alignment["species_interact_uniref_sto"] = write_multimer_a3ms(pair_ids, uniref_sto_alignments,
                                                                                   outdir, 'species_interact_uniref_sto')
                    print(f"species_interact_uniref_sto: {len(pair_ids)} pairs")

                if len(uniprot_sto_alignments) > 0:
                    pair_ids = Species_interact_v2.get_interactions(uniprot_sto_alignments)
                    alignment["species_interact_uniprot_sto"] = write_multimer_a3ms(pair_ids, uniprot_sto_alignments,
                                                                                    outdir, 'species_interact_uniprot_sto')
                    print(f"species_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "string_interact":
                if len(uniref_a3m_alignments) > 0:
                    pair_ids = runners['string_interact'].get_interactions(uniref_a3m_alignments)
                    alignment["string_interact_uniref_a3m"] = write_multimer_a3ms(pair_ids, uniref_a3m_alignments,
                                                                                  outdir, 'string_interact_uniref_a3m')
                    print(f"string_interact_uniref_a3m: {len(pair_ids)} pairs")

                if len(uniref_sto_alignments) > 0:
                    pair_ids = runners['string_interact'].get_interactions(uniref_sto_alignments)
                    alignment["string_interact_uniref_sto"] = write_multimer_a3ms(pair_ids, uniref_sto_alignments,
                                                                                  outdir, 'string_interact_uniref_sto')
                    print(f"string_interact_uniref_sto: {len(pair_ids)} pairs")

                if len(uniprot_sto_alignments) > 0:
                    pair_ids = runners['string_interact'].get_interactions(uniprot_sto_alignments)
                    alignment["string_interact_uniprot_sto"] = write_multimer_a3ms(pair_ids, uniprot_sto_alignments,
                                                                                   outdir, 'string_interact_uniprot_sto')
                    print(f"string_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "uniclust_oxmatch":
                if len(uniclust_a3m_alignments) > 0:
                    pair_ids = UNICLUST_oxmatch_v2.get_interactions(uniclust_a3m_alignments)
                    alignment["uniclust_oxmatch_a3m"] = write_multimer_a3ms(pair_ids, uniclust_a3m_alignments,
                                                                            outdir, 'uniclust_oxmatch_a3m')
                    print(f"uniclust_oxmatch_a3m: {len(pair_ids)} pairs")

            elif method == "uniprot_distance":
                if len(uniref_a3m_alignments) > 0:
                    pair_ids = UNIPROT_distance_v2.get_interactions(uniref_a3m_alignments)
                    alignment["uniprot_distance_uniref_a3m"] = write_multimer_a3ms(pair_ids, uniref_a3m_alignments,
                                                                                outdir, 'uniprot_distance_uniref_a3m')
                    print(f"uniprot_distance_uniref_a3m: {len(pair_ids)} pairs")

                if len(uniref_sto_alignments) > 0:
                    pair_ids = UNIPROT_distance_v2.get_interactions(uniref_sto_alignments)
                    alignment["uniprot_distance_uniref_sto"] = write_multimer_a3ms(pair_ids, uniref_sto_alignments,
                                                                                outdir, 'uniprot_distance_uniref_sto')
                    print(f"uniprot_distance_uniref_sto: {len(pair_ids)} pairs")

                if len(uniprot_sto_alignments) > 0:
                    pair_ids = UNIPROT_distance_v2.get_interactions(uniprot_sto_alignments)
                    alignment["uniprot_distance_uniprot_sto"] = write_multimer_a3ms(pair_ids, uniprot_sto_alignments,
                                                                                 outdir, 'uniprot_distance_uniprot_sto')
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
            self.pdb_interact_runner = PDB_interact_v2(self.params['uniprot2pdb_mapping_file'])
            self.pdb_interact_runner.load_data()
            self.runners['pdb_interact'] = self.pdb_interact_runner

        if "string_interact" in self.methods:
            self.string_interact_runner = STRING_interact_v2(self.params['string2uniprot_map'])
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
