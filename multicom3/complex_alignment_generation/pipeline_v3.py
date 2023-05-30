import argparse
import os
import sys
import time
from multiprocessing import Pool

from multicom3.common.util import is_dir, is_file, read_option_file, makedir_if_not_exists
from multicom3.common.util import makedir_if_not_exists
from multicom3.complex_alignment_generation.pdb_interact_v3 import PDB_interact_v3
from multicom3.complex_alignment_generation.species_interact_v3 import Species_interact_v3
from multicom3.complex_alignment_generation.string_interact_v3 import STRING_interact_v3
from multicom3.complex_alignment_generation.uniclust_oxmatch_v3 import UNICLUST_oxmatch_v3
from multicom3.complex_alignment_generation.uniprot_distance_v3 import UNIPROT_distance_v3
from multicom3.monomer_alignment_generation.alignment import *
from tqdm import tqdm

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


def write_concatenated_alignment(paired_rows, alignments):
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
        seen_seqs[alignments[i].main_id] = [alignments[i].main_seq]

    pair_id_count = 0
    for pair_index in range(1, len(list(paired_rows[:, 0]))):

        row_indices = list(paired_rows[pair_index, :])
        seqs = []
        headers = []

        for j in range(len(alignments)):
            index = row_indices[j]
            header = ""
            if index == -1:
                header = f'placeholder{pair_id_count}'
                seq = '-' * len(alignments[j].main_seq)
            else:
                seq = alignments[j].seqs[index]
                header, start, end = parse_header(alignments[j].headers[index])
                header = f"{header}_{start}-{end}"
            
            headers += [header]
            seqs += [seq]            
            filter_pair_ids[f'id_{j + 1}'] += [header]
            filter_pair_ids[f'index_{j + 1}'] += [pair_id_count + 1]

        concatenated_header = _prepare_header(headers)

        # save the information
        sequences_to_write[concatenated_header] = seqs

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


def write_multimer_a3ms(pair_ids, alignments, outdir, method, is_homomers=False):
    outdir = outdir + '/' + method

    makedir_if_not_exists(outdir)

    sequences_full, sequences_monomers, pair_ids = write_concatenated_alignment(pair_ids, alignments)

    pair_ids.to_csv(f"{outdir}/{method}_interact.csv", index=False)

    complex_alignment_file = f"{outdir}/{method}.a3m"
    with open(complex_alignment_file, "w") as of:
        write_a3m(sequences_full, of)

    # save the alignment files
    monomer_alignment_files = []
    for monomer_id in sequences_monomers:
        mon_alignment_file = f"{outdir}/{monomer_id}_con.a3m"
        with open(mon_alignment_file, "w") as of:
            write_a3m(sequences_monomers[monomer_id], of)
        monomer_alignment_files += [mon_alignment_file]

    if is_homomers:
        max_seq_num = 50000
        per_max_seq_num = int(max_seq_num - len(pair_ids)) / len(sequences_monomers)

        unpaired_sequences = {}
        homomers_sequences = {}
        for alignment, monomer_id in zip(alignments, sequences_monomers):
            homomers_sequences[monomer_id] = {'headers': [], 'seqs': []}

            unpaired_sequences[monomer_id] = {}
            paired_headers = [header for header in sequences_monomers[monomer_id]]
            paired_sequences = [sequences_monomers[monomer_id][header] for header in sequences_monomers[monomer_id]]
            for seqindx, seq in enumerate(alignment.seqs):
                if seq not in paired_sequences and alignment.ids[seqindx] not in paired_headers:
                    unpaired_sequences[monomer_id][alignment.ids[seqindx]] = seq
                    if len(unpaired_sequences[monomer_id]) >= per_max_seq_num:
                        break

        seqlen = len(alignments[0].main_seq)
        for monomer_id in unpaired_sequences:
            add_count = 0
            for header in unpaired_sequences[monomer_id]:
                homomers_sequences[monomer_id]['headers'] += [header]
                homomers_sequences[monomer_id]['seqs'] += [unpaired_sequences[monomer_id][header]]
                add_count += 1
            for other_monomer_id in homomers_sequences:
                if other_monomer_id == monomer_id:
                    continue
                for i in range(add_count):
                    homomers_sequences[other_monomer_id]['headers'] += ['placeholder']
                    homomers_sequences[other_monomer_id]['seqs'] += ['-' * seqlen]
        # print(homomers_sequences)

        for monomer_id in homomers_sequences:
            mon_alignment_file = f"{outdir}/{monomer_id}_con.a3m"
            with open(mon_alignment_file, "a") as of:
                for i in range(len(homomers_sequences[monomer_id]['headers'])):
                    of.write(f">{homomers_sequences[monomer_id]['headers'][i]}\n"
                             f"{homomers_sequences[monomer_id]['seqs'][i]}\n")

    return {'aln_file': complex_alignment_file, 'pair_ids': pair_ids, 'monomer_files': monomer_alignment_files}


def concatenate_alignments(inparams):
    runners, alignment, methods, hhfilter, is_homomers = inparams

    outdir = alignment['outdir']
    print(f"Concatenating {' and '.join([alignment[chain]['name'] for chain in alignment if chain != 'outdir'])}")
    try:

        makedir_if_not_exists(outdir)

        uniref_sto_alignments = []
        uniclust_a3m_alignments = []
        uniref_a3m_alignments = []
        uniprot_sto_alignments = []
        colabfold_alignments = []

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
            colabfold_alignments += [alignment[chain]["colabfold_a3m"]]

        for method in methods:
            if method == "pdb_interact":
                if len(uniref_a3m_alignments) > 0:
                    pair_ids = runners['pdb_interact'].get_interactions_v2(uniref_a3m_alignments, is_homomers)
                    alignment["pdb_interact_uniref_a3m"] = write_multimer_a3ms(pair_ids, uniref_a3m_alignments,
                                                                               outdir, 'pdb_interact_uniref_a3m',
                                                                               is_homomers)
                    print(f"pdb_interact_uniref_a3m: {len(pair_ids)} pairs")

                if len(uniref_sto_alignments) > 0:
                    pair_ids = runners['pdb_interact'].get_interactions_v2(uniref_sto_alignments, is_homomers)
                    alignment["pdb_interact_uniref_sto"] = write_multimer_a3ms(pair_ids, uniref_sto_alignments,
                                                                               outdir, 'pdb_interact_uniref_sto',
                                                                               is_homomers)
                    print(f"pdb_interact_uniref_sto: {len(pair_ids)} pairs")

                if len(uniprot_sto_alignments) > 0:
                    pair_ids = runners['pdb_interact'].get_interactions_v2(uniprot_sto_alignments, is_homomers)
                    alignment["pdb_interact_uniprot_sto"] = write_multimer_a3ms(pair_ids, uniprot_sto_alignments,
                                                                                outdir, 'pdb_interact_uniprot_sto',
                                                                                is_homomers)
                    print(f"pdb_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "species_interact":
                if len(uniref_a3m_alignments) > 0:
                    pair_ids = Species_interact_v3.get_interactions_v2(uniref_a3m_alignments)
                    alignment["species_interact_uniref_a3m"] = write_multimer_a3ms(pair_ids, uniref_a3m_alignments,
                                                                                   outdir,
                                                                                   'species_interact_uniref_a3m',
                                                                                   is_homomers)
                    print(f"species_interact_uniref_a3m: {len(pair_ids)} pairs")

                if len(uniref_sto_alignments) > 0:
                    pair_ids = Species_interact_v3.get_interactions_v2(uniref_sto_alignments)
                    alignment["species_interact_uniref_sto"] = write_multimer_a3ms(pair_ids, uniref_sto_alignments,
                                                                                   outdir,
                                                                                   'species_interact_uniref_sto',
                                                                                   is_homomers)
                    print(f"species_interact_uniref_sto: {len(pair_ids)} pairs")

                if len(uniprot_sto_alignments) > 0:
                    pair_ids = Species_interact_v3.get_interactions_v2(uniprot_sto_alignments)
                    alignment["species_interact_uniprot_sto"] = write_multimer_a3ms(pair_ids, uniprot_sto_alignments,
                                                                                    outdir,
                                                                                    'species_interact_uniprot_sto',
                                                                                    is_homomers)
                    print(f"species_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "string_interact":
                if len(uniref_a3m_alignments) > 0:
                    pair_ids = runners['string_interact'].get_interactions_v2(uniref_a3m_alignments)
                    alignment["string_interact_uniref_a3m"] = write_multimer_a3ms(pair_ids, uniref_a3m_alignments,
                                                                                  outdir, 'string_interact_uniref_a3m',
                                                                                  is_homomers)
                    print(f"string_interact_uniref_a3m: {len(pair_ids)} pairs")

                if len(uniref_sto_alignments) > 0:
                    pair_ids = runners['string_interact'].get_interactions_v2(uniref_sto_alignments)
                    alignment["string_interact_uniref_sto"] = write_multimer_a3ms(pair_ids, uniref_sto_alignments,
                                                                                  outdir, 'string_interact_uniref_sto',
                                                                                  is_homomers)
                    print(f"string_interact_uniref_sto: {len(pair_ids)} pairs")

                if len(uniprot_sto_alignments) > 0:
                    pair_ids = runners['string_interact'].get_interactions_v2(uniprot_sto_alignments)
                    alignment["string_interact_uniprot_sto"] = write_multimer_a3ms(pair_ids, uniprot_sto_alignments,
                                                                                   outdir,
                                                                                   'string_interact_uniprot_sto',
                                                                                   is_homomers)
                    print(f"string_interact_uniprot_sto: {len(pair_ids)} pairs")

            elif method == "uniclust_oxmatch":
                if len(uniclust_a3m_alignments) > 0:
                    pair_ids = UNICLUST_oxmatch_v3.get_interactions_v2(uniclust_a3m_alignments)
                    alignment["uniclust_oxmatch_a3m"] = write_multimer_a3ms(pair_ids, uniclust_a3m_alignments,
                                                                            outdir, 'uniclust_oxmatch_a3m', is_homomers)
                    print(f"uniclust_oxmatch_a3m: {len(pair_ids)} pairs")

            elif method == "uniprot_distance":
                if len(uniref_a3m_alignments) > 0:
                    pair_ids = UNIPROT_distance_v3.get_interactions_v2(uniref_a3m_alignments)
                    alignment["uniprot_distance_uniref_a3m"] = write_multimer_a3ms(pair_ids, uniref_a3m_alignments,
                                                                                   outdir,
                                                                                   'uniprot_distance_uniref_a3m',
                                                                                   is_homomers)
                    print(f"uniprot_distance_uniref_a3m: {len(pair_ids)} pairs")

                if len(uniref_sto_alignments) > 0:
                    pair_ids = UNIPROT_distance_v3.get_interactions_v2(uniref_sto_alignments)
                    alignment["uniprot_distance_uniref_sto"] = write_multimer_a3ms(pair_ids, uniref_sto_alignments,
                                                                                   outdir,
                                                                                   'uniprot_distance_uniref_sto',
                                                                                   is_homomers)
                    print(f"uniprot_distance_uniref_sto: {len(pair_ids)} pairs")

                if len(uniprot_sto_alignments) > 0:
                    pair_ids = UNIPROT_distance_v3.get_interactions_v2(uniprot_sto_alignments)
                    alignment["uniprot_distance_uniprot_sto"] = write_multimer_a3ms(pair_ids, uniprot_sto_alignments,
                                                                                    outdir,
                                                                                    'uniprot_distance_uniprot_sto',
                                                                                    is_homomers)
                    print(f"uniprot_distance_uniprot_sto: {len(pair_ids)} pairs")

        os.system(f"touch {outdir}/DONE")

    except Exception as e:
        print(e)
        return None


class Complex_alignment_concatenation_pipeline:

    def __init__(self, params, run_methods, multiprocess=False, process_num=1):

        self.params = params

        self.methods = run_methods

        print("Using methods:")
        print(self.methods)

        self.runners = {}

        if "pdb_interact" in self.methods:
            self.pdb_interact_runner = PDB_interact_v3(self.params['uniprot2pdb_mapping_file'],
                                                       self.params['complexes_list'])
            self.pdb_interact_runner.load_data()
            self.runners['pdb_interact'] = self.pdb_interact_runner

        if "string_interact" in self.methods:
            self.string_interact_runner = STRING_interact_v3(self.params['string2uniprot_map'])
            self.string_interact_runner.load_data(750)
            self.runners['string_interact'] = self.string_interact_runner

        self.multiprocess = multiprocess
        self.process_num = process_num

    def concatenate(self, alignments, hhfilter, is_homomers):
        res_alignments = []
        if self.multiprocess:
            concatenate_list = []
            for alignment in alignments:
                concatenate_list.append([self.runners, alignment, self.methods, hhfilter, is_homomers])
            pool = Pool(processes=self.process_num)
            res_alignments = pool.map(concatenate_alignments, concatenate_list)
            pool.close()
            pool.join()
        else:
            for alignment in alignments:
                res_alignments += [
                    concatenate_alignments([self.runners, alignment, self.methods, hhfilter, is_homomers])]
        return res_alignments
