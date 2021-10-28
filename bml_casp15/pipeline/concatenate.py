import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_dir, is_file, read_option_file, makedir_if_not_exists
from bml_casp15.complex.geno_dist import Geno_interact
from bml_casp15.complex.pdb_interact import PDB_interact
from bml_casp15.complex.species_interact import Species_interact
from bml_casp15.complex.string_interact import STRING_interact
from bml_casp15.complex.uniclust_oxmatch import UNICLUST_oxmatch
from bml_casp15.complex.uniprot_distance import UNIPROT_distance
from bml_casp15.alignment.alignment import *
from bml_casp15.complex.complex import write_concatenated_alignment
from bml_casp15.common.util import makedir_if_not_exists


def write_dimer_a3ms(pair_ids, aln_1, aln_2, outdir):

    makedir_if_not_exists(outdir)

    target_header, sequences_full, sequences_monomer_1, sequences_monomer_2 = write_concatenated_alignment(pair_ids,
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
    print(pair_ids)

    complex_ailgnment_file = f"{outdir}/{aln_1.main_id}_{aln_2.main_id}.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)

    return complex_ailgnment_file


class Concatenate:

    def __init__(self, params, concatenate_methods=["geno_dist", "pdb_interact", "species_interact",
                                                    "string_interact", "uniclust_oxmatch", "uniprot_distance"]):
        self.params = params
        self.methods = concatenate_methods

    def concatenate(self, alignments):

        Geno_interact_runner, pdb_interact_runner, string_interact_runner = None, None, None

        for alignment in alignments:

            chain1_a3ms = alignment["chain1"]

            chain2_a3ms = alignment["chain2"]

            print(f"Concatenating {chain1_a3ms['name']} and {chain2_a3ms['name']}")
            uniref_a3m_aln1, uniref_a3m_aln2 = None, None
            if chain1_a3ms["uniref_a3m"] is not None:
                with open(chain1_a3ms["uniref_a3m"]) as f:
                    uniref_a3m_aln1 = Alignment.from_file(f, format="a3m")

            if chain2_a3ms["uniref_a3m"] is not None:
                with open(chain2_a3ms["uniref_a3m"]) as f:
                    uniref_a3m_aln2 = Alignment.from_file(f, format="a3m")

            uniref_sto_aln1, uniref_sto_aln2 = None, None
            if chain1_a3ms["uniref_sto"] is not None:
                with open(chain1_a3ms["uniref_sto"]) as f:
                    uniref_sto_aln1 = Alignment.from_file(f, format="stockholm")

            if chain2_a3ms["uniref_sto"] is not None:
                with open(chain2_a3ms["uniref_sto"]) as f:
                    uniref_sto_aln2 = Alignment.from_file(f, format="stockholm")

            uniclust_a3m_aln1, uniclust_a3m_aln2 = None, None
            if chain1_a3ms["uniclust_a3m"] is not None:
                with open(chain1_a3ms["uniclust_a3m"]) as f:
                    uniclust_a3m_aln1 = Alignment.from_file(f, format="a3m")

            if chain2_a3ms["uniclust_a3m"] is not None:
                with open(chain2_a3ms["uniclust_a3m"]) as f:
                    uniclust_a3m_aln2 = Alignment.from_file(f, format="a3m")

            for method in self.methods:

                if method == "geno_dist":

                    if Geno_interact_runner is None:
                        Geno_interact_runner = Geno_interact(self.params['uniprot_to_embl_table'],
                                                             self.params['ena_genome_location_table'])

                    if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:

                        pair_ids = Geno_interact_runner.get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)

                        alignment["geno_dist_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1, uniref_a3m_aln2,
                                                         alignment['outdir'] + '/geno_dist_a3m')

                    if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                        pair_ids = Geno_interact_runner.get_interactions(uniref_sto_aln1, uniref_sto_aln2)

                        alignment["geno_dist_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1, uniref_sto_aln2,
                                                         alignment['outdir'] + '/geno_dist_sto')

                elif method == "pdb_interact":

                    if pdb_interact_runner is None:
                        pdb_interact_runner = PDB_interact(self.params['uniprot2pdb_mapping_file'],
                                                                        self.params['dimers_list'])
                        pdb_interact_runner.load_data()

                    if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                        pair_ids = pdb_interact_runner.get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)

                        alignment["pdb_interact_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1, uniref_a3m_aln2,
                                                         alignment['outdir'] + '/pdb_interact_a3m')

                    if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                        pair_ids = pdb_interact_runner.get_interactions(uniref_sto_aln1, uniref_sto_aln2)

                        alignment["pdb_interact_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1, uniref_sto_aln2,
                                                         alignment['outdir'] + '/pdb_interact_sto')

                elif method == "species_interact":

                    if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                        pair_ids = Species_interact.get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)

                        alignment["species_interact_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1, uniref_a3m_aln2,
                                                         alignment['outdir'] + '/species_interact_a3m')

                    if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                        pair_ids = Species_interact.get_interactions(uniref_sto_aln1, uniref_sto_aln2)

                        alignment["species_interact_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1, uniref_sto_aln2,
                                                         alignment['outdir'] + '/species_interact_sto')

                elif method == "string_interact":

                    if string_interact_runner is None:
                        string_interact_runner = STRING_interact(self.params['string2uniprot_map'])
                        string_interact_runner.load_data()

                    if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                        pair_ids = string_interact_runner.get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)

                        alignment["string_interact_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1, uniref_a3m_aln2,
                                                         alignment['outdir'] + '/string_interact_a3m')

                    if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                        pair_ids = string_interact_runner.get_interactions(uniref_sto_aln1, uniref_sto_aln2)

                        alignment["string_interact_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1, uniref_sto_aln2,
                                                         alignment['outdir'] + '/string_interact_sto')

                elif method == "uniclust_oxmatch":

                    if uniclust_a3m_aln1 is not None and uniclust_a3m_aln2 is not None:
                        pair_ids = UNICLUST_oxmatch.get_interactions(uniclust_a3m_aln1, uniclust_a3m_aln2)

                        alignment["uniclust_oxmatch_a3m"] = write_dimer_a3ms(pair_ids, uniclust_a3m_aln1, uniclust_a3m_aln2,
                                                         alignment['outdir'] + '/uniclust_oxmatch_a3m')

                elif method == "uniprot_distance":

                    if uniref_a3m_aln1 is not None and uniref_a3m_aln2 is not None:
                        pair_ids = UNIPROT_distance.get_interactions(uniref_a3m_aln1, uniref_a3m_aln2)

                        alignment["uniprot_distance_a3m"] = write_dimer_a3ms(pair_ids, uniref_a3m_aln1, uniref_a3m_aln2,
                                                         alignment['outdir'] + '/uniprot_distance_a3m')

                    if uniref_sto_aln1 is not None and uniref_sto_aln2 is not None:
                        pair_ids = UNIPROT_distance.get_interactions(uniref_sto_aln1, uniref_sto_aln2)

                        alignment["uniprot_distance_sto"] = write_dimer_a3ms(pair_ids, uniref_sto_aln1, uniref_sto_aln2,
                                                         alignment['outdir'] + '/uniprot_distance_sto')

        return alignments
