import re
from collections import OrderedDict, defaultdict
from copy import copy
import numpy as np
import pandas as pd
from bml_casp15.alignment import alignment


class Geno_interact:

    def __init__(self, uniprot_to_embl_table, ena_genome_location_table):
        self.uniprot_to_embl_table = uniprot_to_embl_table
        self.ena_genome_location_table = ena_genome_location_table

    def extract_cds_ids(self, alignment):

        def _split_annotation_string(annotation_string):
            """
            reformats the ENA annotation string as a list of
            [(read,cds)] tuples
            """
            full_annotation = [
                tuple(x.split(":")) for x in
                annotation_string.split(",")
            ]  # list of lists in format [read,cds]

            return full_annotation

        def _remove_redundant_cds(uniprot_and_genome_cds):
            """
            Removes CDSs that have hits to multiple genomes
            Returns a list of tuples (Uniprot_AC, CDS)
            """

            filtered_uniprot_and_cds = []
            for uniprot_ac, genome_and_cds in uniprot_and_genome_cds:

                count_reads = defaultdict(list)

                for genome, cds in genome_and_cds:
                    count_reads[cds].append(genome)

                # check how many genomes are associated with a particular CDS,
                # only keep CDSs that can be matched to one genome
                for cds, genomes in count_reads.items():
                    if len(genomes) == 1:
                        filtered_uniprot_and_cds.append((uniprot_ac, cds))

            return filtered_uniprot_and_cds

        # store IDs in set for faster membership checking
        target_ids = set(alignment.ids)

        # initialize list of list of (uniprot,[(genome,CDS)]) tuples
        # eg [(uniprot1,[(genome1,cds1),(genome2,cds2)])]
        # TODO: I know this is ugly but I need to save the uniprot id
        # for later use, and I need to save the genome so that I can remove
        # cds's that are hit to multiple genomes
        genome_and_cds = []

        # read through the data table line-by-line to improve speed
        with open(self.uniprot_to_embl_table) as f:
            for line in f:
                # ena_data is formatted as 'genome1:cds1,genome2:cds2'
                uniprot_ac, _, ena_data = line.rstrip().split(" ")

                if uniprot_ac in target_ids:
                    genome_and_cds.append((
                        uniprot_ac, _split_annotation_string(ena_data)
                    ))

        # clean the uniprot to embl hits
        # returns a list of uniprot_ac, cds_id tuples
        filtered_uniprot_and_cds = _remove_redundant_cds(genome_and_cds)

        return filtered_uniprot_and_cds

    def extract_embl_annotation(self, uniprot_and_cds):

        # initialize list of list
        # will contain [[CDS,genome_id,uniprot_ac,gene_start,gene_end]]
        embl_cds_annotation = []

        # convert cds_target_ids to set to speed up membership checking
        cds_target_ids = [x for _, x in uniprot_and_cds]
        cds_target_set = set(cds_target_ids)

        # dictionary of uniprot acs to save for each cds
        cds_to_uniprot = {y: x for x, y in uniprot_and_cds}

        # extract the annotation
        # note: we don't use the Uniprot AC from this file
        # as the mapping is sometimes ambiguous
        # ie, a CDS can map to multiple Uniprot ACs
        with open(self.ena_genome_location_table) as inf:
            for line in inf:

                cds_id, genome_id, _, start, end = (
                    line.rstrip().split("\t")
                )

                # if this row of the table contains a CDS id in our set
                # save the information
                if cds_id in cds_target_set:
                    embl_cds_annotation.append([
                        cds_id, genome_id, cds_to_uniprot[cds_id], start, end
                    ])

        genome_location_table = pd.DataFrame(embl_cds_annotation, columns=[
            "cds", "genome_id", "uniprot_ac", "gene_start", "gene_end"
        ])

        # write the annotation
        return genome_location_table

    def add_full_header(table, alignment):
        new_df = pd.DataFrame()
        for _, row in table.iterrows():
            row_copy = copy(row).to_frame().transpose()

            # for each full_id that corresponds to that uniprot AC
            for full_id in alignment.headers[row_copy.uniprot_ac.values[0]]:
                # create a new row and save that full_id
                #print(full_id)
                row_copy = row_copy.assign(full_id=full_id)
                row_copy = row_copy.assign(short_id=full_id.split()[0])
                new_df = pd.concat([new_df, row_copy])
        return new_df

    def find_possible_partners(gene_location_table_1, gene_location_table_2):

        possible_partners_list = []

        # drop any rows that are missing location information for the CDS
        gene_location_table_1.dropna(axis=0, inplace=True)
        gene_location_table_2.dropna(axis=0, inplace=True)

        # after dropping NAs, coerce remaining locations to integers
        gene_location_table_1[["gene_start", "gene_end"]] = gene_location_table_1[["gene_start", "gene_end"]].astype(
            int)
        gene_location_table_2[["gene_start", "gene_end"]] = gene_location_table_2[["gene_start", "gene_end"]].astype(
            int)

        # drop duplicate rows - speeds up calculation
        gene_location_table_1.drop_duplicates(inplace=True)
        gene_location_table_2.drop_duplicates(inplace=True)

        # group the gene location tables by genome ID
        location_groups_1 = gene_location_table_1.groupby("genome_id")
        location_groups_2 = gene_location_table_2.groupby("genome_id")

        # iterate over EMBL genomes found in the first alignment
        for genome in location_groups_1.groups.keys():

            # if the same genome is found in the second alignment
            if genome in location_groups_2.groups.keys():

                # extract the entries corresponding to the current genome
                gene_location_subset_1 = location_groups_1.get_group(genome)
                gene_location_subset_2 = location_groups_2.get_group(genome)

                # compare all pairs of CDSs
                # that originate from the current genome
                # find the distances between all pairs
                for _, first_cds in gene_location_subset_1.iterrows():
                    for _, second_cds in gene_location_subset_2.iterrows():
                        distance_between_genes = get_distance(
                            (  # location of the first cds
                                first_cds["gene_start"],
                                first_cds["gene_end"]
                            ),
                            (  # location of the second cds
                                second_cds["gene_start"],
                                second_cds["gene_end"]
                            ),
                        )

                        # get the uniprot ids corresponding to the two CDSs
                        uniprot_id_1 = first_cds["short_id"]
                        uniprot_id_2 = second_cds["short_id"]

                        possible_partners_list.append(
                            (uniprot_id_1, uniprot_id_2, distance_between_genes)
                        )

        possible_partners = pd.DataFrame(
            possible_partners_list,
            columns=["uniprot_id_1", "uniprot_id_2", "distance"]
        )

        return possible_partners

    def best_reciprocal_matching(possible_partners):
        """
        Amongst all possible pairings of CDSs
        for two monomer proteins, finds
        those where both sequences are closest
        on the genome to each other
        Parameters
        ----------
        possible_partners : pd.DataFrame
            Columns : uniprot_id_1, uniprot_id_2, distance_between_genes
            Generated by the find_possible_partners function.
            Gives all pairs of uniprot identifiers from alignment 1
            with every uniprot identifier from alignment 2 that is found in
            the same genome, and the number of nucleotides between their
            corresponding coding DNA sequences (CDS)
        Returns
        -------
        pd.DataFrame
            Columns: uniprot_id_1, uniprot_id_2, distance_between_genes
            All pairings of uniprot identifiers that are reciprocally the
            closest to one another in a genome.
        """
        # initialize list to store the matches
        id_pairing_list = []

        # group the table by first and second uniprot identifier
        id_group_1 = possible_partners.groupby("uniprot_id_1")
        id_group_2 = possible_partners.groupby("uniprot_id_2")

        # loop through all sequences in first alignment, and find the closest reciprocal
        # partner for each
        for uniprot_id_1 in id_group_1.groups.keys():

            # get the table that corresponds to the current uniprot id
            id_subset_1 = id_group_1.get_group(uniprot_id_1)

            # what is the closest sequence in second alignment (w.r.t. to genome distance)?
            _index_of_closest = id_subset_1["distance"].idxmin()
            closest_to_uniprot_1 = id_subset_1.loc[_index_of_closest]["uniprot_id_2"]

            # get all of the rows that contain the closest hit in the second alignment
            id_subset_2 = id_group_2.get_group(closest_to_uniprot_1)

            # find the closest sequence in first alignment to the above sequence in second alignment
            _index_of_closest = id_subset_2["distance"].idxmin()
            closest_to_uniprot_2 = id_subset_2.loc[_index_of_closest]["uniprot_id_1"]

            # check if matched sequences are reciprocally the closest on the genome
            if closest_to_uniprot_2 == uniprot_id_1:
                distance = id_subset_1["distance"].min()
                id_pairing_list.append((uniprot_id_1, closest_to_uniprot_1, distance))

        # convert the data to a pandas dataframe
        id_pairing = pd.DataFrame(
            id_pairing_list,
            columns=["uniprot_id_1", "uniprot_id_2", "distance"]
        )

        return id_pairing

    def load_monomer_info(self, alignment):

        cds_ids = self.extract_cds_ids(alignment)

        #print(cds_ids)
        # extract genome location information from ENA

        genome_location_table = self.extract_embl_annotation(cds_ids)

        #print(genome_location_table)

        genome_location_table = Geno_interact.add_full_header(genome_location_table, alignment)

        #print(genome_location_table)
        return genome_location_table

    def get_interactions(self, alignment1, alignment2):

        gene_location_table_1 = self.load_monomer_info(alignment1)
        gene_location_table_2 = self.load_monomer_info(alignment2)

        if gene_location_table_1.empty or gene_location_table_2.empty:
            id_pairing = pd.DataFrame({'id_1': [], 'id_2': []})
        else:
            # find all possible matches
            possible_partners = Geno_interact.find_possible_partners(
                gene_location_table_1, gene_location_table_2
            )

            # find the best reciprocal matches
            id_pairing_unfiltered = Geno_interact.best_reciprocal_matching(possible_partners)

            # filter best reciprocal matches by genome distance threshold
            if params["genome_distance_threshold"]:
                distance_threshold = int(params["genome_distance_threshold"])
                id_pairing = id_pairing_unfiltered[id_pairing_unfiltered.distance.astype('int64') < distance_threshold]
            else:
                id_pairing = id_pairing_unfiltered

            id_pairing.loc[:, "id_1"] = id_pairing.loc[:, "uniprot_id_1"]
            id_pairing.loc[:, "id_2"] = id_pairing.loc[:, "uniprot_id_2"]

        return id_pairing
