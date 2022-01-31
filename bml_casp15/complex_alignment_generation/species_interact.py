import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *


class Species_interact:

    def read_species_annotation_table(annotation):

        SPECIES_ANNOTATION_COLUMNS = ["OS", "Tax"]

        data = annotation

        # initialize the column to extract the species information from
        annotation_column = None
        current_num_annotations = 0

        # Determine whether to extract based on the "OS" field
        # or the "Tax" field. Generally, OS is for Uniprot
        for column in SPECIES_ANNOTATION_COLUMNS:
            # if this column contains more non-null values
            if column not in data:
                continue

            num_annotations = sum(data[column].notnull())
            if num_annotations > current_num_annotations:
                # use that column to extract data
                annotation_column = column
                current_num_annotations = num_annotations

        # if we did not find an annotation column, return an error
        if annotation_column is None:
            return None

        # creates a new column called species with the species annotations
        data.loc[:, "species"] = data.loc[:, annotation_column]

        return data[["id", "name", "species"]]

    def extract_header_annotation(alignment):
        columns = [
            ("GN", "gene"),
            ("OS", "organism"),
            ("PE", "existence_evidence"),
            ("SV", "sequence_version"),
            ("n", "num_cluster_members"),
            ("Tax", "taxon"),
            ("RepID", "representative_member")
        ]

        col_to_descr = OrderedDict(columns)
        regex = re.compile("\s({})=".format(
            "|".join(col_to_descr.keys()))
        )

        res = []
        for seq_id in alignment.ids:
            full_header = alignment.headers[seq_id][0]
            anno = None
            if ("GS" in alignment.annotation and
                    seq_id in alignment.annotation["GS"] and
                    "DE" in alignment.annotation["GS"][seq_id]):
                anno = alignment.annotation["GS"][seq_id]["DE"]
            else:
                split = full_header.split(maxsplit=1)
                if len(split) == 2:
                    _, anno = split

            # extract info from line if we got one
            if anno is not None:
                # do split on known field names o keep things
                # simpler than a gigantic full regex to match
                # (some fields are allowed to be missing)
                pairs = re.split(regex, anno)
                pairs = ["id", seq_id, "name"] + pairs

                # create feature-value map
                feat_map = dict(zip(pairs[::2], pairs[1::2]))
                res.append(feat_map)
            else:
                res.append({"id": seq_id})

        df = pd.DataFrame(res)
        return df.reindex(
            ["id", "name"] + list(col_to_descr.keys()),
            axis=1
        )

    def most_similar_by_organism(similarities, id_to_organism):

        # merge the two data frames
        data = similarities.merge(id_to_organism, on="id")

        # find the most similar in every organism
        most_similar_in_species = data.sort_values(by="identity_to_query").groupby("species").last()
        most_similar_in_species["species"] = most_similar_in_species.index
        most_similar_in_species = most_similar_in_species.reset_index(drop=True)

        return most_similar_in_species

    def find_paralogs(target_id, id_to_organism, similarities, identity_threshold):
        """
        Finds all the sequences in the alignment that originate
        from the same species as the target_id, if those sequences
        are below the identity threshold (ie, are diverged from
        the query)
        Parameters
        ----------
        target_id : str
            Full identifier of the query sequence
        similarities : pd.DataFrame
            The contents of identities.csv
        id_to_organism :  pd.DataFrame
            The contents of annotation.csv
        identity_threshold : float
            Sequences above this identity to the query are not considered paralogs
        Returns
        -------
        pd.DataFrame
            with columns id, species, identity_to_query
            Entries are paralogs found in the same genome as the query id
        """

        # get all the rows that have an id that contains the
        # query id. This includes the focus sequence and its hit to
        # itself in the database.
        annotation_data = similarities.merge(id_to_organism, on="id")
        contains_annotation = [target_id in x for x in annotation_data.id]
        query_hits = annotation_data.loc[contains_annotation, :]
        # get the species annotation for the query sequence
        query_species = list(query_hits.species.dropna())

        # get all rows that are from the query species
        paralogs = annotation_data.query("species == @query_species")

        # confirm that paralogs are below the similarity threshold
        # ie, are diverged in sequence space from the query
        paralogs = paralogs.query("identity_to_query < @identity_threshold")
        return paralogs

    def filter_best_reciprocal(alignment, paralogs, most_similar_in_species, allowed_error=0.02):
        """
        Takes in a dictionary of the best hit to each genome
        Removes sequences that are not the best reciprocal hit to the query sequence
        Parameters
        ----------
        alignment : str
            Path to sequence alignment file
        paralogs : pd.DataFrame
            Rows correspond to paralogs to the query sequence
            Created by find_paralogs() function
        most_similar_in_species : pd.DataFrame
            Contains the id, species name, and percent identity to query
            for each sequence that was the best hit to the query in its
            respective species
        allowed_error : float
            In order for a sequence to be filtered out of the alignment,
            it must be more identitical to a paralog sequence than the
            target sequence by at least this amount
        Returns
        -------
        pd.DataFrame
            Contains the id, species name, and percenty identity to query
            for each sequence that was the best reciprocal hit to the query sequence
        """
        # Create an n_paralogs x n_sequences ndarray
        # where entry i,j is percent identity of paralog i to sequence j
        # note the identity here will be different than for the unfiltered alignment

        # initialize the matrix

        identity_mat = np.zeros((len(paralogs), len(alignment.ids)), dtype=float)

        for idx, paralog_id in enumerate(paralogs.id):
            # calculate the % identity of every seq in the alignment to current paralog
            identities = alignment.identities_to(alignment[alignment.id_to_index[paralog_id]])
            # save the results
            identity_mat[idx, :] = identities['identity_to_query']

        indices_to_keep = []
        # for every sequence in the alignment that is the most similar to the query
        # in its respective species...
        for index, row in most_similar_in_species.iterrows():

            # get the index of that sequence in the alignment.
            alignment_index = alignment.id_to_index[row.id]

            # Keep sequences if they are the best reciprocal hit -
            # i.e., that sequence is not more similar to any paralog
            # than it is to the query sequence
            if np.all(identity_mat[:, alignment_index] < row.identity_to_query + allowed_error):
                indices_to_keep.append(index)

        return most_similar_in_species.loc[indices_to_keep, :]

    def load_monomer_info(alignment):

        annotation = Species_interact.extract_header_annotation(alignment)

        annotation_table = Species_interact.read_species_annotation_table(annotation)

        if annotation_table is None:
            return None, None

        similarities = alignment.identities_to(alignment.main_seq)

        most_similar_in_species = Species_interact.most_similar_by_organism(similarities, annotation_table)

        paralogs = Species_interact.find_paralogs(alignment.main_id, annotation_table, similarities, 0.95)

        most_similar_in_species = Species_interact.filter_best_reciprocal(alignment, paralogs, most_similar_in_species)

        return annotation_table, most_similar_in_species

    def get_interactions(alignment1, alignment2):

        annotation_table_1, most_similar_in_species_1 = Species_interact.load_monomer_info(alignment1)
        annotation_table_2, most_similar_in_species_2 = Species_interact.load_monomer_info(alignment2)

        if annotation_table_1 is None or annotation_table_2 is None:
            species_intersection = pd.DataFrame({'id_1': [], 'id_2': []})
        else:
            species_intersection = most_similar_in_species_1.merge(
                most_similar_in_species_2,
                how="inner",  # takes the intersection
                on="species",  # merges on species identifiers
                suffixes=("_1", "_2")
            )

        return species_intersection
