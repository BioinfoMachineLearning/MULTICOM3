from collections import OrderedDict, defaultdict
import re
import pandas as pd
from alignment import read_fasta, read_a3m, write_a3m, cal_identities_to_target
from copy import copy

def read_species_annotation_table(annotation):
    """
    Reads in the annotation.csv file and decides which column
    contains the species information - this differs for uniref
    and uniprot alignments. Adds a new column called "species"
    with the correct annotations.
    Note: Uniprot and Uniref databases can have different
    annotations even for the same sequence.
    Parameters
    ----------
    annotation_file : str
        path to annotation file
    Returns
    -------
    pd.DataFrame
        an annotation dataframe with columns id, species, and annotation
    """
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


def extract_header_annotation(a3m_file):
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
    seqs = read_a3m(a3m_file)
    for seq in seqs:
        seq_id = None
        anno = None
        split = seq.split(maxsplit=1)
        if len(split) == 2:
            seq_id, anno = split
        else:
            seq_id = seq
            anno = None
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
    """
    For each species in the alignment, finds the sequence
    from that species that is most similar to the target sequence
    Parameters
    ----------
    similarities : pd.DataFrame
        The contents of identities.csv
    id_to_organism :  pd.DataFrame
        The contents of annotation.csv

    Returns
    -------
    pd.DataFrame
        With columns id, species, identity_to_query.
        Where each row is the sequence in a particular species
        that was the most similar to the target sequence.

    """

    # merge the two data frames
    data = similarities.merge(id_to_organism, on="id")

    # find the most similar in every organism
    ttt = data.sort_values(by="identity_to_query").groupby("species")
    most_similar_in_species = data.sort_values(by="identity_to_query").groupby("species").last()
    most_similar_in_species["species"] = most_similar_in_species.index
    most_similar_in_species = most_similar_in_species.reset_index(drop=True)

    return most_similar_in_species

# def find_paralogs(target_id, id_to_organism, similarities, identity_threshold):
#     """
#     Finds all the sequences in the alignment that originate
#     from the same species as the target_id, if those sequences
#     are below the identity threshold (ie, are diverged from
#     the query)
#     Parameters
#     ----------
#     target_id : str
#         Full identifier of the query sequence
#     similarities : pd.DataFrame
#         The contents of identities.csv
#     id_to_organism :  pd.DataFrame
#         The contents of annotation.csv
#     identity_threshold : float
#         Sequences above this identity to the query are not considered paralogs
#     Returns
#     -------
#     pd.DataFrame
#         with columns id, species, identity_to_query
#         Entries are paralogs found in the same genome as the query id
#     """
#
#     # get all the rows that have an id that contains the
#     # query id. This includes the focus sequence and its hit to
#     # itself in the database.
#     annotation_data = similarities.merge(id_to_organism, on="id")
#     contains_annotation = [target_id in x for x in annotation_data.id]
#     query_hits = annotation_data.loc[contains_annotation , :]
#     # get the species annotation for the query sequence
#     query_species = list(query_hits.species.dropna())
#
#     # get all rows that are from the query species
#     paralogs = annotation_data.query("species == @query_species")
#
#     # confirm that paralogs are below the similarity threshold
#     # ie, are diverged in sequence space from the query
#     paralogs = paralogs.query("identity_to_query < @identity_threshold")
#     return paralogs


# def filter_best_reciprocal(alignment, paralogs, most_similar_in_species, allowed_error=0.02):
#     """
#     Takes in a dictionary of the best hit to each genome
#     Removes sequences that are not the best reciprocal hit to the query sequence
#     Parameters
#     ----------
#     alignment : str
#         Path to sequence alignment file
#     paralogs : pd.DataFrame
#         Rows correspond to paralogs to the query sequence
#         Created by find_paralogs() function
#     most_similar_in_species : pd.DataFrame
#         Contains the id, species name, and percent identity to query
#         for each sequence that was the best hit to the query in its
#         respective species
#     allowed_error : float
#         In order for a sequence to be filtered out of the alignment,
#         it must be more identitical to a paralog sequence than the
#         target sequence by at least this amount
#     Returns
#     -------
#     pd.DataFrame
#         Contains the id, species name, and percenty identity to query
#         for each sequence that was the best reciprocal hit to the query sequence
#     """
#     with open(alignment, "r") as inf:
#         ali = Alignment.from_file(inf)
#
#     # Create an n_paralogs x n_sequences ndarray
#     # where entry i,j is percent identity of paralog i to sequence j
#     # note the identity here will be different than for the unfiltered alignment
#
#     # initialize the matrix
#     identity_mat = np.zeros((len(paralogs), len(ali.ids)), dtype=float)
#
#     for idx, paralog_id in enumerate(paralogs.id):
#         # calculate the % identity of every seq in the alignment to current paralog
#         identities = ali.identities_to(ali[ali.id_to_index[paralog_id]])
#
#         # save the results
#         identity_mat[idx, :] = identities
#
#     indices_to_keep = []
#     # for every sequence in the alignment that is the most similar to the query
#     # in its respective species...
#     for index, row in most_similar_in_species.iterrows():
#
#         # get the index of that sequence in the alignment.
#         alignment_index = ali.id_to_index[row.id]
#
#         # Keep sequences if they are the best reciprocal hit -
#         # i.e., that sequence is not more similar to any paralog
#         # than it is to the query sequence
#         if np.all(identity_mat[:, alignment_index] < row.identity_to_query + allowed_error):
#
#             indices_to_keep.append(index)
#
#     return most_similar_in_species.loc[indices_to_keep, :]

def write_concatenated_alignment(id_pairing, alignment_1, alignment_2, name_1, name_2):
    """
    Concatenate monomer alignments into a complex alignment
    and output to file.

    Parameters
    ----------
    id_pairing : pd.DataFrame
        dataframe with columns id_1 and id_2
        indicating the pairs of sequences to be concatenated
    alignment_1 : str
        Path to alignment file for first monomer alignment
    alignment_2 : str
        Path to alignment file for second monomer alignment
    target_sequence_1 : str
        Target sequence identifier for first monomer alignment
    target_sequence_2 : str
        Target sequence identifier for second monomer alignment
    Returns
    -------
    str
        Header of the concatenated target sequence
    int
        Index of target sequence in the alignment
    Alignment
        the full concatenated alignment
    Alignment
        An alignment of the first monomer sequences with
        only the sequences contained in the concatenated
        alignment
    Alignment
        An alignment of the second monomer sequences with
        only the sequences contained in the concatenated
        aligment
    """

    def _prepare_header(id1, id2):
        # id1_id2
        header_format = "{}_{}"

        concatenated_header = header_format.format(id1, id2)

        return concatenated_header

    sequences_to_write = []  # list of (header,seq1,seq2) tuples

    # load the monomer alignments
    ali_1 = read_a3m(alignment_1, True)
    ali_2 = read_a3m(alignment_2, True)

    # prepare the target sequence
    # target_sequences = (ali_1, ali_2)

    # Target header must end with /1-range for correct focus mode
    # length = len(target_sequences[0]) + len(target_sequences[1])

    target_header = "{}_{}".format(name_1, name_2)

    # store target sequence for writing

    sequences_to_write.append((target_header, ali_1[name_1], ali_2[name_2]))

    # the target sequence is the first in the output file
    target_seq_idx = 0

    # create other headers and sequences
    for id1, id2 in zip(id_pairing.id_1, id_pairing.id_2):

        # prepare the concatenated header
        concatenated_header = _prepare_header(id1, id2)

        # save the information
        sequences_to_write.append(
            (
                concatenated_header,
                ali_1[id1],
                ali_2[id2]
            )
        )

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

    return target_header, target_seq_idx, sequences_full, sequences_monomer_1, sequences_monomer_2


def retrieve_sequence_ids(fileobj, regex=None):
#     """
#     Returns all identifiers in a FASTA alignment;
#     extracts ids based on the given regular expressions
#
#     Note: if multiple regular expressions match the
#     FASTA file header, will return the string extracted
#     by the FIRST match
#
#     Parameters
#     ----------
#     fileobj : file-like object
#         FASTA alignment file
#     regex : list of str, optional (default: None)
#         Regular expression strings to extract sequence ids;
#         if None uses list of default regular expressions
#         to extract Uniprot and UniRef identifiers
#     Returns
#     -------
#     list of str
#         Sequence ids
#     dict
#         Map from sequence id to list of str containing
#         the full sequence headers corresponding to that
#         sequence id
#     """
     if regex is None:
         regex = [
             # example: >UniRef100_H6SNJ6/11-331
             "^Uni\w+\_(\w+).*/",

             # example: >tr|Q1NYN0|Q1NYN0_9FLAO
             "^\w+\|(\w+)\|\w+\/",

             # example: >NQO8_THET8/1-365
             "^(\w+).*/.*$",

             # example: >Q60019|NQO8_THET8/1-365
             "^\w+\|\w+\|(\w+)",
         ]

     sequence_ids = []
     id_to_full_header = defaultdict(list)

     for current_id, _ in read_fasta(fileobj):
         for pattern in regex:
             m = re.match(pattern, current_id)
             # require a non-None match and at least one extracted pattern
             if m and len(m.groups()) > 0:
                 # this extracts the parenthesized match
                 sequence_ids.append(m.group(1))
                 id_to_full_header[m.group(1)].append(current_id)
                 break

     return sequence_ids, id_to_full_header


def get_distance(annotation_1, annotation_2):
    """
    Compute distance between two CDS locations on
    the same genome

    Parameters
    ----------
    annotation_1 : tuple of (int, int)
        genome start and end location for the first CDS
    annotation_2 : tuple of (int, int)
        genome start and end location for the first CDS

    Returns
    -------
    int
        Distance between gene 1 and gene 2 on the
        ENA genome
    """

    # extract the two locations from the annotation
    # sort each so that smaller genome position is first
    location_1 = sorted(annotation_1)
    location_2 = sorted(annotation_2)

    # sort the two annotations so that the one with
    # an earlier start site is first
    x, y = sorted((location_1, location_2))

    # if not overlapping, calculate the distance
    if x[0] <= x[1] < y[0]:
        return y[0] - x[1]

    # if overlapping, return 0
    return 0


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


def find_possible_partners(gene_location_table_1, gene_location_table_2):
    """
    Constructs a dataframe of all possible sequence pairings
    between the two monomer alignments. Best reciprocal matches
    will then be selected using best_reciprocal_matching().
    Parameters
    ----------
    gene_location_table_1: pd.DataFrame
        Dataframe of gene locations for the first
        protein, generated by extract_embl_annotation
    gene_location_table_2: pd.DataFrame
        Dataframe of gene locations for the second
        protein, generated by extract_embl_annotation
    Returns
    -------
    pd.DataFrame
        Columns: uniprot_id_1, uniprot_id_2, distance_between_genes
        Gives all pairs of uniprot identifiers from alignment 1
        with every uniprot identifier from alignment 2 that is found in
        the same genome, and the number of nucleotides between their
        corresponding coding DNA sequences (CDS)
    """

    possible_partners_list = []

    # drop any rows that are missing location information for the CDS
    gene_location_table_1.dropna(axis=0, inplace=True)
    gene_location_table_2.dropna(axis=0, inplace=True)

    # after dropping NAs, coerce remaining locations to integers
    gene_location_table_1[["gene_start", "gene_end"]] = gene_location_table_1[["gene_start", "gene_end"]].astype(int)
    gene_location_table_2[["gene_start", "gene_end"]] = gene_location_table_2[["gene_start", "gene_end"]].astype(int)

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


def plot_distance_distribution(id_pair_to_distance, outfile):
    """
    plots the distribution of genome distances. This is designed
    to run on the output of best_reciprocal_matching().

    Parameters
    ----------
    id_pair_to_distance : pd.DataFrame
        with columns uniprot_id_1, uniprot_id_2, distance_between_genes
    outfile : str
        path to file to save plot
    """

    distances = list(id_pair_to_distance["distance"])
    distances = sorted(distances)

    # Make sure that there is at least one non-nan distance
    # in the data frame
    if len(distances) == 0:
        raise ValueError(
            "No valid distances provided"
        )

    cdf = range(len(distances))

    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.gca()
    ax1.set_xscale("log")
    ax1.set_xlim(xmin=1, xmax=max(distances))
    ax1.set_ylabel("Number of sequences")
    ax1.set_xlabel("Genome distance (bases)")
    ax1.plot(distances, cdf)

    plt.savefig(outfile)


def extract_cds_ids(alignment_file,
                    uniprot_to_embl_table):
    """
    Extracts mapping from set of Uniprot IDs to EMBL
    Coding DNA sequence (CDS) from precomputed ID mapping table.
    Will only include CDSs that can be mapped unambiguously
    to one EMBL genome.

    Parameters
    ----------
    alignment_file : str
        Path to alignment with sequences for which IDs
        should be retrieved
    uniprot_to_embl_table : str
        Path to uniprot to embl mapping database
    Returns
    -------
    list of (str, str)
        A list of Uniprot ac, CDS id pairs
        the CDS id(s) corresponding to each
        Uniprot AC. Uniprot ACs may be repeated
        if there were multiple CDS hits.
    """

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

    # extract identifiers from sequence alignment
    with open(alignment_file) as f:
        sequence_id_list, _ = retrieve_sequence_ids(f)

    # store IDs in set for faster membership checking
    target_ids = set(sequence_id_list)

    # initialize list of list of (uniprot,[(genome,CDS)]) tuples
    # eg [(uniprot1,[(genome1,cds1),(genome2,cds2)])]
    # TODO: I know this is ugly but I need to save the uniprot id
    # for later use, and I need to save the genome so that I can remove
    # cds's that are hit to multiple genomes
    genome_and_cds = []

    # read through the data table line-by-line to improve speed
    with open(uniprot_to_embl_table) as f:
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


def extract_embl_annotation(uniprot_and_cds,
                            ena_genome_location_table):
    """
    Reads genomic location information
    for all input EMBL coding DNA sequences (CDS) corresponding
    to the amino cid sequences in the alignment; creates
    a pd.DataFrame with the following columns:
    cds_id, genome_id, uniprot_ac, gene_start, gene_end
    Each row is a unique CDS. Uniprot ACs may be repeated if one
    Uniprot AC hits multiple CDS.

    Parameters
    ----------
    uniprot_and_cds : list of (str, str)
        A list of uniprot ac, CDS id pairs for which to extract
        genome location
    ena_genome_location_table : str
        Path to ENA genome location database table, which is a
        a tsv file with the following columns:
        cds_id, genome_id, uniprot_ac, genome_start, genome_end
    genome_location_filename : str
        File to write containing CDS location info for
        target sequences
    Returns
    -------
    pd.DataFrame
        Columns: cds, genome_id, uniprot_ac, gene_start, gene_end
        Each row is a unique CDS. Uniprot ACs may be repeated if one
        Uniprot AC hits multiple CDS.
    """

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
    with open(ena_genome_location_table) as inf:
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


def add_full_header(table, alignment_file):
    """
    Add column called full_id
    with the header in the sequence alignment
    that corresponds to the uniprot_AC in the
    genome location table
    Parameters
    ----------
    table : pd.DataFrame
        Columns: cds, genome_id, uniprot_ac, gene_start, gene_end
        Each row is a unique CDS. Uniprot ACs may be repeated if one
        Uniprot AC hits multiple CDS.
    alignment_file : str
        Path to sequence alignment
    Returns
    -------
    pd.DataFrame
        Same as above but with a "full_id"
        column
    """

    with open(alignment_file) as inf:
        _, id_to_header = retrieve_sequence_ids(inf)

    new_df = pd.DataFrame()

    for _, row in table.iterrows():
        row_copy = copy(row).to_frame().transpose()

        # for each full_id that corresponds to that uniprot AC
        for full_id in id_to_header[row_copy.uniprot_ac.values[0]]:
            # create a new row and save that full_id
            row_copy = row_copy.assign(full_id=full_id)
            row_copy = row_copy.assign(short_id=full_id.split()[0])
            new_df = pd.concat([new_df, row_copy])
    return new_df
