import re
from collections import namedtuple, OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd

HMMER_PREFIX_WARNING = "# WARNING: seq names have been made unique by adding a prefix of"


class DefaultOrderedDict(OrderedDict):
    """
    Source:
    http://stackoverflow.com/questions/36727877/inheriting-from-defaultddict-and-ordereddict
    Answer by http://stackoverflow.com/users/3555845/daniel
    Maybe this one would be better?
    http://stackoverflow.com/questions/6190331/can-i-do-an-ordered-default-dict-in-python
    """

    def __init__(self, default_factory=None, **kwargs):
        OrderedDict.__init__(self, **kwargs)
        self.default_factory = default_factory

    def __missing__(self, key):
        result = self[key] = self.default_factory()
        return result


def read_fasta(fileobj):
    current_sequence = ""
    current_id = None

    for line in fileobj:
        if line.startswith(">"):
            if current_id is not None:
                yield current_id, current_sequence

            current_id = line.rstrip()[1:]
            current_sequence = ""

        elif not line.startswith(";"):
            current_sequence += line.rstrip()

    yield current_id, current_sequence


def write_fasta(sequences, fileobj):
    for seq_id in sequences:
        fileobj.write(f">{seq_id}\n")
        fileobj.write(sequences[seq_id] + "\n")


def write_aln(sequences, fileobj):
    for seq_id in sequences:
        fileobj.write(sequences[seq_id] + "\n")


ID_EXTRACTION_REGEX = [
    # example: >UniRef100_H6SNJ6/11-331
    "^Uni\w+\_(\w+).*/",

    # example: >UniRef100_H6SNJ6
    "^Uni\w+\_(\w+).*",

    # example: >Q60019|NQO8_THET8/1-365
    "^\w+\|(\w+)\|\w+",

    # example: >tr|Q1NYN0|Q1NYN0_9FLAO
    "^\w+\|(\w+)\|\w+\/",

    # example: MGYP000191463914/8-90
    "^MGYP(\w+).*/",

    #
    # # example: >NQO8_THET8/1-365
    # "^(\w+).*/.*$",
]


def retrieve_sequence_ids(ids, regex=None):
    if regex is None:
        regex = ID_EXTRACTION_REGEX

    sequence_ids = []
    headers = []
    id_to_full_header = defaultdict(list)

    for current_id in ids:
        for pattern in regex:
            m = re.match(pattern, current_id)
            # require a non-None match and at least one extracted pattern
            if m and len(m.groups()) > 0:
                # this extracts the parenthesized match
                sequence_ids.append(m.group(1))
                id_to_full_header[m.group(1)].append(current_id)
                headers += [current_id]
                break

    return sequence_ids, id_to_full_header, headers


# Holds information of a parsed Stockholm alignment file
StockholmAlignment = namedtuple(
    "StockholmAlignment",
    ["seqs", "gf", "gc", "gs", "gr"]
)


def read_stockholm(fileobj, read_annotation=False, raise_hmmer_prefixes=True):
    name_to_sequence = OrderedDict()
    gf = DefaultOrderedDict(list)
    gc = DefaultOrderedDict(str)
    gs = DefaultOrderedDict(lambda: DefaultOrderedDict(list))
    gr = DefaultOrderedDict(lambda: DefaultOrderedDict(str))

    # line counter within current alignment (can be more than one per file)
    i = 0

    # read alignment
    for line in fileobj:
        if i == 0 and not line.startswith("# STOCKHOLM 1.0"):
            raise ValueError(
                "Not a valid Stockholm alignment: "
                "Header missing. {}".format(line.rstrip())
            )

        if raise_hmmer_prefixes and line.startswith(HMMER_PREFIX_WARNING):
            raise ValueError(
                "HMMER added identifier prefixes to alignment because of non-unique "
                "sequence identifiers. Either some sequence identifier is present "
                "twice in the sequence database, or your target sequence identifier is "
                "the same as an identifier in the database. In the first case, please fix "
                "your sequence database. In the second case, please choose a different "
                "sequence identifier for your target sequence that does not overlap with "
                "the sequence database."
            )

        # annotation lines
        if line.startswith("#"):
            if read_annotation:
                if line.startswith("#=GF"):
                    # can have multiple lines for the same feature
                    _, feat, val = line.rstrip().split(maxsplit=2)
                    gf[feat].append(val)
                elif line.startswith("#=GC"):
                    # only one line with the same GC label
                    _, feat, seq = line.rstrip().split(maxsplit=2)
                    gc[feat] += seq
                elif line.startswith("#=GS"):
                    # can have multiple lines for the same feature
                    _, seq_id, feat, val = line.rstrip().split(maxsplit=3)
                    gs[seq_id][feat] = val
                elif line.startswith("#=GR"):
                    # per sequence, only one line with a certain GR feature
                    _, seq_id, feat, seq = line.rstrip().split()
                    gr[seq_id][feat] += seq

            i += 1

        # actual alignment lines
        else:
            splitted = line.rstrip().split(maxsplit=2)
            # there might be empty lines, so check for valid split
            if len(splitted) == 2:
                name, sequence = splitted
                if name not in name_to_sequence:
                    name_to_sequence[name] = ''
                name_to_sequence[name] += sequence
            i += 1

    seqs = OrderedDict()
    query = ''
    keep_columns = []
    for seq_index, seqid in enumerate(name_to_sequence.keys()):
        sequence = name_to_sequence[seqid]
        if seq_index == 0:
            # Gather the columns with gaps from the query
            query = sequence
            keep_columns = [i for i, res in enumerate(query) if res != '-']

        # Remove the columns with gaps in the query from all sequences.
        aligned_sequence = ''.join([sequence[c] for c in keep_columns])

        seqs[seqid] = aligned_sequence

    yield StockholmAlignment(seqs, gf, gc, gs, gr)


def read_a3m(fileobj, inserts="first"):
    seqs = OrderedDict()

    for i, (seq_id, seq) in enumerate(read_fasta(fileobj)):
        # remove any insert gaps that may still be in alignment
        # (just to be sure)
        seq = seq.replace(".", "")

        if inserts == "first":
            # define "spacing" of uppercase columns in
            # final alignment based on target sequence;
            # remaining columns will be filled with insert
            # gaps in the other sequences
            if i == 0:
                uppercase_cols = [
                    j for (j, c) in enumerate(seq)
                    if (c == c.upper() or c == "-")
                ]
                gap_template = np.array(["."] * len(seq))
                filled_seq = seq
            else:
                uppercase_chars = [
                    c for c in seq if c == c.upper() or c == "-"
                ]
                filled = np.copy(gap_template)
                filled[uppercase_cols] = uppercase_chars
                filled_seq = "".join(filled)

        elif inserts == "delete":
            # remove all lowercase letters and insert gaps .;
            # since each sequence must have same number of
            # uppercase letters or match gaps -, this gives
            # the final sequence in alignment
            seq = "".join([c for c in seq if c == c.upper() and c != "."])
        else:
            raise ValueError(
                "Invalid option for inserts: {}".format(inserts)
            )

        seqs[seq_id] = filled_seq

    return seqs


def write_a3m(sequences, fileobj):
    for seq_id in sequences:
        fileobj.write(f">{seq_id}\n")
        fileobj.write(sequences[seq_id] + "\n")


def detect_format(fileobj):
    for i, line in enumerate(fileobj):
        # must be first line of Stockholm file by definition
        if i == 0 and line.startswith("# STOCKHOLM 1.0"):
            return "stockholm"

        # This indicates a FASTA file
        if line.startswith(">"):
            return "fasta"

        # Skip comment lines and empty lines for FASTA detection
        if line.startswith(";") or line.rstrip() == "":
            continue

        # Arriving here means we could not detect format
        return None


def sequences_to_matrix(sequences):
    """
    Transforms a list of sequences into a
    numpy array.
    Parameters
    ----------
    sequences : list-like (str)
        List of strings containing aligned sequences
    Returns
    -------
    numpy.array
        2D array containing sequence alignment
        (first axis: sequences, second axis: columns)
    """
    if len(sequences) == 0:
        raise ValueError("Need at least one sequence")

    N = len(sequences)
    L = len(next(iter(sequences)))
    matrix = np.empty((N, L), dtype=np.str)

    for i, seq in enumerate(sequences):
        if len(seq) != L:
            raise ValueError(
                "Sequences have differing lengths: i={} L_0={} L_i={}".format(
                    i, L, len(seq)
                )
            )

        matrix[i] = np.array(list(seq))

    return matrix


class Alignment:

    def __init__(self, ids, seqs, annotation=None):

        if len(ids) == 0:
            raise ValueError("Alignment file is empty!")

        self.main_id = ids[0]
        self.main_seq = seqs[0]
        self.L = len(seqs[0])

        self.ids, self.headers_dict, self.headers = retrieve_sequence_ids(ids[1:])
        if len(self.ids) == 0:
            self.ids = ids[1:]
            self.headers = defaultdict(list)
            for id in self.ids:
                self.headers_dict[id].append(id)
                self.headers += [id]

        self.seqs = seqs[1:]
        self.N = len(self.seqs)

        if len(self.ids) != self.N:
            print(ids)
            raise ValueError(
                f"Number of sequence IDs and length of alignment do not match: {self.N} and {len(self.ids)}"
            )

        self.id_to_index = {id_: i for i, id_ in enumerate(self.ids)}
        self.matrix = sequences_to_matrix(seqs)

        self.annotation = annotation
        if len(self.annotation) > 1:
            annotation_ids = deepcopy(list(self.annotation["GS"].keys()))
            for annotation_id in annotation_ids:
                for pattern in ID_EXTRACTION_REGEX:
                    m = re.match(pattern, annotation_id)
                    if m and len(m.groups()) > 0:
                        self.annotation["GS"][m.group(1)] = self.annotation["GS"].pop(annotation_id)
                        break

    @classmethod
    def from_dict(cls, seq_dict, annotation=None):
        return cls(list(seq_dict.keys()), list(seq_dict.values()), annotation)

    @classmethod
    def from_file(cls, fileobj, format="fasta",
                  a3m_inserts="first", raise_hmmer_prefixes=True):

        annotation = {}
        # read in sequence alignment from file

        if format == "fasta":
            seqs = OrderedDict()
            for seq_id, seq in read_fasta(fileobj):
                seqs[seq_id] = seq
        elif format == "stockholm":
            # only reads first Stockholm alignment contained in file
            ali = next(
                read_stockholm(
                    fileobj, read_annotation=True,
                    raise_hmmer_prefixes=raise_hmmer_prefixes
                )
            )
            seqs = ali.seqs
            annotation["GF"] = ali.gf
            annotation["GC"] = ali.gc
            annotation["GS"] = ali.gs
            annotation["GR"] = ali.gr
        elif format == "a3m":
            seqs = read_a3m(fileobj, inserts=a3m_inserts)
        else:
            raise ValueError("Invalid alignment format: {}".format(format))

        return cls.from_dict(seqs, annotation)

    def __getitem__(self, index):
        """
        .. todo::
            eventually this should allow fancy indexing and offer the functionality of select()
        """
        if index in self.id_to_index:
            return self.seqs[self.id_to_index[index]]
        elif index in range(self.N):
            return self.seqs[index]
        else:
            raise KeyError(
                "Not a valid index for sequence alignment: {}".format(index)
            )

    def __len__(self):
        return self.N

    def count(self, char, axis="pos", normalize=True):
        """
        Count occurrences of a character in the sequence
        alignment.
        .. note::
            The counts are raw counts not adjusted for
            sequence redundancy.
        Parameters
        ----------
        char : str
            Character which is counted
        axis : {"pos", "seq"}, optional (default="pos")
            Count along positions or sequences
        normalize : bool, optional (default=True)
            Normalize count for length of axis (i.e. relative count)
        Returns
        -------
        np.array
            Vector containing counts of char along the axis
        Raises
        ------
        ValueError
            Upon invalid axis specification
        """
        if axis == "pos":
            naxis = 0
        elif axis == "seq":
            naxis = 1
        else:
            raise ValueError("Invalid axis: {}".format(axis))

        c = np.sum(self.matrix == char, axis=naxis)
        if normalize:
            c = c / self.matrix.shape[naxis]

        return c

    def cal_sequence_identity_from_seq(seq1, seq2):
        common_count = 0
        for i in range(len(seq1)):
            char1 = seq1[i:i + 1]
            char2 = seq2[i:i + 1]
            if char1 != '-' and char1 == char2:
                common_count += 1
        seqid = float(common_count) / float(len(seq1))
        return seqid

    def identities_to(self, trg):
        seqids = []
        for seq in self.seqs:
            seqids += [Alignment.cal_sequence_identity_from_seq(seq, trg)]
        return pd.DataFrame({"id": self.ids, "identity_to_query": seqids})

    def write(self, fileobj, format="fasta"):
        """
        Write an alignment to a file.
        Parameters
        ----------
        fileobj : file-like object
            File to which alignment is saved
        format : {"fasta", "aln", "a3m"}
            Output format for alignment
        width : int
            Column width for fasta alignment
        Raises
        ------
        ValueError
            Upon invalid file format specification
        """
        if format == "fasta":
            write_fasta(self.seqs, fileobj)
        elif format == "a3m":
            write_a3m(self.seqs, fileobj)
        elif format == "aln":
            write_aln(self.seqs, fileobj)
        else:
            raise ValueError(
                "Invalid alignment format: {}".format(format)
            )
