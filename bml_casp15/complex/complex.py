from collections import OrderedDict, defaultdict
import re
import pandas as pd
from bml_casp15.alignment.alignment import *
from copy import copy


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
    for id1, id2 in zip(id_pairing.id_1, id_pairing.id_2):

        # prepare the concatenated header
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

    return target_header, sequences_full, sequences_monomer_1, sequences_monomer_2