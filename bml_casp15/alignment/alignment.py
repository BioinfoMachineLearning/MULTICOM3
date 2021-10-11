from collections import OrderedDict
import re
import pandas as pd


def read_fasta(fileobj):
    """
    Generator function to read a FASTA-format file
    (includes aligned FASTA, A2M, A3M formats)
    Parameters
    ----------
    fileobj : file-like object
        FASTA alignment file
    Returns
    -------
    generator of (str, str) tuples
        Returns tuples of (sequence ID, sequence)
    """
    current_sequence = ""
    current_id = None

    for line in fileobj:
        # Start reading new entry. If we already have
        # seen an entry before, return it first.
        if line.startswith(">"):
            if current_id is not None:
                yield current_id, current_sequence

            current_id = line.rstrip()[1:]
            current_sequence = ""

        elif not line.startswith(";"):
            current_sequence += line.rstrip()

    # Also do not forget last entry in file
    yield current_id, current_sequence


def read_a3m(a3m_file, only_id=False, max_gap_fraction=0.9, prefilter=None):

    seqs = OrderedDict()
    contents = open(a3m_file, "r").readlines()
    current_id = None
    current_seq = ""

    for line in contents:
        if line.startswith(">"):
            if current_id is not None:
                add = True
                if prefilter is not None and current_id.split('_')[-1].startswith(prefilter):
                    add = False

                if current_seq.count('-') / float(len(current_seq)) > max_gap_fraction:
                    add = False

                if add:
                    seqs[current_id] = current_seq

            current_id = line[1:].rstrip('\n')
            if only_id:
                current_id = current_id.split()[0]
        else:
            current_seq = line.rstrip('\n')

    if current_id is not None:
        seqs[current_id] = current_seq

    return seqs


def write_a3m(sequences, fileobj):
    for seq_id in sequences:
        fileobj.write(">{}\n".format(seq_id))
        fileobj.write(sequences[seq_id] + "\n")


def clean_a3m(infile, outfile):
    contents = open(infile).readlines()
    contents_new = []
    for content in contents:
        if content.startswith('>'):
            contents_new += [content]
        else:
            contents_new += [c for c in content if (c == c.upper() or c == "-" or c == '\n')]

    f = open(outfile, 'w')
    f.writelines(contents_new)
    f.close()


def align_fasta(clustalo_program, fastafile):
    results = os.popen(f"{clustalo_program} -i {fastafile}").read()
    results = results.split("\n")

    alignments = {}
    targetname = ""
    seq = ""
    for i in range(len(results)):
        line = results[i].replace('\n', '').replace('\r', '')
        if line[0:1] == ">":
            if len(targetname) > 0:
                alignments[targetname] = seq
            targetname = line[1:]
            seq = ""
        else:
            seq += line

    if len(targetname) > 0:
        alignments[targetname] = seq

    return alignments


def cal_sequence_identity_from_fasta(clustalo_program, fasta1, fasta2):
    contents1 = open(fasta1, 'r').readlines()
    if len(contents1) < 2:
        die(f"The format in {fasta1} is wrong!\n")

    name1 = contents1[0].replace('\n', '').replace('\r', '').split()[0]
    seq1 = contents1[1].replace('\n', '').replace('\r', '')

    contents2 = open(fasta2, 'r').readlines()
    if len(contents2) < 2:
        die(f"The format in {fasta2} is wrong!\n")

    name2 = contents2[0].replace('\n', '').replace('\r', '').split()[0]
    seq2 = contents2[1].replace('\n', '').replace('\r', '')

    with open(f"{name1[1:]}_{name2[1:]}.fasta", "w") as f:
        f.write(f"{name1}\n{seq1}\n{name2}\n{seq2}")

    alignments = align_fasta(clustalo_program, f"{name1[1:]}_{name2[1:]}.fasta")

    if len(alignments.keys()) < 2:
        die(f"The output for alignments are wrong! {clustalo_program} -i {name1[1:]}_{name2[1:]}.fasta\n")

    aln1 = alignments[name1[1:]]
    aln2 = alignments[name2[1:]]

    seqid = cal_sequence_identity_from_seq(aln1, aln2)

    os.system(f"rm {name1[1:]}_{name2[1:]}.fasta")

    return seqid


def cal_sequence_identity_from_seq(seq1, seq2):
    common_count = 0

    for i in range(len(seq1)):
        char1 = seq1[i:i + 1]
        char2 = seq2[i:i + 1]
        if char1 != '-' and char1 == char2:
            common_count += 1

    seqid = float(common_count) / float(len(seq1))

    return seqid


def cal_identities_to_target(alns, trg):
    seqids = []
    for alnid in alns:
        seqid = cal_sequence_identity_from_seq(alns[alnid], trg)
        seqids += [seqid]
    return pd.DataFrame({"id": [key.split()[0] for key in alns.keys()], "identity_to_query": seqids})
