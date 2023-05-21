import os, sys
import pandas as pd
import numpy as np

PDB_CHAIN_IDS_UNRELAX = 'BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
PDB_MAX_CHAINS = len(PDB_CHAIN_IDS)  # := 62.


def read_qa_txt_as_df(infile):
    models = []
    scores = []
    for line in open(infile):
        line = line.rstrip('\n')
        contents = line.split()
        if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or \
                contents[0] == "END":
            continue
        model, score = line.split()
        models += [model]
        scores += [float(score)]
    df = pd.DataFrame({'model': models, 'score': scores})
    df = df.sort_values(by=['score'], ascending=False)
    df.reset_index(inplace=True)
    return df


class FastaChain:
    def __init__(self, sequence, description):
        self.sequence = sequence
        self.description = description


def parse_fasta(fasta_string):
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith('>'):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append('')
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line

    return sequences, descriptions


def make_chain_id_map(sequences, descriptions):
    """Makes a mapping from PDB-format chain ID to sequence and description."""
    if len(sequences) != len(descriptions):
        raise ValueError('sequences and descriptions must have equal length. '
                         f'Got {len(sequences)} != {len(descriptions)}.')
    if len(sequences) > PDB_MAX_CHAINS:
        raise ValueError('Cannot process more chains than the PDB format supports. '
                         f'Got {len(sequences)} chains.')
    chain_id_map = {}
    chain_id_seq_map = {}
    for chain_id, sequence, description in zip(
            PDB_CHAIN_IDS, sequences, descriptions):
        chain_id_map[chain_id] = FastaChain(sequence=sequence, description=description)
        chain_id_seq_map[chain_id] = sequence
    return chain_id_map, chain_id_seq_map

def make_chain_id_map_unrelaxed(sequences, descriptions):
    """Makes a mapping from PDB-format chain ID to sequence and description."""
    if len(sequences) != len(descriptions):
        raise ValueError('sequences and descriptions must have equal length. '
                         f'Got {len(sequences)} != {len(descriptions)}.')
    if len(sequences) > PDB_MAX_CHAINS:
        raise ValueError('Cannot process more chains than the PDB format supports. '
                         f'Got {len(sequences)} chains.')
    chain_id_map = {}
    chain_id_seq_map = {}
    for chain_id, sequence, description in zip(
            PDB_CHAIN_IDS_UNRELAX, sequences, descriptions):
        chain_id_map[chain_id] = FastaChain(sequence=sequence, description=description)
        chain_id_seq_map[chain_id] = sequence
    return chain_id_map, chain_id_seq_map


def complete_result(outputdir):
    complete = True
    for i in range(0, 5):
        model = f'{outputdir}/ranked_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def parse_pdb_row(row, param):
    result = ""
    if param == "anum":
        result = row[6:6 + 5]
    elif param == "aname":
        result = row[12:12 + 4]
    elif param == "altloc":
        result = row[16:16 + 1]
    elif param == "rname":
        result = row[17:17 + 3]
    elif param == "rnum":
        result = row[22:22 + 5]
    elif param == "chain":
        result = row[21:21 + 1]
    elif param == "x":
        result = row[30:30 + 8]
    elif param == "y":
        result = row[38:38 + 8]
    elif param == "z":
        result = row[46:46 + 8]
    elif param == "bfactor":
        result = row[60:66]
    else:
        print(f"Invalid row[{row}] or parameter[{param}]")
        return None
    return result.lstrip().rstrip()


def extract_pdb(pdb, newpdb, start, end):
    if start > end:
        die(f"wrong index <start:{start}, end:{end}>\n")

    contents = open(pdb, "r").readlines()
    new_PDBlines = []
    for line in contents:
        if line.find('ATOM') < 0:
            continue
        this_rnum = int(parse_pdb_row(line, "rnum"))
        if this_rnum >= start and this_rnum <= end:
            new_PDBlines += [line]

    # (c) Reindex Chain. Assumptions: non-standard residues removed, alternative locations removed, one model, one chain.
    resCounter = 0
    atomCounter = 0
    prevrNum = "XX"
    with open(newpdb, "w") as f:
        for line in new_PDBlines:
            if line.find('ATOM') < 0:
                continue
            this_rnum = parse_pdb_row(line, "rnum")
            if prevrNum != this_rnum:
                prevrNum = this_rnum
                resCounter = resCounter + 1
            atomCounter = atomCounter + 1
            rnum_string = format(str(resCounter), '>4')
            anum_string = format(str(atomCounter), '>5')
            row = line[0:6] + anum_string + line[11:11 + 5] + " " + \
                  line[17:17 + 3] + " " + " " + rnum_string + " " + line[27:]
            f.write(row)
        f.write("END\n")

def get_stoichiometry_from_fasta(chain_id_map, sequences):
    stoichiometry_counts = []
    processed_sequences = []
    for chain_id in chain_id_map:
        chain_sequence = chain_id_map[chain_id].sequence
        if chain_sequence in processed_sequences:
            continue
        stoichiometry_counts += [len([index for index, sequence in enumerate(sequences) if sequence == chain_sequence])]
        processed_sequences += [chain_sequence]
    
    stoichiometry = ""
    for chain_id, count in zip(PDB_CHAIN_IDS, stoichiometry_counts):
        stoichiometry += f"{chain_id}{count}"
    return stoichiometry

