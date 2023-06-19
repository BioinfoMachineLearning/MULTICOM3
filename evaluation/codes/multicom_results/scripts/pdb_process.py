import os
import sys
import argparse
from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser
import warnings
import numpy as np

warnings.filterwarnings('ignore')
#'TER', 'END'
def delete_elements(fhandle, element='TER'):
    for line in fhandle:
        if line.startswith(element):
                continue
        yield line
    
def pad_line(line):
    """Helper function to pad line to 80 characters in case it is shorter"""
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    if '\n' not in line[:81]:
        return line
    return line[:81]  # 80 + newline character

#starting_resid = -1
def renumber_residues(fhandle, starting_resid):
    """Resets the residue number column to start from a specific number.
    """
    _pad_line = pad_line
    prev_resid = None  # tracks chain and resid
    resid = starting_resid - 1  # account for first residue
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        line = _pad_line(line)
        if line.startswith(records):
            line_resuid = line[17:27]
            if line_resuid != prev_resid:
                prev_resid = line_resuid
                resid += 1
                if resid > 9999:
                    emsg = 'Cannot set residue number above 9999.\n'
                    sys.stderr.write(emsg)
                    sys.exit(1)

            yield line[:22] + str(resid).rjust(4) + line[26:]

        else:
            yield line

#starting_value = -1
def renumber_atom_serials(fhandle, starting_value):
    """Resets the atom serial number column to start from a specific number.
    """

    # CONECT 1179  746 1184 1195 1203
    fmt_CONECT = "CONECT{:>5s}{:>5s}{:>5s}{:>5s}{:>5s}" + " " * 49 + "\n"
    char_ranges = (slice(6, 11), slice(11, 16),
                   slice(16, 21), slice(21, 26), slice(26, 31))

    serial_equiv = {'': ''}  # store for conect statements

    serial = starting_value
    records = ('ATOM', 'HETATM')
    for line in fhandle:
        if line.startswith(records):
            serial_equiv[line[6:11].strip()] = serial
            yield line[:6] + str(serial).rjust(5) + line[11:]
            serial += 1
            if serial > 99999:
                emsg = 'Cannot set atom serial number above 99999.\n'
                sys.stderr.write(emsg)
                sys.exit(1)

        elif line.startswith('ANISOU'):
            # Keep atom id as previous atom
            yield line[:6] + str(serial - 1).rjust(5) + line[11:]

        elif line.startswith('CONECT'):
            # 6:11, 11:16, 16:21, 21:26, 26:31
            serials = [line[cr].strip() for cr in char_ranges]

            # If not found, return default
            new_serials = [str(serial_equiv.get(s, s)) for s in serials]
            conect_line = fmt_CONECT.format(*new_serials)

            yield conect_line
            continue

        elif line.startswith('MODEL'):
            serial = starting_value
            yield line

        elif line.startswith('TER'):
            yield line[:6] + str(serial).rjust(5) + line[11:]
            serial += 1

        else:
            yield line

def alter_chain(fhandle, chain_id):
    """Sets the chain identifier column in all ATOM/HETATM records to a value.
    """
    _pad_line = pad_line
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        if line.startswith(records):
            line = _pad_line(line)
            yield line[:21] + chain_id + line[22:]
            # print(line[:21] + chain_id + line[22:])
        else:
            yield line
            # print(line)

def write_pdb_file(new_pdb, pdb_file):
    if os.path.exists(pdb_file):
        os.remove(pdb_file)
    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                open(pdb_file, 'a').write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        open(pdb_file, 'a').write(''.join(_buffer))
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

def get_sequence_from_pdb(pdb_file):
    fh = open(pdb_file, 'r')
    sequence_list = []
    for record in SeqIO.parse(fh, 'pdb-atom'):
        sequence_list.append(str(record.seq))
    return sequence_list

def remove_head(fhandle):
    """Remove head information.
    """
    records = ('ATOM', 'TER', 'END')
    for line in fhandle:
        if not line.startswith(records):
            continue
        yield line

def delete_residues(fhandle, residue_range, step=1):
    """Deletes residues within a certain numbering range.
    """
    prev_res = None
    res_counter = -1
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        if line.startswith(records):

            res_id = line[21:26]  # include chain ID
            if res_id != prev_res:
                prev_res = res_id
                res_counter += 1

            if int(line[22:26]) in residue_range and res_counter % step == 0:
                continue

        yield line

def keep_residues(fhandle, residue_range):
    """Deletes residues within a certain numbering range.
    """
    prev_res = None
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        if line.startswith(records):

            res_id = line[21:26]  # include chain ID
            if res_id != prev_res:
                prev_res = res_id

            if int(line[22:26]) not in residue_range:
                continue

        yield line

def splitComplex2chains(src_pdb_file, dst_pdb_path = None):
    line_list = open(src_pdb_file, 'r').readlines()
    dir_name = os.path.dirname(src_pdb_file)
    name = os.path.basename(src_pdb_file).split('.')[0]
    dst_pdb_files = []
    monomer_list=[]
    line_last = ''
    count = 0
    for line in line_list:
        monomer_list.append(line)
        if 'TER' in line:
            monomer_list.append('END\n')
            if dst_pdb_path is None and os.path.exists(dst_pdb_path):
                monomer_pdb_file = f'{dir_name}/{name}_{count}.pdb'
            else:
                monomer_pdb_file = f'{dst_pdb_path}/{name}_{count}.pdb'

            open(monomer_pdb_file, 'w').write(''.join(monomer_list))
            dst_pdb_files.append(monomer_pdb_file)
            monomer_list = []
            count += 1
        elif line == line_list[-1] and 'TER' not in line_list[-2]:
            monomer_list.append('END\n')
            if dst_pdb_path is None and os.path.exists(dst_pdb_path):
                monomer_pdb_file = f'{dir_name}/{name}_{count}.pdb'
            else:
                monomer_pdb_file = f'{dst_pdb_path}/{name}_{count}.pdb'

            open(monomer_pdb_file, 'w').write(''.join(monomer_list))
            dst_pdb_files.append(monomer_pdb_file)
            monomer_list = []
        else:
            line_last = line
    return dst_pdb_files


def get_heavyatom_dist_from_pdbfile(pdb_file, length):
    seq_name = os.path.basename(pdb_file).split('.')[0]
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(seq_name, pdb_file)

    model = structure[0]

    chain_id = list(model.child_dict.keys())
    chain = model[chain_id[0]]

    h_dist_map = np.zeros((length, length))

    L = len(chain)
    print(L, length)
    if L != length:
        print('got length error of %s, should be %s, got %s'%(pdb_file, length, L))
        return h_dist_map

    chain_dict = list(chain.child_dict.keys())
    for i in range(0, length):
        for j in range(i, length):
            if i == j:
                continue
            res_i = chain[chain_dict[i]]
            res_j = chain[chain_dict[j]]
            dist_list =[]
            for atom_i in res_i:
                for atom_j in res_j:
                    if ('C' in atom_i.name or 'N' in atom_i.name or 'O' in atom_i.name or 'S' in atom_i.name) and \
                        ('C' in atom_j.name or 'N' in atom_j.name or 'O' in atom_j.name or 'S' in atom_j.name):
                        dist_list.append(atom_i-atom_j)
                    else:
                        continue
            min_dist = np.min(dist_list) 
            h_dist_map[i, j] = min_dist
    h_dist_map += h_dist_map.T
    return h_dist_map

def get_realdist_from_pdbfile(pdb_file, length):
    seq_name = os.path.basename(pdb_file).split('.')[0]
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(seq_name, pdb_file)
    model = structure[0]
    chain_id = list(model.child_dict.keys())
    chain = model[chain_id[0]]

    cbcb_dist_mat = np.zeros((length, length))
    L = len(chain)
    print(L, length)
    if L != length:
        print('got length error of %s, should be %s, got %s'%(pdb_file, length, L))
        return cbcb_dist_mat
    chain_dict = list(chain.child_dict.keys())
    for i in range(0, length):
        for j in range(i, length):
            if i == j:
                continue
            residue_i = chain[chain_dict[i]]
            residue_j = chain[chain_dict[j]]
            try:
                n_i  = residue_i["N"]
                ca_i = residue_i["CA"]
                c_i  = residue_i["C"]
                n_j  = residue_j["N"]
                ca_j = residue_j["CA"]
                c_j  = residue_i["C"]
            except KeyError as e:
                continue
            try:
                cb_i = residue_i["CB"]
            except KeyError as e:
                # recreate Cb given N,Ca,C
                cb_i = ca_i.copy()
                b = ca_i.get_vector() - n_i.get_vector()
                c = c_i.get_vector() - ca_i.get_vector()
                a = np.cross(b._ar, c._ar)
                # the calculation is quote from trRosetta
                cb_i.set_coord(-0.58273431 * a + 0.56802827 * b._ar - 0.54067466 * c._ar + ca_i.get_vector()._ar)
            try:
                cb_j = residue_j["CB"]
            except KeyError as e:
                cb_j = ca_j.copy()
                b = ca_j.get_vector() - n_j.get_vector()
                c = c_j.get_vector() - ca_j.get_vector()
                a = np.cross(b._ar, c._ar)
                cb_j.set_coord(-0.58273431 * a + 0.56802827 * b._ar - 0.54067466 * c._ar + ca_j.get_vector()._ar)

                cbcb_dist = cb_i - cb_j
                cbcb_dist_mat[i, j] = cbcb_dist
                cbcb_dist_mat[j, i] = cbcb_dist
    cbcb_dist_mat += cbcb_dist_mat.T
    return cbcb_dist_mat

#process complex_pdb_file 
def process_pdbfile(in_file, out_file):
    fhandle = open(in_file, 'r')
    fhandle = delete_elements(fhandle, 'TER')
    fhandle = alter_chain(fhandle, 'A')
    fhandle = renumber_residues(fhandle, 1)
    fhandle = renumber_atom_serials(fhandle, 1)
    write_pdb_file(fhandle, out_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='pdb process')
    parser.add_argument('-p', '--pdb_file', help="pdb file", type=str, required=True)
    parser.add_argument('-o', '--out_file', help="output file", type=str, required=True)

    args = parser.parse_args()  
    pdb_file = args.pdb_file
    out_file = args.out_file
    fhandle = open(pdb_file, 'r')

    fhandle = delete_elements(fhandle, 'TER')
    fhandle = alter_chain(fhandle, 'A')
    fhandle = renumber_residues(fhandle, 1)
    fhandle = renumber_atom_serials(fhandle, 1)
    write_pdb_file(fhandle, out_file)
