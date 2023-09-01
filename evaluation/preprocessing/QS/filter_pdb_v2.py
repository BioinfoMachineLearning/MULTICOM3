from itertools import chain
import os, sys
import numpy as np
import argparse
from pdb_process import *
from bml_casp15.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs

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

#starting_resid = -1
def renumber_residues_with_order(fhandle, reorder_indices):
    """Resets the residue number column to start from a specific number.
    """
    _pad_line = pad_line
    prev_resid = None  # tracks chain and resid
    counter = -1
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        line = _pad_line(line)
        if line.startswith(records):
            line_resuid = line[17:27]
            if line_resuid != prev_resid:
                prev_resid = line_resuid
                counter += 1
                if reorder_indices[counter] > 9999:
                    emsg = 'Cannot set residue number above 9999.\n'
                    sys.stderr.write(emsg)
                    sys.exit(1)

            yield line[:22] + str(reorder_indices[counter]).rjust(4) + line[26:]

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


def reindex_pdb_file(in_file, out_file, keep_indices, reorder_indices = []):
    # print(keep_indices)
    fhandle = open(in_file, 'r')
    # fhandle = renumber_residues(fhandle, 1)
    # fhandle = renumber_atom_serials(fhandle, 1)
    fhandle = keep_residues(fhandle, keep_indices)
    if len(reorder_indices) > 0:
        fhandle = renumber_residues_with_order(fhandle, reorder_indices)
    else:
        fhandle = renumber_residues(fhandle, 1)
    fhandle = renumber_atom_serials(fhandle, 1)
    write_pdb_file(fhandle, out_file)


# s input str
# c search char
def find_all(s, c):
    idx = s.find(c)
    while idx != -1:
        yield idx
        idx = s.find(c, idx + 1)


def align_sequences(clustalw_program, casp_seq, native_seq, pdb_seq, outfile, workdir):

    with open(workdir + '/' + outfile, 'w') as fw:
        if len(pdb_seq) == 0:
            fw.write(f"%CASP\n{casp_seq}\n%NATIVE\n{native_seq}")
        else:    
            fw.write(f"%CASP\n{casp_seq}\n%NATIVE\n{native_seq}\n%PDB\n{pdb_seq}")

    cmd = f"{clustalw_program}  -MATRIX=BLOSUM -TYPE=PROTEIN -INFILE={workdir}/{outfile} -OUTFILE={workdir}/{outfile}.out >/dev/null 2>&1"
    os.system(cmd)

    casp_align_seq = ""
    native_align_seq = ""
    pdb_align_seq = ""
    for line in open(workdir + '/' + outfile + '.out'):
        if line[0:4] == "CASP":
            casp_align_seq += line.split()[1].rstrip('\n')
        if line[0:6] == "NATIVE":
            native_align_seq += line.split()[1].rstrip('\n')
        if line[0:3] == "PDB":
            pdb_align_seq += line.split()[1].rstrip('\n')

    return casp_align_seq, native_align_seq, pdb_align_seq


def cal_sequence_identity_from_seq(seq1, seq2):
    common_count = 0
    for i in range(len(seq1)):
        char1 = seq1[i:i + 1]
        char2 = seq2[i:i + 1]
        if char1 != '-' and char2 != 'X' and char1 == char2:
            common_count += 1
    seqid = float(common_count) / float(len([1 for i in range(len(seq2)) if seq2[i] != '-' and seq2[i] != 'X']))
    return seqid


def cal_sequence_identity(clustalw_program, seq1, seq2, workdir):
    with open(workdir + '/tmp.fasta', 'w') as fw:
        fw.write(f"%SEQ1\n{seq1}\n%SEQ2\n{seq2}")

    cmd = f"{clustalw_program}  -MATRIX=BLOSUM -TYPE=PROTEIN -INFILE={workdir}/tmp.fasta -OUTFILE={workdir}/tmp.fasta.out >/dev/null 2>&1"
    os.system(cmd)

    seq1 = ""
    seq2 = ""
    for line in open(workdir + '/tmp.fasta.out'):
        if line[0:4] == "SEQ1":
            seq1 += line.split()[1].rstrip('\n')
        if line[0:4] == "SEQ2":
            seq2 += line.split()[1].rstrip('\n')

    return cal_sequence_identity_from_seq(seq1, seq2)


def get_sequence(inpdb):
    """Enclosing logic in a function to simplify code"""

    seq_to_res_mapping = []
    res_codes = [
        # 20 canonical amino acids
        ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
        ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
        ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
        ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
        ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
        # Non-canonical amino acids
        # ('MSE', 'M'), ('SOC', 'C'),
        # Canonical xNA
        ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
        ('  T', 'T'),
    ]

    three_to_one = dict(res_codes)
    # _records = set(['ATOM  ', 'HETATM'])
    _records = set(['ATOM  '])

    sequence = []
    read = set()
    for line in open(inpdb):
        line = line.strip()
        if line[0:6] in _records:
            resn = line[17:20]
            resi = line[22:26]
            icode = line[26]
            r_uid = (resn, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
            else:
                continue
            aa_resn = three_to_one.get(resn, 'X')
            sequence.append(aa_resn)
            seq_to_res_mapping += [int(resi)]

    return {'sequence': ''.join(sequence), 'mapping': seq_to_res_mapping}


def split_pdb(complex_pdb, outdir):
    pdbs = {}
    pre_chain = None
    i = 0
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if pre_chain is None:
            pre_chain = chain_name
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            pdbs[chain_name] = outdir + '/' + chain_name + '.pdb'
            # fw.write(line[:21] + ' ' + line[22:])
            fw.write(line)
        elif chain_name == pre_chain:
            # fw.write(line[:21] + ' ' + line[22:])
            fw.write(line)
        else:
            fw.close()
            i = i + 1
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            pdbs[chain_name] = outdir + '/' + chain_name + '.pdb'
            # fw.write(line[:21] + ' ' + line[22:])
            fw.write(line)
            pre_chain = chain_name
    fw.close()
    return pdbs


def remove_head(fhandle):
    """Remove head information.
    """
    records = ('ATOM', 'TER', 'END')
    for line in fhandle:
        if not line.startswith(records):
            continue
        yield line

def remove_pdb_file_head(in_file, out_file):
    fhandle = open(in_file, 'r')
    fhandle = remove_head(fhandle)
    write_pdb_file(fhandle, out_file)


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


def combine_atom_file(atom_file_list, complex_file):
    atom = []
    for file in atom_file_list:
        atom += open(file, 'r').readlines()[:-1]
    with open(complex_file, 'w') as myfile:
        myfile.write(''.join(atom))
        myfile.write('END\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-f", "--fastadir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-t", "--truepdbdir", help="output folder for the pdb files", type=str, required=True)
    parser.add_argument("-tb", "--tarballdir", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-tmp", "--tmpdir", type=str, required=True)
    parser.add_argument("-c", "--clustalw", type=str, required=True)

    args = parser.parse_args()

    makedir_if_not_exists(args.tmpdir)

    makedir_if_not_exists(args.outdir)

    pdboutdir = args.outdir + '/pdb_filtered'
    makedir_if_not_exists(pdboutdir)

    predoutdir = args.outdir + '/pred_filtered'
    makedir_if_not_exists(predoutdir)

    for true_pdb in os.listdir(args.truepdbdir):
        if true_pdb.find('.pdb') < 0:
            continue

        print(f"processing {true_pdb}")
        targetname = true_pdb[0:true_pdb.find('.pdb')]
        tmpdir = args.tmpdir + '/' + targetname
        makedir_if_not_exists(tmpdir)

        makedir_if_not_exists(tmpdir + '/native')
        true_pdb_chain_files = split_pdb(args.truepdbdir + '/' + true_pdb, tmpdir + '/native')
        
        fasta_file = args.fastadir + '/' + targetname + '.fasta'
        with open(fasta_file) as f:
            input_fasta_str = f.read()
        input_seqs, input_descs = parse_fasta(input_fasta_str)
        chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                           descriptions=input_descs)

        print(chain_id_map)
        native_filtered_pdbs = []
        predoutdir_target = predoutdir + '/' + targetname
        makedir_if_not_exists(predoutdir_target)

        for ori_pred_pdb_file in os.listdir(args.tarballdir + '/' + targetname):
            ori_pred_pdb_file = f'{args.tarballdir}/{targetname}/{ori_pred_pdb_file}'
            if os.path.isdir(ori_pred_pdb_file):
                continue
            print(f"processing {ori_pred_pdb_file}")
            pred_name = os.path.basename(ori_pred_pdb_file).split('.')[0]
            pred_pdb_file = f'{tmpdir}/{pred_name}.pdb'
            remove_pdb_file_head(ori_pred_pdb_file, pred_pdb_file)

            model_workdir = tmpdir + '/' + pred_name
            makedir_if_not_exists(model_workdir)
            pred_pdb_chain_files = split_pdb(pred_pdb_file, model_workdir)

            # check chain ids
            pass_check = True
            if len(true_pdb_chain_files) < len(pred_pdb_chain_files):
                print(f"The chains of {true_pdb}:{len(true_pdb_chain_files)} doesn't match with chains in the {ori_pred_pdb_file}:{len(pred_pdb_chain_files)}!")
                continue

            # for chain_id in true_pdb_chain_files:
            #     if chain_id not in pred_pdb_chain_files:
            #         print(f"chain {chain_id} in the native structure cannot find in the {ori_pred_pdb_file}")
            #         pass_check = False
            #         break

            for chain_id in pred_pdb_chain_files:
                if chain_id not in true_pdb_chain_files:
                    print(f"chain {chain_id} in the native structure cannot find in the {ori_pred_pdb_file}")
                    pass_check = False
                    break

            if not pass_check:
                continue

            # check chain order
            pass_order = True
            for chain_id in pred_pdb_chain_files:
                if chain_id not in chain_id_map:
                    pass_order = False
                    break

                seqid = cal_sequence_identity(args.clustalw, 
                                        get_sequence(pred_pdb_chain_files[chain_id])['sequence'], 
                                        chain_id_map[chain_id].sequence, 
                                        model_workdir)
                if seqid < 0.2:
                    print(chain_id)
                    print(get_sequence(pred_pdb_chain_files[chain_id])['sequence'])
                    print(chain_id_map[chain_id].sequence)
                    print(seqid)
                    pass_order = False
                    break
            

            if not pass_order:
                print(f"The chain order of {ori_pred_pdb_file} doesn't match with the casp fasta file!")
                continue         

            true_pdb_file_reindex_list = []
            pred_pdb_file_reindex_list = []
            for chain_id in true_pdb_chain_files:
                true_pdb_file = true_pdb_chain_files[chain_id]

                if chain_id not in pred_pdb_chain_files:
                    
                    casp_seq = chain_id_map[chain_id].sequence

                    native_seq_res = get_sequence(true_pdb_chain_files[chain_id])
                   
                    casp_align_seq, native_align_seq, _ = align_sequences(args.clustalw, 
                                                                        casp_seq, native_seq_res['sequence'], "", 
                                                                        f"{pred_name}_{chain_id}.fasta", model_workdir)

                    index_casp_temp = list(find_all(casp_align_seq, '-'))  # default the 0 is the best paired
                    index_native_temp = list(find_all(native_align_seq, '-'))

                    indexDel = index_casp_temp + index_native_temp
                    indices_keep = [i for i in range(len(casp_align_seq)) if i not in indexDel]

                    valid_char_counter = -1
                    native_indices_keep = []
                    for i in range(len(casp_align_seq)):
                        if native_align_seq[i] != '-':
                            valid_char_counter += 1
                        if i in indices_keep:
                            native_indices_keep += [native_seq_res['mapping'][valid_char_counter]]

                    true_pdb_file_reindex_temp = f'{model_workdir}/{targetname}_{chain_id}.pdb'
                    reindex_pdb_file(true_pdb_file, true_pdb_file_reindex_temp, native_indices_keep)
                    true_pdb_file_reindex_list.append(true_pdb_file_reindex_temp)

                else:
                    pred_pdb_file = pred_pdb_chain_files[chain_id]
                    
                    casp_seq = chain_id_map[chain_id].sequence

                    native_seq_res = get_sequence(true_pdb_chain_files[chain_id])
                    pdb_seq_res = get_sequence(pred_pdb_chain_files[chain_id])

                    casp_align_seq, native_align_seq, pdb_align_seq = align_sequences(args.clustalw, 
                                                                        casp_seq, native_seq_res['sequence'], pdb_seq_res['sequence'], 
                                                                        f"{pred_name}_{chain_id}.fasta", model_workdir)

                    index_casp_temp = list(find_all(casp_align_seq, '-'))  # default the 0 is the best paired
                    index_native_temp = list(find_all(native_align_seq, '-'))

                    indexDel = index_casp_temp + index_native_temp
                    indices_keep = [i for i in range(len(casp_align_seq)) if i not in indexDel]

                    valid_char_counter = -1
                    native_indices_keep = []
                    for i in range(len(casp_align_seq)):
                        if native_align_seq[i] != '-':
                            valid_char_counter += 1
                        if i in indices_keep:
                            native_indices_keep += [native_seq_res['mapping'][valid_char_counter]]

                    true_pdb_file_reindex_temp = f'{model_workdir}/{targetname}_{chain_id}.pdb'
                    reindex_pdb_file(true_pdb_file, true_pdb_file_reindex_temp, native_indices_keep)
                    true_pdb_file_reindex_list.append(true_pdb_file_reindex_temp)

                    pdb_indices_keep = []
                    pdb_indices_order = []
                    pdb_indices_counter = 0
                    valid_char_counter = -1
                    for i in range(len(casp_align_seq)):
                        if native_align_seq[i] != '-':
                            pdb_indices_counter += 1
                        if pdb_align_seq[i] != '-':
                            valid_char_counter += 1
                        if i in indices_keep and pdb_align_seq[i] != '-':
                            pdb_indices_keep += [pdb_seq_res['mapping'][valid_char_counter]]
                            pdb_indices_order += [pdb_indices_counter]

                    pred_pdb_file_reindex_temp = f'{model_workdir}/{pred_name}_{chain_id}.pdb'
                    reindex_pdb_file(pred_pdb_file, pred_pdb_file_reindex_temp, pdb_indices_keep, pdb_indices_order)
                    pred_pdb_file_reindex_list.append(pred_pdb_file_reindex_temp) 

            true_pdb_file_reindex = f'{model_workdir}/{targetname}_filtered.pdb'
            pred_pdb_file_reindex = f'{model_workdir}/{pred_name}_filtered.pdb'

            combine_atom_file(true_pdb_file_reindex_list, true_pdb_file_reindex)
            combine_atom_file(pred_pdb_file_reindex_list, pred_pdb_file_reindex)
            print(f"### save final filter pdb in {true_pdb_file_reindex}")
            print(f"### save final filter pdb in {pred_pdb_file_reindex}")

            native_filtered_pdbs += [true_pdb_file_reindex]

            os.system(f"cp {pred_pdb_file_reindex} {predoutdir_target}/{pred_name}")

        if len(native_filtered_pdbs) > 0:
            pdb_filtered = native_filtered_pdbs[0]

            for i in range(1, len(native_filtered_pdbs)):
                cmd = f"diff {native_filtered_pdbs[0]} {native_filtered_pdbs[i]}"
                result = os.popen(cmd).read().strip()
                if len(result) > 0:
                    raise Exception(f"The contents between {native_filtered_pdbs[0]} and {native_filtered_pdbs[i]} are different!")
            
            os.system(f"cp {native_filtered_pdbs[0]} {pdboutdir}")
        else:
            print("Failed to filter the native structure")

        
