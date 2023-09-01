import argparse, os
import sys
import logging
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO

def int_to_chain(i, base=62):
    """
    int_to_chain(int,int) -> str

    Converts a positive integer to a chain ID. Chain IDs include uppercase
    characters, numbers, and optionally lowercase letters.

    i = a positive integer to convert
    base = the alphabet size to include. Typically 36 or 62.
    """
    if i < 0:
        raise ValueError("positive integers only")
    if base < 0 or 62 < base:
        raise ValueError("Invalid base")

    quot = int(i) // base
    rem = i % base
    if rem < 26:
        letter = chr(ord("A") + rem)
    elif rem < 36:
        letter = str(rem - 26)
    else:
        letter = chr(ord("a") + rem - 36)
    if quot == 0:
        return letter
    else:
        return int_to_chain(quot - 1, base) + letter


class OutOfChainsError(Exception): pass


def rename_chains(structure):
    """Renames chains to be one-letter chains

    Existing one-letter chains will be kept. Multi-letter chains will be truncated
    or renamed to the next available letter of the alphabet.

    If more than 62 chains are present in the structure, raises an OutOfChainsError

    Returns a map between new and old chain IDs, as well as modifying the input structure
    """
    next_chain = 0  #
    # single-letters stay the same
    chainmap = {c.id: c.id for c in structure.get_chains() if len(c.id) == 1}
    for o in structure.get_chains():
        if len(o.id) != 1:
            if o.id[0] not in chainmap:
                chainmap[o.id[0]] = o.id
                o.id = o.id[0]
            else:
                c = int_to_chain(next_chain)
                while c in chainmap:
                    next_chain += 1
                    c = int_to_chain(next_chain)
                    if next_chain >= 62:
                        raise OutOfChainsError()
                chainmap[c] = o.id
                o.id = c
    return chainmap

def clean_pdb(inpdb, outpdb):
    with open(outpdb, 'w') as fw:
        for line in open(inpdb):
            if line.startswith('ATOM') or line.startswith('TER'):
                fw.write(line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-t", "--targetlist", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-p", "--pdblist", help="output folder for the pdb files", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()

    os.chdir(args.outdir)

    print('Downloading pdb files')

    for targetname, pdbcode in zip(open(args.targetlist), open(args.pdblist)):
        # pdb_file = name.lower() + '.pdb1'
        pdb_file = pdbcode.strip() + '.pdb1'
        if not os.path.exists(pdb_file):
            os.system(f'wget https://files.rcsb.org/download/{pdb_file}')

        if os.path.exists(pdb_file):
            clean_pdb(pdb_file, targetname.strip() + '.pdb')
        else:
            ciffile = f"{pdbcode.strip()}.cif"
            pdbfile = f"{targetname.strip()}.pdb"
            if not os.path.exists(ciffile):
                os.system(f'wget https://files.rcsb.org/download/{ciffile}')

            strucid = ciffile[:4] if len(ciffile) > 4 else "1xxx"

            # Read file
            parser = MMCIFParser()
            structure = parser.get_structure(strucid, ciffile)

            # rename long chains
            try:
                chainmap = rename_chains(structure)
            except OutOfChainsError:
                logging.error("Too many chains to represent in PDB format")
                sys.exit(1)

            for new, old in chainmap.items():
                if new != old:
                    logging.info("Renaming chain {0} to {1}".format(old, new))

            # Write PDB
            io = PDBIO()
            io.set_structure(structure)
            # TODO What happens with large structures?
            io.save(pdbfile)
