import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir


def reindex_pdb(inpdb, outpdb):
    resCounter = 0
    atomCounter = 0
    prevrNum = "XX"
    prevchain = "XX"
    contents = []
    for line in open(inpdb):
        if not line.startswith('ATOM'):
            continue
        rnum = line[22:27]
        chain = line[21:22]
        add_ter = False
        if prevchain != chain:
            if prevchain != "XX":
                contents += ["TER\n"]
            prevchain = chain
            resCounter = 1
            atomCounter = 0
            prevrNum = rnum
        elif prevrNum != rnum:
            prevrNum = rnum
            resCounter += 1
        atomCounter += 1
        rnum_string = "{:>4}".format(resCounter)
        anum_string = "{:>5}".format(atomCounter)
        if chain == "0":
            row_chain = "a"
        else:
            row_chain = chain
        row = f"{line[:6]}{anum_string}{line[11:21]}{row_chain}{rnum_string}{line[26:]}"
        contents += [row]
    contents += ["TER\n"]
    with open(outpdb, 'w') as fw:
        fw.writelines("".join(contents))


add_local_script = '/home/bml_casp15/TS_run/scripts/add_local_score.pl'
pdb2seq = '/home/bml_casp15/BML_CASP15/utils/pdb2seq.pl'

def add_local_scores(inpdb, outdir, outpdb):
    tmpdir = outdir + '/tmp'
    chain_contents = {}
    for line in open(inpdb):
        if line.startswith('ATOM'):
            chain_id = line[21]
            if chain_id not in chain_contents:
                chain_contents[chain_id] = [line]
            else:
                chain_contents[chain_id] += [line]

    new_chain_contents = {}
    for chain_id in chain_contents:
        chain_workdir = tmpdir + '/' + chain_id
        makedir_if_not_exists(chain_workdir)
        makedir_if_not_exists(chain_workdir + '/in')
        with open(f"{chain_workdir}/in/{chain_id}.pdb", 'w') as fw:
            fw.writelines(chain_contents[chain_id])
        cmd = f"perl {pdb2seq} {chain_workdir}/in/{chain_id}.pdb"
        print(cmd)
        seq = os.popen(f"perl {pdb2seq} {chain_workdir}/in/{chain_id}.pdb").read()[0].strip()
        with open(chain_workdir + '/' + chain_id + '.fasta', 'w') as fw:
            fw.write(f">{chain_id}\n{seq}")
        makedir_if_not_exists(chain_workdir + '/out')
        os.system(f"perl {add_local_script} {chain_id} {chain_workdir}/{chain_id}.fasta {chain_workdir}/in {chain_workdir}/out "
                  f"/home/bml_casp15/TS_run/scripts")
        new_chain_contents[chain_id] = open(f"{chain_workdir}/out/{chain_id}.pdb").readlines()

    with open(outpdb, 'w') as fw:
        for chain_id in new_chain_contents:
            fw.writelines(new_chain_contents[chain_id])
            fw.write("TER\n")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inpdb', type=is_file, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    args = parser.parse_args()
    makedir_if_not_exists(args.outdir)
    reindex_pdb(args.inpdb, args.outdir + '/final_pre_score.pdb')
    add_local_scores(args.outdir + '/final_pre_score.pdb', args.outdir, args.outdir + '/final.pdb')
