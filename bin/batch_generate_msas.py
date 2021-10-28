import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer
from bml_casp15.common.util import is_dir, is_file, read_option_file, makedir_if_not_exists
from bml_casp15.alignment.alignment import *

def run_hhblits(inparams):

    fasta, outdir, hhblits_binary, database = inparams

    hhblits_runner = hhblits.HHBlits(binary_path=hhblits_binary, databases=[database])

    return hhblits_runner.query(fasta, outdir)


def run_jackhmmer(inparams):

    fasta, outdir, jackhmmer_binary, database = inparams

    jackhmmer_runner = jackhmmer.Jackhmmer(binary_path=jackhmmer_binary, database_path=database)

    return jackhmmer_runner.query(fasta, outdir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--heterodimer_fasta', type=is_file, required=True)
    parser.add_argument('--hhblits', type=is_file, required=True)
    parser.add_argument('--jackhmmer', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_dir, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    # Three kinds of msas should be generated
    # 1. HHblits on UniRef30
    # 2. Jackhmmer on UniRef30
    # 3. HHblits on UniClust30


    fastadir = args.outdir + '/fasta'
    makedir_if_not_exists(fastadir)

    contents = open(args.heterodimer_fasta).readlines()

    targets = []

    while len(contents) > 0:

        name = contents.pop(0)
        name = name.replace('\n', '').replace('\r', '')

        seq = contents.pop(0)
        seq = seq.replace('\n', '').replace('\r', '')

        if name.find('>') != 0:
            die("fasta file is not in fasta format.\n")

        name = name[1:]
        if name.find('.') > 0:
            die("sequence name can't include .\n")

        with open(fastadir + '/' + name + '.fasta', "w") as f:
            f.write(f">{name}\n{seq}")

        targets += [name]

    # 1. HHblits on UniRef30
    outdir = args.outdir + '/hhblits_uniref30'
    makedir_if_not_exists(outdir)

    process_list = []

    for target in targets:

        process_list.append([fastadir + '/' + target + '.fasta', outdir,
                             args.hhblits, params['uniref_db_dir'] + '/' + params['uniref_db']])

    pool = Pool(processes=7)
    results = pool.map(run_hhblits, process_list)
    pool.close()
    pool.join()

    # 2. Jackhmmer on UniRef30

    outdir = args.outdir + '/jackhmmer_uniref30'
    makedir_if_not_exists(outdir)

    process_list = []

    for target in targets:

        process_list.append([fastadir + '/' + target + '.fasta', outdir, args.jackhmmer, params['uniref90_fasta']])

    pool = Pool(processes=7)
    results = pool.map(run_jackhmmer, process_list)
    pool.close()
    pool.join()


    # 3. HHblits on UniClust30
    outdir = args.outdir + '/hhblits_uniclust30'
    makedir_if_not_exists(outdir)

    process_list = []

    for target in targets:

        process_list.append([fastadir + '/' + target + '.fasta', outdir,
                             args.hhblits, params['uniclust_db_dir'] + '/' + params['uniclust_db']])

    pool = Pool(processes=7)
    results = pool.map(run_hhblits, process_list)
    pool.close()
    pool.join()

