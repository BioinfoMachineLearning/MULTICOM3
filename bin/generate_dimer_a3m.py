import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.pipeline.concatenate import Concatenate
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer

def generate_a3ms_for_single_seq(inparams):

    name, fasta, outdir, params = inparams
    uniref30 = params['uniref_db_dir'] + '/' + params['uniref_db']
    uniclust30 = params['uniclust_db_dir'] + '/' + params['uniclust_db']
    uniref90_fasta = params['uniref90_fasta']
    hhblits_binary = params['hhblits_program']
    jackhmmer_binary = params['jackhmmer_program']

    os.chdir(outdir)
    uniref30_a3m = f"{name}_uniref30.a3m"
    uniclust30_a3m = f"{name}_uniclust30.a3m"
    uniref90_sto = f"{name}_uniref90.sto"

    if not os.path.exists(uniref30_a3m):
        hhblits_runner1 = hhblits.HHBlits(binary_path=hhblits_binary, databases=[uniref30])
        uniref30_result = hhblits_runner1.query(fasta, outdir)
        os.system(f"mv {uniref30_result['a3m']} {name}_uniref30.a3m")

    if not os.path.exists(uniclust30_a3m):
        hhblits_runner2 = hhblits.HHBlits(binary_path=hhblits_binary, databases=[uniclust30])
        uniclust30_result = hhblits_runner2.query(fasta, outdir)
        os.system(f"mv {uniclust30_result['a3m']} {name}_uniclust30.a3m")

    if not os.path.exists(uniref90_sto):
        jackhmmer_runner = jackhmmer.Jackhmmer(binary_path=jackhmmer_binary, database_path=uniref90_fasta)
        uniref90_result = jackhmmer_runner.query(fasta, outdir)
        os.system(f"mv {uniref90_result['sto']} {name}_uniref90.sto")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--fastadir', type=is_dir, required=True)
    parser.add_argument('--dimerlist', type=is_file, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    makedir_if_not_exists(args.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    process_list = []
    print("Start to generate alignments for targets")
    alignments = []
    with open(args.dimerlist) as f:
        for line in f:
            chain1, chain2 = line.rstrip('\n').split()

            if not os.path.exists(f"{args.fastadir}/{chain1}.fasta") or not os.path.exists(f"{args.fastadir}/{chain2}.fasta"):
                print(f"Cannot find fasta file for {chain1} or {chain2}!")
                continue

            outdir = args.output_dir + f"/{chain1}_{chain2}"
            makedir_if_not_exists(outdir)

            os.system(f"cp {args.fastadir}/{chain1}.fasta {outdir}")
            os.system(f"cp {args.fastadir}/{chain2}.fasta {outdir}")

            makedir_if_not_exists(outdir + '/' + chain1)
            process_list.append([chain1, f"{outdir}/{chain1}.fasta", outdir + '/' + chain1, params])

            makedir_if_not_exists(outdir + '/' + chain2)
            process_list.append([chain2, f"{outdir}/{chain2}.fasta", outdir + '/' + chain2, params])

            chain1_a3ms = {'name': chain1,
                           'uniref_a3m': f"{outdir}/{chain1}/{chain1}_uniref30.a3m",
                           'uniref_sto': f"{outdir}/{chain1}/{chain1}_uniref90.sto",
                           'uniclust_a3m': f"{outdir}/{chain1}/{chain1}_uniclust30.a3m"}

            chain2_a3ms = {'name': chain2,
                           'uniref_a3m': f"{outdir}/{chain2}/{chain2}_uniref30.a3m",
                           'uniref_sto': f"{outdir}/{chain2}/{chain2}_uniref90.sto",
                           'uniclust_a3m': f"{outdir}/{chain2}/{chain2}_uniclust30.a3m"}

            alignment = {'chain1': chain1_a3ms, 'chain2': chain2_a3ms, 'outdir': outdir}
            alignments += [alignment]

    pool = Pool(processes=5)
    results = pool.map(generate_a3ms_for_single_seq, process_list)
    pool.close()
    pool.join()

    print("Start to concatenate alignments for dimers")

    dimer_concatenate = Concatenate(params)

    alignments = dimer_concatenate.concatenate(alignments)



