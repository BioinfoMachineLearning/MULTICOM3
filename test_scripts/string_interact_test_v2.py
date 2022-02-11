import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer
from bml_casp15.common.util import is_dir, is_file, read_option_file, makedir_if_not_exists
from bml_casp15.complex_alignment_generation import string_interact_v2
from bml_casp15.monomer_alignment_generation.alignment import *
from bml_casp15.complex_alignment_generation.pipeline_v2 import write_concatenated_alignment
import pathlib

def run_hhblits(inparams):
    fasta, outdir, hhblits_binary, database = inparams

    hhblits_runner = hhblits.HHBlits(binary_path=hhblits_binary, databases=[database])

    outfile = outdir + '/' + pathlib.Path(fasta).stem + '.a3m'

    if os.path.exists(outfile):
        return {'a3m': outfile}

    return hhblits_runner.query(fasta, outfile)

def run_jackhmmer(inparams):

    fasta, outdir, jackhmmer_binary, database = inparams

    jackhmmer_runner = jackhmmer.Jackhmmer(binary_path=jackhmmer_binary, database_path=database)

    outfile = outdir + '/' + pathlib.Path(fasta).stem + '.sto'

    if os.path.exists(outfile):
        return {'sto': outfile}

    return jackhmmer_runner.query(fasta, outfile)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)
    parser.add_argument('--fasta1', type=is_file, required=True)
    parser.add_argument('--fasta2', type=is_file, required=True)
    parser.add_argument('--hhblits', type=is_file, required=True)
    parser.add_argument('--jackhmmer', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_dir, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    # test hhblits
    outdir = args.outdir + '/hhblits'
    makedir_if_not_exists(outdir)

    process_list = []
    process_list.append([args.fasta1, outdir, args.hhblits, params['uniref_db_dir'] + '/' + params['uniref_db']])
    process_list.append([args.fasta2, outdir, args.hhblits, params['uniref_db_dir'] + '/' + params['uniref_db']])
    pool = Pool(processes=2)
    results = pool.map(run_hhblits, process_list)
    pool.close()
    pool.join()

    with open(results[0]['a3m']) as f:
        aln_1 = Alignment.from_file(f, format="a3m")

    with open(results[1]['a3m']) as f:
        aln_2 = Alignment.from_file(f, format="a3m")

    string_interact_runner = string_interact_v2.STRING_interact_v2(params['string2uniprot_map'])
    string_interact_runner.load_data()
    pair_ids = string_interact_runner.get_interactions([aln_1, aln_2])

    sequences_full, sequences_monomers, pair_ids = write_concatenated_alignment(pair_ids, [aln_1, aln_2])

    pair_ids.to_csv(f"{outdir}/string_interact.csv", index=False)
    print(pair_ids)

    complex_ailgnment_file = f"{outdir}/string_interact.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)

    # save the alignment files
    for monomer_id in sequences_monomers:
        mon_alignment_file = f"{outdir}/{monomer_id}_con.a3m"
        with open(mon_alignment_file, "w") as of:
            write_a3m(sequences_monomers[monomer_id], of)

    # test jackhmmer
    outdir = args.outdir + '/jackhmmer'
    makedir_if_not_exists(outdir)

    process_list = []
    process_list.append([args.fasta1, outdir, args.jackhmmer, params['uniref90_fasta']])
    process_list.append([args.fasta2, outdir, args.jackhmmer, params['uniref90_fasta']])
    pool = Pool(processes=2)
    results = pool.map(run_jackhmmer, process_list)
    pool.close()
    pool.join()

    with open(results[0]['sto']) as f:
        aln_1 = Alignment.from_file(f, format="stockholm")

    with open(results[1]['sto']) as f:
        aln_2 = Alignment.from_file(f, format="stockholm")

    pair_ids = string_interact_runner.get_interactions([aln_1, aln_2])

    sequences_full, sequences_monomers, pair_ids = write_concatenated_alignment(pair_ids, [aln_1, aln_2])

    pair_ids.to_csv(f"{outdir}/string_interact.csv", index=False)
    print(pair_ids)

    complex_ailgnment_file = f"{outdir}/string_interact.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)

    # save the alignment files
    for monomer_id in sequences_monomers:
        mon_alignment_file = f"{outdir}/{monomer_id}_con.a3m"
        with open(mon_alignment_file, "w") as of:
            write_a3m(sequences_monomers[monomer_id], of)

