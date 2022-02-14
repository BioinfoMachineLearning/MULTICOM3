import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer
from bml_casp15.common.util import is_dir, is_file, read_option_file, makedir_if_not_exists
from bml_casp15.complex_alignment_generation.uniclust_oxmatch_v2 import UNICLUST_oxmatch_v2
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
    parser.add_argument('--fastas', type=str, required=True)
    parser.add_argument('--hhblits', type=is_file, required=True)
    parser.add_argument('--jackhmmer', type=is_file, required=True)
    parser.add_argument('--outdir', type=is_dir, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    # test hhblits
    outdir = args.outdir + '/hhblits'
    makedir_if_not_exists(outdir)

    fasta_paths = args.fastas.split(',')

    process_list = []
    for fasta_path in fasta_paths:
        process_list.append([fasta_path, outdir, args.hhblits, params['uniclust_db_dir'] + '/' + params['uniclust_db']])

    pool = Pool(processes=5)
    results = pool.map(run_hhblits, process_list)
    pool.close()
    pool.join()

    alignments = []
    for result in results:
        with open(result['a3m']) as f:
            alignments += [Alignment.from_file(f, format="a3m")]

    pair_ids = UNICLUST_oxmatch_v2.get_interactions(alignments)

    sequences_full, sequences_monomers, pair_ids = write_concatenated_alignment(pair_ids, alignments)

    pair_ids.to_csv(f"{outdir}/oxmatch_interact.csv", index=False)
    print(pair_ids)

    complex_ailgnment_file = f"{outdir}/oxmatch_interact.a3m"
    with open(complex_ailgnment_file, "w") as of:
        write_a3m(sequences_full, of)

    # save the alignment files
    for monomer_id in sequences_monomers:
        mon_alignment_file = f"{outdir}/{monomer_id}_con.a3m"
        with open(mon_alignment_file, "w") as of:
            write_a3m(sequences_monomers[monomer_id], of)