import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.tool.foldseek import *
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.quaternary_structure_refinement.util import split_pdb
from bml_casp15.common.protein import parse_fasta, make_chain_id_map
from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('indir', None, 'Path to multimer fastas')
flags.DEFINE_string('outdir', None, 'Output directory')
FLAGS = flags.FLAGS

def search_templates(params, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_program = params['foldseek_program']
    foldseek_pdb_database = params['foldseek_pdb_database_dir']
    foldseek_runner = Foldseek(binary_path=foldseek_program,
                               databases=[foldseek_pdb_database])
    return foldseek_runner.query(pdb=inpdb, outdir=outdir, progressive_threshold=2000, tmscore_threshold=0.5)


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.outdir)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    for i in range(5):
        pdb = f'{FLAGS.indir}/qa{i+1}.pdb'
        if not os.path.exists(pdb):
            raise Exception(f"Cannot find {pdb}")

        current_work_dir = f"{FLAGS.outdir}/qa{i+1}.pdb"
        makedir_if_not_exists(current_work_dir)

        chain_pdbs = split_pdb(start_pdb, current_work_dir)

        template_results = []

        out_template_dir = current_work_dir + '/templates'
        makedir_if_not_exists(out_template_dir)

        for chain_id in chain_pdbs:
            if chain_id not in chain_id_map:
                raise Exception("Multimer fasta file and model doesn't match!")
            monomer_work_dir = current_work_dir + '/' + chain_id_map[chain_id].description
            makedir_if_not_exists(monomer_work_dir)
            os.system(f"mv {chain_pdbs[chain_id]} {monomer_work_dir}/{chain_id_map[chain_id].description}.pdb")
            foldseek_res = search_templates(params, pdb, monomer_work_dir + '/foldseek')
            template_results[chain_id] = foldseek_res

        for restype in ['local_alignment', 'global_alignment']:
            prev_df = None
            for chain_idx, (chain_id, template_result) in enumerate(zip(chain_id_map, template_results)):
                chain_template_res = template_result[restype]

                pdbcodes = []
                for i in range(chain_template_res):
                    pdbcodes += [chain_template_res.loc[i, 'target'][0:4]]
                chain_template_res['pdbcode'] = pdbcodes

                if prev_df is None:
                    prev_df = chain_template_res
                else:
                    prev_df = prev_df.merge(chain_template_res, how="inner", on='pdbcode',
                                            suffixes=(str(chain_idx), str(chain_idx + 1)))
            prev_df.to_csv(f'{current_work_dir}/complex_templates_evalue.csv')
            for i in range(len(prev_df)):
                template = prev_df.loc[i, 'pdbcode']
                raw_pdb = f"{params['foldseek_pdb_raw_database_dir']}/{template[0:4].lower()}.pdb1.gz"
                if not os.path.exists(raw_pdb):
                    print(f"Cannot find {raw_pdb}")
                    continue
                os.system(f"cp {raw_pdb} {out_template_dir}")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'indir',
        'outdir'
    ])
    app.run(main)