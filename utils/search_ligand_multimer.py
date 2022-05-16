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
    foldseek_pdb_database = params['foldseek_pdb_database']
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
        pdb = f'{FLAGS.indir}/qa{i + 1}.pdb'
        if not os.path.exists(pdb):
            raise Exception(f"Cannot find {pdb}")

        current_work_dir = f"{FLAGS.outdir}/qa{i + 1}.pdb"
        makedir_if_not_exists(current_work_dir)

        chain_pdbs = split_pdb(pdb, current_work_dir)
        print(chain_pdbs)

        template_results = []

        out_template_dir = current_work_dir + '/templates'
        makedir_if_not_exists(out_template_dir)

        for chain_id in chain_pdbs:
            if chain_id not in chain_id_map:
                raise Exception("Multimer fasta file and model doesn't match!")
            monomer_work_dir = current_work_dir + '/' + chain_id_map[chain_id].description
            makedir_if_not_exists(monomer_work_dir)
            os.system(f"mv {chain_pdbs[chain_id]} {monomer_work_dir}/{chain_id_map[chain_id].description}.pdb")
            foldseek_res = search_templates(params, f"{monomer_work_dir}/{chain_id_map[chain_id].description}.pdb",
                                            monomer_work_dir + '/foldseek')
            template_results += [foldseek_res]

            monomer_template_dir = monomer_work_dir + '/templates'
            makedir_if_not_exists(monomer_template_dir)
            for j in range(len(foldseek_res['all_alignment'])):
                template = foldseek_res['all_alignment'].loc[j, 'target']
                raw_pdb = f"{params['foldseek_pdb_raw_database_dir']}/{template[0:4].lower()}.pdb1.gz"
                # print(os.path.exists(raw_pdb))
                # if not os.path.exists(raw_pdb):
                #     print(f"Cannot find {raw_pdb}")
                #     continue
                os.system(f"cp {raw_pdb} {monomer_template_dir}")

        for restype in ['local_alignment', 'global_alignment']:
            prev_df = None
            for chain_idx in range(len(template_results)):
                curr_df = template_results[chain_idx][restype]
                print(curr_df)

                pdbcodes = []
                for j in range(len(curr_df)):
                    pdbcodes += [curr_df.loc[j, 'target'][0:4]]

                curr_df = curr_df.add_suffix(str(chain_idx + 1))
                curr_df['pdbcode'] = pdbcodes

                if prev_df is None:
                    prev_df = curr_df
                else:
                    prev_df = prev_df.merge(curr_df, how="inner", on='pdbcode')

                print(prev_df)

            if restype == 'local_alignment':
                min_evalues = []
                for j in range(len(prev_df)):
                    min_evalue = np.min(np.array([prev_df.loc[j, f"evalue{k + 1}"] for k in range(len(chain_id_map))]))
                    min_evalues += [min_evalue]
                prev_df['min_evalue'] = min_evalues
                prev_df = prev_df.sort_values(by='min_evalue')
                prev_df.reset_index(inplace=True, drop=True)
                prev_df.head(1000).to_csv(f'{current_work_dir}/complex_templates_evalue.csv')
            else:
                max_tmscores = []
                for j in range(len(prev_df)):
                    max_tmscore = np.max(np.array([prev_df.loc[j, f"evalue{j + 1}"] for k in range(len(chain_id_map))]))
                    max_tmscores += [max_tmscore]
                prev_df['max_tmscore'] = max_tmscores
                prev_df = prev_df.sort_values(by=['max_tmscore'], ascending=False)
                prev_df.reset_index(inplace=True, drop=True)
                prev_df.head(1000).to_csv(f'{current_work_dir}/complex_templates_tmscore.csv')

            for j in range(len(prev_df)):
                template = prev_df.loc[j, 'pdbcode']
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
