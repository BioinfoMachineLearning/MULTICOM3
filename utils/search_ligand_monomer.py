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
from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('indir', None, 'Path to multimer fastas')
flags.DEFINE_string('outdir', None, 'Output directory')
FLAGS = flags.FLAGS


def search_templates(params, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_program = params['foldseek_program']
    foldseek_pdb_database = params['foldseek_rcsb_pdb_database']
    foldseek_runner = Foldseek(binary_path=foldseek_program,
                               databases=[foldseek_pdb_database])
    return foldseek_runner.query(pdb=inpdb, outdir=outdir, tmscore_threshold=0.5)


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.outdir)

    for i in range(5):
        pdb = f'{FLAGS.indir}/qa{i+1}.pdb'
        if not os.path.exists(pdb):
            raise Exception(f"Cannot find {pdb}")

        pdb_outdir = f"{FLAGS.outdir}/qa{i+1}.pdb"
        makedir_if_not_exists(pdb_outdir)
        template_dir = pdb_outdir + '/templates'
        makedir_if_not_exists(template_dir)

        foldseek_res = search_templates(params, pdb, pdb_outdir)

        for j in range(len(foldseek_res['all_alignment'])):
            template = foldseek_res['all_alignment'].loc[j, 'target']
            raw_pdb = f"{params['foldseek_rcsb_pdb_raw_database_dir']}/pdb{template[0:4].lower()}.ent.gz"
            if not os.path.exists(raw_pdb):
                print(f"Cannot find {raw_pdb}")
                continue
            os.system(f"cp {raw_pdb} {template_dir}")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'indir',
        'outdir'
    ])
    app.run(main)