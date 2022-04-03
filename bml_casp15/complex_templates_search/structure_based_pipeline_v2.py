import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
import pathlib
from bml_casp15.tool.foldseek import *


def search_templates_foldseek(foldseek_program, databases, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_runner = Foldseek(binary_path=foldseek_program, databases=databases)
    return foldseek_runner.query_with_tmalign(pdb=inpdb, outdir=outdir)


def create_template_df(templates):
    row_list = []
    for i in range(len(templates)):
        target = templates.loc[i, 'target']
        qaln = templates.loc[i, 'qaln']
        qstart = int(templates.loc[i, 'qstart'])
        qend = int(templates.loc[i, 'qend'])
        taln = templates.loc[i, 'taln']
        tstart = templates.loc[i, 'tstart']
        tend = templates.loc[i, 'tend']
        evalue = float(templates.loc[i, 'evalue'])
        aln_len = int(templates.loc[i, 'alnlen'])
        if target.find('.atom.gz') > 0:
            target = target.replace('.atom.gz', '')
            pdbcode = target[0:4]
        row_dict = dict(template=target.split()[0],
                        tpdbcode=pdbcode,
                        aln_temp=taln,
                        tstart=tstart,
                        tend=tend,
                        aln_query=qaln,
                        qstart=qstart,
                        qend=qend,
                        tmscore=evalue,
                        aligned_length=aln_len)
        row_list += [row_dict]
    if len(row_list) == 0:
        return pd.DataFrame(columns=['template', 'tpdbcode', 'aln_temp', 'tstart', 'tend',
                                     'aln_query', 'qstart', 'qend', 'evalue', 'aligned_length'])
    return pd.DataFrame(row_list)


class Complex_structure_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

    def search(self, monomers, outdir, is_homodimer=False):

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        prev_df = None
        for i, monomer in enumerate(monomers):
            monomer_work_dir = outdir + '/' + pathlib.Path(monomer).stem
            makedir_if_not_exists(monomer_work_dir)
            foldseek_df = search_templates_foldseek(
                foldseek_program=self.params['foldseek_program'],
                databases=[self.params['foldseek_pdb_database']],
                inpdb=monomer,
                outdir=monomer_work_dir + '/foldseek')
            curr_df = create_template_df(foldseek_df)
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, how="inner", on='tpdbcode', suffixes=(str(i), str(i + 1)))

        avg_tmscores = []
        for i in range(len(prev_df)):
            avg_tmscore = 0
            for j in range(len(monomers)):
                avg_tmscore += float(prev_df.loc[i, f"tmscore{j + 1}"])
            avg_tmscore = avg_tmscore / len(monomers)
            avg_tmscores += [avg_tmscore]

        prev_df['avg_tmscore'] = avg_tmscores

        prev_df_sorted = prev_df.sort_values(by='avg_tmscore', ascending=False)

        prev_df_sorted.head(int(self.params['template_count'])).to_csv(outdir + "/structure_templates.csv", index=False)

        print("The structure based template searching for dimers has finished!")

        os.system(f"rm -rf {outdir}/tmscore")
