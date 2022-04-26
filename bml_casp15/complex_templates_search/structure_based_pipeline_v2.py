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
from bml_casp15.complex_templates_search.sequence_based_pipeline import assess_hhsearch_hit, PrefilterError
from bml_casp15.complex_templates_search.parsers import TemplateHit
from bml_casp15.monomer_structure_refinement.util import *


def search_templates_foldseek(foldseek_program, databases, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_runner = Foldseek(binary_path=foldseek_program, databases=databases)
    return foldseek_runner.query_with_tmalign(pdb=inpdb, outdir=outdir)


def create_template_df(templates, query_sequence):
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

        hit = TemplateHit(index=i,
                          name=templates.loc[i, 'target'].split('.')[0],
                          aligned_cols=int(templates.loc[i, 'alnlen']),
                          query=templates.loc[i, 'qaln'],
                          hit_sequence=templates.loc[i, 'taln'],
                          indices_query=build_alignment_indices(templates.loc[i, 'qaln'],
                                                                templates.loc[i, 'qstart']),
                          indices_hit=build_alignment_indices(templates.loc[i, 'taln'],
                                                              templates.loc[i, 'tstart']),
                          sum_probs=0.0)

        try:
            assess_hhsearch_hit(hit=hit, query_sequence=query_sequence)
        except PrefilterError as e:
            msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
            print(msg)
            continue

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

    def search(self, monomer_sequences, monomers_pdbs, outdir, is_homodimer=False):

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        prev_df = None
        for i, monomer in enumerate(monomers_pdbs):
            monomer_work_dir = outdir + '/' + pathlib.Path(monomer).stem
            makedir_if_not_exists(monomer_work_dir)
            foldseek_df = search_templates_foldseek(
                foldseek_program=self.params['foldseek_program'],
                databases=[self.params['foldseek_pdb_database']],
                inpdb=monomer,
                outdir=monomer_work_dir + '/foldseek')
            curr_df = create_template_df(foldseek_df, monomer_sequences[i])
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, how="inner", on='tpdbcode', suffixes=(str(i), str(i + 1)))

        avg_tmscores = []
        for i in range(len(prev_df)):
            avg_tmscore = 0
            for j in range(len(monomers_pdbs)):
                avg_tmscore += float(prev_df.loc[i, f"tmscore{j + 1}"])
            avg_tmscore = avg_tmscore / len(monomers_pdbs)
            avg_tmscores += [avg_tmscore]

        prev_df['avg_tmscore'] = avg_tmscores
        prev_df_sorted = prev_df.sort_values(by='avg_tmscore', ascending=False)
        prev_df_sorted.reset_index(inplace=True, drop=True)

        keep_indices = []
        pdbcodes = []
        cover_chains_in_pdb = {}
        for i in range(len(prev_df_sorted)):
            chain_count = len(set([prev_df_sorted.loc[i, f'template{j+1}'] for j in range(len(monomers_pdbs))]))
            if prev_df_sorted.loc[i, 'tpdbcode'] not in pdbcodes:
                if len(pdbcodes) > 0:
                    max_index = -1
                    max_count = 0
                    for index in cover_chains_in_pdb:
                        if cover_chains_in_pdb[index] > max_count:
                            max_index = index
                            max_count = cover_chains_in_pdb[index]
                    keep_indices += [max_index]

                pdbcodes += [prev_df_sorted.loc[i, 'tpdbcode']]
                cover_chains_in_pdb = {i: chain_count}
            else:
                cover_chains_in_pdb[i] = chain_count

        if len(cover_chains_in_pdb) > 0:
            max_index = -1
            max_count = 0
            for index in cover_chains_in_pdb:
                if cover_chains_in_pdb[index] > max_count:
                    max_index = index
                    max_count = cover_chains_in_pdb[index]
            keep_indices += [max_index]

        prev_df_sorted_filtered = prev_df_sorted.iloc[keep_indices]

        prev_df_sorted_filtered.head(100).to_csv(outdir + "/structure_templates.csv", index=False)

        print("The structure based template searching for dimers has finished!")

        os.system(f"rm -rf {outdir}/tmscore")

        prev_df_sorted_filtered = pd.read_csv(outdir + "/structure_templates.csv")
        cwd = os.getcwd()
        template_dir = outdir + '/templates'
        makedir_if_not_exists(template_dir)
        os.chdir(template_dir)
        for i in range(len(prev_df_sorted_filtered)):
            for j in range(len(monomers_pdbs)):
                template_pdb = prev_df_sorted_filtered.loc[i, f'template{j + 1}'].split()[0]
                os.system(f"cp {self.params['foldseek_pdb_database_dir']}/{template_pdb}.atom.gz .")
                os.system(f"gunzip -f {template_pdb}.atom.gz")
        os.chdir(cwd)
