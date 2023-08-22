import copy
import os
import sys
import time
from multicom3.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
import pathlib
from multicom3.tool.foldseek import *
from multicom3.monomer_templates_concatenation.sequence_based_pipeline import assess_hhsearch_hit, PrefilterError
from multicom3.monomer_templates_concatenation.parsers import TemplateHit
from multicom3.monomer_structure_refinement.util import *


def search_templates_foldseek(foldseek_program, databases, inpdb, outdir):
    makedir_if_not_exists(outdir)
    foldseek_runner = Foldseek(binary_path=foldseek_program, databases=databases)
    return foldseek_runner.query_with_tmalign(pdb=inpdb, outdir=outdir,maxseq=300)


def create_template_df(templates, query_sequence, check_hit=True):
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

        if check_hit:
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
        monomer_template_results = []
        for i, monomer in enumerate(monomers_pdbs):
            monomer_work_dir = outdir + '/' + pathlib.Path(monomer).stem
            makedir_if_not_exists(monomer_work_dir)
            foldseek_df = search_templates_foldseek(
                foldseek_program=self.params['foldseek_program'],
                databases=[self.params['foldseek_pdb_database']],
                inpdb=monomer,
                outdir=monomer_work_dir + '/foldseek')
            curr_df = create_template_df(foldseek_df, monomer_sequences[i])
            monomer_template_results += [curr_df]

            curr_df = curr_df.add_suffix(f"{i + 1}")
            curr_df['tpdbcode'] = curr_df[f'tpdbcode{i + 1}']
            curr_df = curr_df.drop([f'tpdbcode{i + 1}'], axis=1)

            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, how="inner", on='tpdbcode')

        avg_tmscores = []
        for i in range(len(prev_df)):
            avg_tmscore = 0
            for j in range(len(monomers_pdbs)):
                avg_tmscore += float(prev_df.loc[i, f"tmscore{j + 1}"])
            avg_tmscore = avg_tmscore / len(monomers_pdbs)
            avg_tmscores += [avg_tmscore]

        prev_df['avg_tmscore'] = avg_tmscores

        prev_df_sorted = prev_df.sort_values(by=['tpdbcode','avg_tmscore'], ascending=False)
        prev_df_sorted.reset_index(inplace=True, drop=True)

        keep_indices = []
        pdbcodes = []
        cover_chains_in_pdb = {}
        for i in range(len(prev_df_sorted)):
            chain_count = len(set([prev_df_sorted.loc[i, f'template{j + 1}'] for j in range(len(monomers_pdbs))]))
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
        prev_df_sorted_filtered = prev_df_sorted_filtered.sort_values(by='avg_tmscore', ascending=False)
        prev_df_sorted_filtered.reset_index(inplace=True, drop=True)

        if len(keep_indices) < 50:
            print(f"template count is smaller than 50, add monomer templates")
            prev_pd = None
            prev_pd_v2 = None
            for i in range(len(monomer_template_results)):
                row_list = []
                row_index = 0
                for j in range(len(prev_df_sorted_filtered)):
                    row_dict = dict(index=row_index,
                                    template=prev_df_sorted_filtered.loc[j, f'template{i + 1}'],
                                    aln_query=prev_df_sorted_filtered.loc[j, f'aln_query{i + 1}'],
                                    qstart=prev_df_sorted_filtered.loc[j, f'qstart{i + 1}'],
                                    qend=prev_df_sorted_filtered.loc[j, f'qend{i + 1}'],
                                    aln_temp=prev_df_sorted_filtered.loc[j, f'aln_temp{i + 1}'],
                                    tstart=prev_df_sorted_filtered.loc[j, f'tstart{i + 1}'],
                                    tend=prev_df_sorted_filtered.loc[j, f'tend{i + 1}'],
                                    tmscore=prev_df_sorted_filtered.loc[j, f'tmscore{i + 1}'],
                                    aligned_length=prev_df_sorted_filtered.loc[j, f'aligned_length{i + 1}'])
                    row_list += [row_dict]
                    row_index += 1

                seen_templates_sequences = [
                    f"{prev_df_sorted_filtered.loc[j, f'template{i + 1}']}_" \
                    f"{prev_df_sorted_filtered.loc[j, f'aln_temp{i + 1}']}" for j in range(len(prev_df_sorted_filtered))]

                for j in range(len(monomer_template_results[i])):
                    if row_index > 50:
                        break
                    if f"{monomer_template_results[i].loc[j, 'template']}_{monomer_template_results[i].loc[j, 'aln_temp']}" \
                            not in seen_templates_sequences:

                        hit = TemplateHit(index=j,
                                          name=monomer_template_results[i].loc[j, f'template'].split('.')[0],
                                          aligned_cols=int(monomer_template_results[i].loc[j, 'aligned_length']),
                                          query=monomer_template_results[i].loc[j, 'aln_query'],
                                          hit_sequence=monomer_template_results[i].loc[j, 'aln_temp'],
                                          indices_query=build_alignment_indices(
                                              monomer_template_results[i].loc[j, 'aln_query'],
                                              monomer_template_results[i].loc[j, 'qstart']),
                                          indices_hit=build_alignment_indices(
                                              monomer_template_results[i].loc[j, 'aln_temp'],
                                              monomer_template_results[i].loc[j, 'tstart']),
                                          sum_probs=0.0)

                        try:
                            assess_hhsearch_hit(hit=hit, query_sequence=monomer_sequences[i])
                        except PrefilterError as e:
                            msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                            print(msg)
                            continue

                        row_dict = dict(index=row_index,
                                        template=monomer_template_results[i].loc[j, f'template'],
                                        aln_query=monomer_template_results[i].loc[j, f'aln_query'],
                                        qstart=monomer_template_results[i].loc[j, f'qstart'],
                                        qend=monomer_template_results[i].loc[j, f'qend'],
                                        aln_temp=monomer_template_results[i].loc[j, f'aln_temp'],
                                        tstart=monomer_template_results[i].loc[j, f'tstart'],
                                        tend=monomer_template_results[i].loc[j, f'tend'],
                                        tmscore=monomer_template_results[i].loc[j, f'tmscore'],
                                        aligned_length=monomer_template_results[i].loc[j, f'aligned_length'])
                        row_list += [row_dict]
                        row_index += 1

                if len(row_list) == 0:
                    curr_pd = pd.DataFrame(columns=['index', 'template', 'aln_query', 'qstart',
                                                    'qend', 'aln_temp', 'tstart', 'tend', 'tmscore', 'aligned_length'])
                else:
                    curr_pd = pd.DataFrame(row_list)
                curr_pd = curr_pd.add_suffix(f"{i + 1}")
                curr_pd['index'] = curr_pd[f'index{i + 1}']
                curr_pd = curr_pd.drop([f'index{i + 1}'], axis=1)
                if prev_pd is None:
                    prev_pd = curr_pd
                    prev_pd_v2 = curr_pd
                else:
                    prev_pd = prev_pd.merge(curr_pd, how="inner", on='index')
                    prev_pd_v2 = prev_pd_v2.merge(curr_pd, how="outer", on='index')

            prev_pd.to_csv(outdir + "/structure_templates.csv", index=False)
            prev_pd_v2.to_csv(outdir + "/structure_templates_v2.csv", index=False)
        else:
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

        if os.path.exists(outdir + "/structure_templates_v2.csv"):
            prev_df_sorted_filtered = pd.read_csv(outdir + "/structure_templates_v2.csv")
            for i in range(len(prev_df_sorted_filtered)):
                for j in range(len(monomers_pdbs)):
                    template_name = prev_df_sorted_filtered.loc[i, f'template{j + 1}']
                    if not pd.isna(template_name):
                        template_pdb = template_name.split()[0]
                        os.system(f"cp {self.params['foldseek_pdb_database_dir']}/{template_pdb}.atom.gz .")
                        os.system(f"gunzip -f {template_pdb}.atom.gz")
