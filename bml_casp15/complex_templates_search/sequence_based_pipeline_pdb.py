import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
from bml_casp15.complex_templates_search import parsers
from bml_casp15.tool import hhsearch
from bml_casp15.tool import hhalign
import dataclasses
import numpy as np


# Prefilter exceptions.
class PrefilterError(Exception):
    """A base class for template prefilter exceptions."""


class DateError(PrefilterError):
    """An error indicating that the hit date was after the max allowed date."""


class AlignRatioError(PrefilterError):
    """An error indicating that the hit align ratio to the query was too small."""


class DuplicateError(PrefilterError):
    """An error indicating that the hit was an exact subsequence of the query."""


class LengthError(PrefilterError):
    """An error indicating that the hit was too short."""


@dataclasses.dataclass(frozen=False)
class monomer_template_input:
    name: str
    msa_path: str
    hmm_path: str
    seq: str


def create_df(hits, prev_template_hit_indices=None):
    row_list = []
    if prev_template_hit_indices is None:
        for index, hit in enumerate(hits):
            row_dict = dict(index=index,
                            template=hit.name,
                            aligned_length=hit.aligned_cols,
                            aln_temp=hit.hit_sequence,
                            tstart=hit.indices_hit[0] + 1,
                            tend=hit.indices_hit[len(hit.indices_hit) - 1] + 1,
                            aln_query=hit.query,
                            qstart=hit.indices_query[0] + 1,
                            qend=hit.indices_query[len(hit.indices_query) - 1] + 1,
                            sum_probs=hit.sum_probs)
            row_list += [row_dict]
    else:
        for index, hit in zip(prev_template_hit_indices, hits):
            row_dict = dict(index=index,
                            template=hit.name,
                            aligned_length=hit.aligned_cols,
                            aln_temp=hit.hit_sequence,
                            tstart=hit.indices_hit[0] + 1,
                            tend=hit.indices_hit[len(hit.indices_hit) - 1] + 1,
                            aln_query=hit.query,
                            qstart=hit.indices_query[0] + 1,
                            qend=hit.indices_query[len(hit.indices_query) - 1] + 1,
                            sum_probs=hit.sum_probs)
            row_list += [row_dict]

    if len(row_list) == 0:
        return pd.DataFrame(columns=['index', 'template', 'aligned_length', 'aln_temp',
                                     'tstart', 'tend', 'aln_query', 'qstart',
                                     'qend', 'sum_probs'])
    return pd.DataFrame(row_list)


def assess_hhsearch_hit(
        hit: parsers.TemplateHit,
        query_sequence: str,
        max_subsequence_ratio: float = 0.95,
        min_align_ratio: float = 0.1) -> bool:
    aligned_cols = hit.aligned_cols
    align_ratio = aligned_cols / len(query_sequence)

    template_sequence = hit.hit_sequence.replace('-', '')
    length_ratio = float(len(template_sequence)) / len(query_sequence)

    duplicate = (template_sequence in query_sequence and
                 length_ratio > max_subsequence_ratio)

    if align_ratio <= min_align_ratio:
        raise AlignRatioError('Proportion of residues aligned to query too small. '
                              f'Align ratio: {align_ratio}.')

    if duplicate:
        raise DuplicateError('Template is an exact subsequence of query with large '
                             f'coverage. Length ratio: {length_ratio}.')

    if len(template_sequence) < 10:
        raise LengthError(f'Template too short. Length: {len(template_sequence)}.')

    return True


class Complex_sequence_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

        self.template_searcher = hhsearch.HHSearch(
            binary_path='/home/bml_casp15/BML_CASP15/tools/hhsuite-3.2.0/bin/hhsearch',
            databases=[params['pdb_sort90_hhsuite_database']],
            input_format='hmm')

        self.seq_template_searcher = hhalign.HHAlign(binary_path=params['hhalign_program'])

        self.cluster_tsv = params['pdb_sort90_hhsuite_cluster_tsv']

        self.hhmake_program = params['hhmake_program']

        self.atom_dir = params['pdb_sort90_atom_dir']

        self.template_hmm_dir = params['pdb_sort90_hmm_dir']

    def find_matches_between_pdbcodes(self, monomer_code1, monomer_code2):

        cluster_members = set()
        cmd = f"grep {monomer_code2} {self.cluster_tsv}"
        grep_results = os.popen(cmd).read().rstrip('\n').split('\n')
        for grep_result in grep_results:
            pdbcode1, pdbcode2 = grep_result.split('\t')
            cluster_members.add(pdbcode1)
            cluster_members.add(pdbcode2)

        # monomer_code2 is a member in the cluster
        if len(grep_results) == 1 and pdbcode1 != monomer_code2 and pdbcode2 == monomer_code2:
            cmd = f"grep {pdbcode1} {self.cluster_tsv}"
            grep_results = os.popen(cmd).read().rstrip('\n').split('\n')
            for grep_result in grep_results:
                pdbcode1, pdbcode2 = grep_result.split()
                cluster_members.add(pdbcode1)
                cluster_members.add(pdbcode2)

        match_members = [pdbcode for pdbcode in cluster_members if pdbcode[0:4] == monomer_code1[0:4]]

        if len(match_members) == 0:
            return ""

        return match_members[0]

    def align_template(self, src_hmm, trg_hmm, outdir):

        pdb_templates_result = self.seq_template_searcher.query(src_hmm, trg_hmm, outdir)

        pdb_templates_hits = parsers.parse_hhr(hhr_string=pdb_templates_result)

        if len(pdb_templates_hits) > 1:
            print("HHalign returns more than 1 result!")

        return pdb_templates_hits[0]

    def concatenate_templates(self, monomer_inputs, monomer_template_results, outdir, check_hit=True):

        curr_template_hits = []
        for hit in monomer_template_results[0]:
            if check_hit:
                try:
                    assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[0].seq)
                except PrefilterError as e:
                    msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                    print(msg)
                    continue
            curr_template_hits += [hit]

        prev_pd = create_df(curr_template_hits)
        prev_pd = prev_pd.add_suffix("1")
        prev_pd['index'] = prev_pd['index1']
        prev_pd = prev_pd.drop(['index1'], axis=1)
        for i in range(1, len(monomer_template_results)):
            template_count = 0
            seen_templates = []
            curr_template_hits = []
            prev_template_hit_indices = []
            for j in range(len(prev_pd)):
                hit1_name = prev_pd.loc[j, f'template{i}'].split()[0]
                print(f"finding hits for {hit1_name}")
                if template_count > 100:
                    break
                for k in range(len(monomer_template_results[i])):
                    hit2_name = monomer_template_results[i][k].name.split()[0]
                    hit = None
                    if hit1_name[0:4] == hit2_name[0:4]:
                        hit = copy.deepcopy(monomer_template_results[i][k])
                    else:
                        hit_name = self.find_matches_between_pdbcodes(hit1_name, hit2_name)
                        # print(hit1_name)
                        if len(hit_name) > 0:
                            if not os.path.exists(self.template_hmm_dir + '/' + hit_name + '.hmm'):
                                print(f"cannot find {hit_name}.hmm in {self.template_hmm_dir}")
                                continue
                            hit = self.align_template(
                                monomer_inputs[i].hmm_path,
                                self.template_hmm_dir + '/' + hit_name + '.hmm',
                                outdir)

                    if hit is None:
                        continue

                    if check_hit:
                        try:
                            assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[i].seq)
                        except PrefilterError as e:
                            msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                            print(msg)
                            continue

                    if hit.name + hit.hit_sequence not in seen_templates:
                        template_count += 1
                        print(f"found template {hit.name.split()[0]} for {monomer_inputs[i].name} "
                              f"and template {hit1_name} for {monomer_inputs[i - 1].name}, "
                              f"current count: {template_count}")
                        seen_templates += [hit.name + hit.hit_sequence]
                        curr_template_hits += [hit]
                        prev_template_hit_indices += [j]
                        break

            curr_pd = create_df(curr_template_hits, prev_template_hit_indices)
            # print(curr_pd)
            curr_pd = curr_pd.add_suffix(f"{i + 1}")
            curr_pd['index'] = curr_pd[f'index{i + 1}']
            curr_pd = curr_pd.drop([f'index{i + 1}'], axis=1)
            prev_pd = prev_pd.merge(curr_pd, how="inner", on='index')

        max_probs = []
        for i in range(len(prev_pd)):
            sum_probs = []
            for j in range(len(monomer_inputs)):
                sum_prob = float(prev_pd.loc[i, f'sum_probs{j + 1}'])
                sum_probs += [sum_prob]
            max_probs += [np.max(np.array(sum_probs))]
        prev_pd['max_probs'] = max_probs
        prev_pd['tpdbcode'] = [prev_pd.loc[i, 'template1'][0:4] for i in range(len(prev_pd))]
        prev_pd = prev_pd.sort_values(by=['max_probs'], ascending=False)

        return prev_pd

    def filter_same_pdbcodes(self, indf, monomer_count):

        indf_sorted = indf.sort_values(by=['tpdbcode','max_probs'], ascending=False)
        indf_sorted.reset_index(inplace=True, drop=True)

        keep_indices = []
        pdbcodes = []
        cover_chains_in_pdb = {}
        for i in range(len(indf_sorted)):
            chain_count = len(set([indf_sorted.loc[i, f'template{j + 1}'] for j in range(monomer_count)]))
            if indf_sorted.loc[i, 'tpdbcode'] not in pdbcodes:
                if len(pdbcodes) > 0:
                    max_index = -1
                    max_count = 0
                    for index in cover_chains_in_pdb:
                        if cover_chains_in_pdb[index] > max_count:
                            max_index = index
                            max_count = cover_chains_in_pdb[index]
                    keep_indices += [max_index]

                pdbcodes += [indf_sorted.loc[i, 'tpdbcode']]
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

        indf_sorted_filtered = indf_sorted.iloc[keep_indices]
        indf_sorted_filtered = indf_sorted_filtered.sort_values(by='max_probs', ascending=False)
        indf_sorted_filtered.reset_index(inplace=True, drop=True)

        return indf_sorted_filtered

    def search(self, monomer_inputs, outdir):

        monomer_template_results = []

        monomer_names = []

        for monomer_input in monomer_inputs:

            monomer_outdir = outdir + '/' + monomer_input.name

            makedir_if_not_exists(monomer_outdir)

            pdb_hits_out_path = os.path.join(monomer_outdir, 'output.hhr')
            if os.path.exists(pdb_hits_out_path):
                pdb_templates_result = open(pdb_hits_out_path).read()
                monomer_input.msa_path = monomer_outdir + '/' + monomer_input.name + '.a3m'
                monomer_input.hmm_path = monomer_outdir + '/' + monomer_input.name + '.hmm'
            else:
                with open(monomer_input.msa_path) as f:
                    msa_for_templates = f.read()
                    if monomer_input.msa_path.find('.sto') > 0:
                        msa_for_templates = parsers.deduplicate_stockholm_msa(msa_for_templates)
                        msa_for_templates = parsers.remove_empty_columns_from_stockholm_msa(msa_for_templates)
                        msa_for_templates = parsers.convert_stockholm_to_a3m(msa_for_templates)

                        with open(monomer_outdir + '/' + monomer_input.name + '.a3m', 'w') as fw:
                            fw.write(msa_for_templates)

                        monomer_input.msa_path = monomer_outdir + '/' + monomer_input.name + '.a3m'

                    monomer_input.hmm_path = monomer_outdir + '/' + monomer_input.name + '.hmm'
                    os.system(f"{self.hhmake_program} -i {monomer_input.msa_path} -o {monomer_input.hmm_path}")
                    with open(monomer_input.hmm_path) as f:
                        msa_for_templates = f.read()
                        pdb_templates_result = self.template_searcher.query(msa_for_templates, monomer_outdir)

            pdb_template_hits = parsers.parse_hhr(hhr_string=pdb_templates_result)

            pdb_template_hits = sorted(pdb_template_hits, key=lambda x: x.sum_probs, reverse=True)

            monomer_template_results += [pdb_template_hits]

            monomer_names += [monomer_input.name]

        concatenated_pd = self.concatenate_templates(monomer_inputs, monomer_template_results, outdir)
        concatenated_pd = self.filter_same_pdbcodes(concatenated_pd, len(monomer_template_results))
        concatenated_pd_v2 = copy.deepcopy(concatenated_pd)

        if len(concatenated_pd) < 50:

            print(f"template count is smaller than 50, add monomer templates")
            prev_pd = None
            prev_pd_v2 = None
            for i in range(len(monomer_template_results)):
                seen_templates_sequences = [
                    f"{concatenated_pd.loc[j, f'template{i + 1}']}_{concatenated_pd.loc[j, f'aln_temp{i + 1}']}" for j
                    in range(len(concatenated_pd))]
                monomer_template_hits = []
                for hit in monomer_template_results[i]:
                    # print(f"{hit.name}_{hit.hit_sequence}")
                    if f"{hit.name}_{hit.hit_sequence}" in seen_templates_sequences:
                        continue
                    try:
                        assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[i].seq)
                    except PrefilterError as e:
                        msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                        print(msg)
                        continue
                    monomer_template_hits += [hit]

                curr_pd = create_df(monomer_template_hits)
                curr_pd = curr_pd.add_suffix(f"{i + 1}")
                curr_pd['index'] = curr_pd[f'index{i + 1}']
                curr_pd = curr_pd.drop([f'index{i + 1}'], axis=1)
                if prev_pd is None:
                    prev_pd = curr_pd
                    prev_pd_v2 = curr_pd
                else:
                    prev_pd = prev_pd.merge(curr_pd, how="inner", on='index')
                    prev_pd_v2 = prev_pd_v2.merge(curr_pd, how="outer", on='index')

            concatenated_pd = concatenated_pd.append(prev_pd, ignore_index=True)
            concatenated_pd_v2 = concatenated_pd_v2.append(prev_pd_v2)

        concatenated_pd.to_csv(outdir + '/sequence_templates.csv')
        concatenated_pd_v2.to_csv(outdir + '/sequence_templates_v2.csv')

        cwd = os.getcwd()
        template_dir = outdir + '/templates'
        makedir_if_not_exists(template_dir)
        os.chdir(template_dir)
        for i in range(len(concatenated_pd)):
            for j in range(len(monomer_inputs)):
                template_pdb = concatenated_pd.loc[i, f'template{j + 1}'].split()[0]
                os.system(f"cp {self.atom_dir}/{template_pdb}.atom.gz .")
                os.system(f"gunzip -f {template_pdb}.atom.gz")
        os.chdir(cwd)

        concatenated_pd_nocheck = self.concatenate_templates(monomer_inputs, monomer_template_results, outdir, False)
        concatenated_pd_nocheck = self.filter_same_pdbcodes(concatenated_pd_nocheck, len(monomer_template_results))
        if len(concatenated_pd_nocheck) < 50:

            print(f"template count is smaller than 50, add monomer templates")
            prev_pd = None
            for i in range(len(monomer_template_results)):
                seen_templates_sequences = [
                    f"{concatenated_pd_nocheck.loc[j, f'template{i + 1}']}_{concatenated_pd_nocheck.loc[j, f'aln_temp{i + 1}']}"
                    for j
                    in range(len(concatenated_pd_nocheck))]
                monomer_template_hits = []
                for hit in monomer_template_results[i]:
                    # print(f"{hit.name}_{hit.hit_sequence}")
                    if f"{hit.name}_{hit.hit_sequence}" in seen_templates_sequences:
                        continue
                    try:
                        assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[i].seq)
                    except PrefilterError as e:
                        msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                        print(msg)
                        continue
                    monomer_template_hits += [hit]

                curr_pd = create_df(monomer_template_hits)
                curr_pd = curr_pd.add_suffix(f"{i + 1}")
                curr_pd['index'] = curr_pd[f'index{i + 1}']
                curr_pd = curr_pd.drop([f'index{i + 1}'], axis=1)
                if prev_pd is None:
                    prev_pd = curr_pd
                else:
                    prev_pd = prev_pd.merge(curr_pd, how="inner", on='index')

            concatenated_pd_nocheck = concatenated_pd_nocheck.append(prev_pd, ignore_index=True)

        concatenated_pd_nocheck.to_csv(outdir + '/sequence_templates_nocheck.csv')
