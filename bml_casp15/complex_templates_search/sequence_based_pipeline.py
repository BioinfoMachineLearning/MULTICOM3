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


def create_df(hits):
    row_list = []
    for index, hit in enumerate(hits):
        row_dict = dict(index=index,
                        name=hit.name,
                        tpdbcode=hit.name[0:4],
                        aligned_cols=hit.aligned_cols,
                        sum_probs=hit.sum_probs,
                        query=hit.query,
                        hit_sequence=hit.hit_sequence,
                        indices_query='_'.join([str(i) for i in hit.indices_query]),
                        indices_hit='_'.join([str(i) for i in hit.indices_hit]))
        row_list += [row_dict]

    if len(row_list) == 0:
        return pd.DataFrame(columns=['index', 'name', 'tpdbcode', 'aligned_cols', 'sum_probs',
                                     'query', 'hit_sequence', 'indices_query', 'indices_hit'])
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
            binary_path=params['hhsearch_program'],
            databases=[params['pdb70_hhsuite_database']])

        self.single_seq_template_searcher = hhalign.HHAlign(binary_path=params['hhalign_program'])

        self.cluster_tsv = params['pdb70_hhsuite_cluster_tsv']

        self.hhmake_program = params['hhmake_program']

        self.pdb_seqs_dir = params['pdb_seqs_dir']

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

    def align_template(self, src_name, src_a3m, trg_fasta, outdir):

        os.system(f"{self.hhmake_program} -i {src_a3m} -o {outdir}/{src_name}.hmm\n")

        pdb_templates_result = self.single_seq_template_searcher.query(f"{outdir}/{src_name}.hmm", trg_fasta, outdir)

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
        for i in range(1, len(monomer_template_results)):
            template_count = 0
            seen_templates = []
            curr_template_hits = []
            for j in range(len(prev_pd)):
                if i == 1:
                    hit1_name = prev_pd.loc[j, 'name'].split()[0]
                else:
                    hit1_name = prev_pd.loc[j, f'name{i}'].split()[0]
                if template_count > 50:
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
                            if not os.path.exists(self.pdb_seqs_dir + '/' + hit_name + '.fasta'):
                                continue
                            hit = self.align_template(monomer_inputs[i].name,
                                                      monomer_inputs[i].msa_path,
                                                      self.pdb_seqs_dir + '/' + hit_name + '.fasta',
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
                        break

            curr_pd = create_df(curr_template_hits)
            # print(curr_pd)
            prev_pd = prev_pd.merge(curr_pd, how="inner", on='index', suffixes=(str(i), str(i + 1)))

        return prev_pd

    def filter_same_pdbcodes(self, indf):
        pdbcodes = []
        keep_indices = []
        for i in range(len(indf)):
            pdbcode = indf.loc[i, 'tpdbcode1']
            if pdbcode in pdbcodes:
                continue
            keep_indices += [i]
            pdbcodes += [pdbcode]
        return indf.iloc[keep_indices]

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

                    pdb_templates_result = self.template_searcher.query(msa_for_templates, monomer_outdir)

                    with open(pdb_hits_out_path, 'w') as fw:
                        fw.write(pdb_templates_result)

            pdb_template_hits = parsers.parse_hhr(hhr_string=pdb_templates_result)

            pdb_template_hits = sorted(pdb_template_hits, key=lambda x: x.sum_probs, reverse=True)

            monomer_template_results += [pdb_template_hits]

            monomer_names += [monomer_input.name]

        concatenated_pd = self.concatenate_templates(monomer_inputs, monomer_template_results, outdir)

        concatenated_pd = self.filter_same_pdbcodes(concatenated_pd)

        if len(concatenated_pd) < 50:
            print(f"template count is smaller than 50, add monomer templates")
            prev_pd = None
            for i in range(len(monomer_template_results)):
                seen_templates_sequences = [f"{concatenated_pd.loc[j, f'name{j+1}']}_{concatenated_pd.loc[j, f'hit_sequence{j+1}']}" for j in range(len(concatenated_pd))]
                monomer_template_hits = []
                for hit in monomer_template_results[i]:
                    if f"{hit.name.split()[0]}_{hit.hit_sequence}" in seen_templates_sequences:
                        continue
                    try:
                        assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[i].seq)
                    except PrefilterError as e:
                        msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                        print(msg)
                        continue
                    hit.name = hit.name.split()[0]
                    monomer_template_hits += [hit]

                curr_pd = create_df(monomer_template_hits)
                if prev_pd is None:
                    prev_pd = curr_pd
                else:
                    prev_pd = prev_pd.merge(curr_pd, how="inner", on='index', suffixes=(str(i), str(i + 1)))

            concatenated_pd = concatenated_pd.append(prev_pd)

        concatenated_pd.to_csv(outdir + '/sequence_templates.csv')

        concatenated_pd_nocheck = self.concatenate_templates(monomer_inputs, monomer_template_results, outdir, False)

        concatenated_pd_nocheck = self.filter_same_pdbcodes(concatenated_pd_nocheck)

        if len(concatenated_pd_nocheck) < 50:
            print(f"template count is smaller than 50, add monomer templates")
            prev_pd = None
            for i in range(len(monomer_template_results)):
                seen_templates_sequences = [
                    f"{concatenated_pd_nocheck.loc[j, f'name{j + 1}']}_{concatenated_pd_nocheck.loc[j, f'hit_sequence{j + 1}']}" for j
                    in range(len(concatenated_pd_nocheck))]
                monomer_template_hits = []
                for hit in monomer_template_results[i]:
                    if f"{hit.name.split()[0]}_{hit.hit_sequence}" in seen_templates_sequences:
                        continue
                    try:
                        assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[i].seq)
                    except PrefilterError as e:
                        msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                        print(msg)
                        continue
                    hit.name = hit.name.split()[0]
                    monomer_template_hits += [hit]

                curr_pd = create_df(monomer_template_hits)
                if prev_pd is None:
                    prev_pd = curr_pd
                else:
                    prev_pd = prev_pd.merge(curr_pd, how="inner", on='index', suffixes=(str(i), str(i + 1)))

            concatenated_pd_nocheck = concatenated_pd_nocheck.append(prev_pd)

        concatenated_pd_nocheck.to_csv(outdir + '/sequence_templates_nocheck.csv')
