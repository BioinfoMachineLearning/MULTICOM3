import copy
import os
import sys
import time
from multicom3.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
from multicom3.monomer_templates_concatenation import parsers
from multicom3.monomer_templates_search.sequence_based_pipeline_pdb import assess_hhsearch_hit, PrefilterError
from multicom3.tool import hhsearch
from multicom3.tool import hhalign
import dataclasses
import numpy as np

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

        release_date_df = pd.read_csv(params['pdb_release_date_file'])
        self._release_dates = dict(zip(release_date_df['pdbcode'], pdb_release_date_df['release_date']))
        self._max_template_date = datetime.datetime.strptime(params['max_template_date'], '%Y-%m-%d')

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

        src_hmm = os.path.join(outdir, src_name + '.hmm')
        os.system(f"{self.hhmake_program} -i {src_a3m} -o {src_hmm}")

        pdb_templates_result = self.single_seq_template_searcher.query(src_hmm, trg_fasta, outdir)

        pdb_templates_hits = parsers.parse_hhr(hhr_string=pdb_templates_result)

        if len(pdb_templates_hits) > 1:
            print("HHalign returns more than 1 result!")

        return pdb_templates_hits[0]

    def concatenate_templates(self, monomer_inputs, monomer_template_results, outdir, check_hit=True):

        curr_template_hits = []
        for hit in monomer_template_results[0]:
            if check_hit:
                try:
                    assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[0].seq, max_template_date=self._max_template_date, release_dates=self._release_dates)
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
            for j in range(len(prev_pd)):
                hit1_name = prev_pd.loc[j, f'name{i}'].split()[0]
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
                            fasta_file = os.path.join(self.pdb_seq_dir, hit_name + '.fasta')
                            if not os.path.exists(fasta_file):
                                continue
                            hit = self.align_template(monomer_inputs[i].name,
                                                      monomer_inputs[i].msa_path,
                                                      fasta_file,
                                                      outdir)

                    if hit is None:
                        continue

                    if check_hit:
                        try:
                            assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[i].seq, max_template_date=self._max_template_date, release_dates=self._release_dates)
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
        prev_pd = prev_pd.sort_values(by=['max_probs'], ascending=False)
        return prev_pd

    def filter_same_pdbcodes(self, indf, monomer_count):

        indf_sorted = indf.sort_values(by=['tpdbcode1','max_probs'], ascending=False)
        indf_sorted.reset_index(inplace=True, drop=True)

        keep_indices = []
        pdbcodes = []
        cover_chains_in_pdb = {}
        for i in range(len(indf_sorted)):
            chain_count = len(set([indf_sorted.loc[i, f'name{j + 1}'] for j in range(monomer_count)]))
            if indf_sorted.loc[i, 'tpdbcode1'] not in pdbcodes:
                if len(pdbcodes) > 0:
                    max_index = -1
                    max_count = 0
                    for index in cover_chains_in_pdb:
                        if cover_chains_in_pdb[index] > max_count:
                            max_index = index
                            max_count = cover_chains_in_pdb[index]
                    keep_indices += [max_index]

                pdbcodes += [indf_sorted.loc[i, 'tpdbcode1']]
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

            monomer_outdir = os.path.join(outdir, monomer_input.name)

            makedir_if_not_exists(monomer_outdir)

            pdb_hits_out_path = os.path.join(monomer_outdir, 'output.hhr')
            if os.path.exists(pdb_hits_out_path):
                pdb_templates_result = open(pdb_hits_out_path).read()
                monomer_input.msa_path = os.path.join(monomer_outdir, monomer_input.name + '.a3m')
            else:
                with open(monomer_input.msa_path) as f:
                    msa_for_templates = f.read()
                    if monomer_input.msa_path.find('.sto') > 0:
                        msa_for_templates = parsers.deduplicate_stockholm_msa(msa_for_templates)
                        msa_for_templates = parsers.remove_empty_columns_from_stockholm_msa(msa_for_templates)
                        msa_for_templates = parsers.convert_stockholm_to_a3m(msa_for_templates)

                        with open(os.path.join(monomer_outdir, monomer_input.name + '.a3m'), 'w') as fw:
                            fw.write(msa_for_templates)

                        monomer_input.msa_path = os.path.join(monomer_outdir, monomer_input.name + '.a3m')

                    pdb_templates_result = self.template_searcher.query(msa_for_templates, monomer_outdir)

                    with open(pdb_hits_out_path, 'w') as fw:
                        fw.write(pdb_templates_result)

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
                seen_templates_sequences = [f"{concatenated_pd.loc[j, f'name{i+1}']}_{concatenated_pd.loc[j, f'hit_sequence{i+1}']}" for j in range(len(concatenated_pd))]
                monomer_template_hits = []
                for hit in monomer_template_results[i]:
                    if f"{hit.name.split()[0]}_{hit.hit_sequence}" in seen_templates_sequences:
                        continue
                    try:
                        assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[i].seq, max_template_date=self._max_template_date, release_dates=self._release_dates)
                    except PrefilterError as e:
                        msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                        print(msg)
                        continue
                    hit.name = hit.name.split()[0]
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

            concatenated_pd = concatenated_pd.append(prev_pd)
            concatenated_pd_v2 = concatenated_pd_v2.append(prev_pd_v2)

        concatenated_pd.reset_index(inplace=True, drop=True)
        concatenated_pd.to_csv(os.path.join(outdir, 'sequence_templates.csv'))
        concatenated_pd_v2.reset_index(inplace=True, drop=True)
        concatenated_pd_v2.to_csv(os.path.join(outdir, 'sequence_templates_v2.csv'))
