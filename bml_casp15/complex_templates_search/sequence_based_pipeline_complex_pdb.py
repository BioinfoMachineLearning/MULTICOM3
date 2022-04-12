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

def create_df(targetname, hits):
    row_list = []
    for index, hit in enumerate(hits):
        row_dict = dict(index=index,
                        query=targetname,
                        target=hit.name,
                        alnlen=hit.aligned_cols,
                        sum_probs=hit.sum_probs,
                        qaln=hit.query,
                        qstart=hit.indices_query[0]+1,
                        qend=hit.indices_query[len(hit.indices_query)-1]+1,
                        taln=hit.hit_sequence,
                        tstart=hit.indices_hit[0] + 1,
                        tend=hit.indices_hit[len(hit.indices_hit) - 1] + 1)
        row_list += [row_dict]

    if len(row_list) == 0:
        return pd.DataFrame(columns=['query', 'target', 'alnlen', 'sum_probs', 'qaln', 'qstart', 'qend',
                                     'taln', 'tstart', 'tend'])
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
            databases=[params['complex_sort90_hhsuite_database']],
            input_format='hmm')

        self.seq_template_searcher = hhalign.HHAlign(binary_path=params['hhalign_program'])

        self.cluster_tsv = params['complex_sort90_hhsuite_cluster_tsv']

        self.hhmake_program = params['hhmake_program']

        self.atom_dir = params['complex_sort90_atom_dir']

        self.template_hmm_dir = params['complex_sort90_hmm_dir']

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

    def concatenate_templates(self, monomer_inputs, monomer_template_results, outdir):

        curr_template_hits = []
        for hit in monomer_template_results[0]:
            try:
                assess_hhsearch_hit(hit=hit, query_sequence=monomer_inputs[0].seq)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            curr_template_hits += [hit]

        prev_pd = create_df(monomer_inputs[0].name, curr_template_hits)
        for i in range(1, len(monomer_template_results)):
            template_count = 0
            seen_templates = []
            curr_template_hits = []
            for j in range(len(prev_pd)):
                hit1_name = prev_pd.loc[j, 'target'].split()[0]
                print(f"finding hits for {hit1_name}")
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
                            if not os.path.exists(self.template_hmm_dir + '/' + hit_name + '.hmm'):
                                print(f"cannot find {hit_name}.hmm in {self.template_hmm_dir}")
                                continue
                            hit = self.align_template(
                                                      monomer_inputs[i].hmm_path,
                                                      self.template_hmm_dir + '/' + hit_name + '.hmm',
                                                      outdir)

                    if hit is None:
                        continue

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

            curr_pd = create_df(monomer_inputs[i].name, curr_template_hits)
            # print(curr_pd)
            prev_pd = prev_pd.merge(curr_pd, how="inner", on='index', suffixes=(str(i), str(i + 1)))

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

        concatenated_pd.to_csv(outdir + '/sequence_templates.csv')

        cwd = os.getcwd()
        template_dir = outdir + '/templates'
        makedir_if_not_exists(template_dir)
        os.chdir(template_dir)
        for i in range(len(concatenated_pd)):
            for j in range(len(monomer_inputs)):
                template_pdb = concatenated_pd.loc[i, f'target{j+1}'].split()[0]
                os.system(f"cp {self.atom_dir}/{template_pdb}.atom.gz .")
                os.system(f"gunzip -f {template_pdb}.atom.gz")
        os.chdir(cwd)