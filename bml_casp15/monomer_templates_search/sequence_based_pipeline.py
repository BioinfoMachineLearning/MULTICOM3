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
        empty_dict = dict(index=0,
                          name='empty',
                          tpdbcode='empty',
                          aligned_cols=0,
                          sum_probs=0,
                          query='empty',
                          hit_sequence='empty',
                          indices_query=None,
                          indices_hit=None)
        row_list += [empty_dict]
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


class monomer_sequence_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

        self.template_searcher = hhsearch.HHSearch(
            binary_path=params['hhsearch_program'],
            databases=[params['pdb_sort90_hhsuite_database']])

        self.hhmake_program = params['hhmake_program']

        self.pdbdir = params['pdb_sort90_atom_dir']

    def copy_atoms_and_unzip(self, templates, outdir):
        os.chdir(outdir)
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'target']
            os.system(f"cp {self.pdbdir}/{template_pdb} {outdir}")
            os.system(f"gunzip -f {template_pdb}")

    def search(self, fasta_file, a3m, outdir):

        targetname = None
        seq = None

        for line in open(fasta_file):
            line = line.rstrip('\n').strip()
            if line.startswith('>'):
                targetname = line[1:]
            else:
                seq = line

        os.system(f"cp {a3m} {outdir}/{targetname}.a3m")

        os.system(f"{self.hhmake_program} -i {outdir}/{targetname}.a3m -o {outdir}/{targetname}.hmm\n")

        pdb_hits_out_path = os.path.join(outdir, f'pdb_hits.{self.template_searcher.output_format}')
        if os.path.exists(pdb_hits_out_path):
            pdb_templates_result = open(pdb_hits_out_path).read()
        else:
            with open(f"{outdir}/{targetname}.a3m") as f:
                msa_for_templates = f.read()
                pdb_templates_result = self.template_searcher.query(msa_for_templates, outdir)
                with open(pdb_hits_out_path, 'w') as fw:
                    fw.write(pdb_templates_result)

        pdb_template_hits = parsers.parse_hhr(hhr_string=pdb_templates_result)

        pdb_template_hits = sorted(pdb_template_hits, key=lambda x: x.sum_probs, reverse=True)

        curr_template_hits = []
        for hit in pdb_template_hits:
            try:
                assess_hhsearch_hit(hit=hit, query_sequence=seq)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            curr_template_hits += [hit]

        curr_pd = create_df(curr_template_hits)

        curr_pd.to_csv(outdir + '/sequence_templates.csv')

        template_dir = outdir + '/templates'

        makedir_if_not_exists(template_dir)

        self.copy_atoms_and_unzip(templates=curr_pd, outdir=template_dir)
