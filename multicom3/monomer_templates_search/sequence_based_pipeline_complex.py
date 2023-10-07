import copy
import os
import sys
import time
from multicom3.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
from multicom3.monomer_templates_concatenation import parsers
from multicom3.monomer_templates_search.sequence_based_pipeline_pdb import assess_hhsearch_hit
from multicom3.tool import hhsearch
from multicom3.tool import hhalign
import dataclasses


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


class monomer_sequence_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

        self.template_searcher = hhsearch.HHSearch(
            binary_path=params['hhsearch_program'],
            databases=[params['complex_sort90_hhsuite_database']],
            input_format='hmm')

        self.pdbdir = params['complex_sort90_atom_dir']

        self.hhmake_program = params['hhmake_program']

        release_date_df = pd.read_csv(params['pdb_release_date_file'])
        self._release_dates = dict(zip(release_date_df['pdbcode'], release_date_df['release_date']))
        self._max_template_date = datetime.datetime.strptime(params['max_template_date'], '%Y-%m-%d')

    def copy_atoms_and_unzip(self, templates, outdir):
        os.chdir(outdir)
        for i in range(len(templates)):
            template_pdb = templates.loc[i, 'name']
            if template_pdb != "empty":
                template_path = os.path.join(self.pdbdir, template_pdb + '.atom.gz')
                os.system(f"cp {template_path} .")
                os.system(f"gunzip -f {template_pdb}.atom.gz")

    def search(self, targetname, sequence, a3m, outdir):
        makedir_if_not_exists(outdir)
        # pdb_hits_out_path = os.path.join(outdir, f'pdb_hits.{self.template_searcher.output_format}')
        pdb_hits_out_path = os.path.join(outdir, f'output.hhr')
        if os.path.exists(pdb_hits_out_path):
            pdb_templates_result = open(pdb_hits_out_path).read()
        else:
            trg_a3m = os.path.join(outdir, targetname + '.a3m')
            trg_hmm = os.path.join(outdir, targetname + '.hmm')
            with open(a3m) as f:
                msa_for_templates = f.read()
                if a3m.find('.sto') > 0:
                    msa_for_templates = parsers.deduplicate_stockholm_msa(msa_for_templates)
                    msa_for_templates = parsers.remove_empty_columns_from_stockholm_msa(msa_for_templates)
                    msa_for_templates = parsers.convert_stockholm_to_a3m(msa_for_templates)

                    with open(trg_a3m, 'w') as fw:
                        fw.write(msa_for_templates)
                else:
                    os.system(f"cp {a3m} {trg_a3m}")

            os.system(f"{self.hhmake_program} -i {trg_a3m} -o {trg_hmm}")
            with open(trg_hmm) as f:
                msa_for_templates = f.read()
                pdb_templates_result = self.template_searcher.query(msa_for_templates, outdir)
                # print(pdb_templates_result)
                # with open(pdb_hits_out_path, 'w') as fw:
                #     fw.write(pdb_templates_result)

        pdb_template_hits = parsers.parse_hhr(hhr_string=pdb_templates_result)

        pdb_template_hits = sorted(pdb_template_hits, key=lambda x: x.sum_probs, reverse=True)

        curr_template_hits = []
        for hit in pdb_template_hits:
            try:
                assess_hhsearch_hit(hit=hit, query_sequence=sequence, max_template_date=self._max_template_date, release_dates=self._release_dates)
            except PrefilterError as e:
                msg = f'hit {hit.name.split()[0]} did not pass prefilter: {str(e)}'
                print(msg)
                continue
            curr_template_hits += [hit]

        curr_pd = create_df(curr_template_hits)

        curr_pd.to_csv(os.path.join(outdir, 'sequence_templates.csv'))

        template_dir = os.path.join(outdir, 'templates')

        makedir_if_not_exists(template_dir)

        self.copy_atoms_and_unzip(templates=curr_pd, outdir=template_dir)
