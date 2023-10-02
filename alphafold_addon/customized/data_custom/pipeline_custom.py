# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Functions for building the input features for the AlphaFold model."""

import os
from typing import Any, Mapping, MutableMapping, Optional, Sequence, Union
from absl import logging
from alphafold.common import residue_constants
from alphafold.data_custom import custom_params
from alphafold.data_custom import msa_identifiers
from alphafold.data_custom import parsers
from alphafold.data_custom import templates
from alphafold.data_custom.tools import hhblits
from alphafold.data_custom.tools import hhsearch
from alphafold.data_custom.tools import hmmsearch
from alphafold.data_custom.tools import jackhmmer
import numpy as np

# Internal import (7716).

FeatureDict = MutableMapping[str, np.ndarray]
TemplateSearcher = Union[hhsearch.HHSearch, hmmsearch.Hmmsearch]

def read_fasta(fileobj):
    current_sequence = ""
    current_id = None
    for line in fileobj:
        if line.startswith(">"):
            if current_id is not None:
                yield current_id, current_sequence
            current_id = line.rstrip()[1:]
            current_sequence = ""
        elif not line.startswith(";"):
            current_sequence += line.rstrip()
    yield current_id, current_sequence


def combine_a3ms(a3ms, outa3m):
    with open(outa3m, 'w') as fw:
        query_name, query_seq = None, None
        for a3m in a3ms:
            with open(a3m, 'r') as fileobj:
                for i, (seq_id, seq) in enumerate(read_fasta(fileobj)):
                    if i == 0:
                        if query_name is None and query_seq is None:
                            query_name = seq_id
                            query_seq = seq
                            fw.write(f">{seq_id}\n{seq}\n")
                        elif query_name != seq_id or query_seq != seq:
                            raise ValueError("The input a3ms don't have the same query name or query sequences")
                    else:
                        fw.write(f">{seq_id}\n{seq}\n")


def make_sequence_features(
        sequence: str, description: str, num_res: int) -> FeatureDict:
    """Constructs a feature dict of sequence features."""
    features = {}
    features['aatype'] = residue_constants.sequence_to_onehot(
        sequence=sequence,
        mapping=residue_constants.restype_order_with_x,
        map_unknown_to_x=True)
    features['between_segment_residues'] = np.zeros((num_res,), dtype=np.int32)
    features['domain_name'] = np.array([description.encode('utf-8')],
                                       dtype=np.object_)
    features['residue_index'] = np.array(range(num_res), dtype=np.int32)
    features['seq_length'] = np.array([num_res] * num_res, dtype=np.int32)
    features['sequence'] = np.array([sequence.encode('utf-8')], dtype=np.object_)
    return features


def make_msa_features(msas: Sequence[parsers.Msa], msa_output_dir: str, msa_save_path: str, filter=True) -> FeatureDict:
    """Constructs a feature dict of MSA features."""
    if not msas:
        raise ValueError('At least one MSA must be provided.')

    int_msa = []
    deletion_matrix = []
    species_ids = []
    seen_desc = []
    seen_sequences = []
    for msa_index, msa in enumerate(msas):
        if not msa:
            raise ValueError(f'MSA {msa_index} must contain at least one sequence.')
        for sequence_index, sequence in enumerate(msa.sequences):
            if filter and sequence in seen_sequences:
                continue
            seen_sequences += [sequence]
            seen_desc += [msa.descriptions[sequence_index]]
            int_msa.append([residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence])
            deletion_matrix.append(msa.deletion_matrix[sequence_index])
            identifiers = msa_identifiers.get_identifiers(msa.descriptions[sequence_index])
            species_ids.append(identifiers.species_id.encode('utf-8'))

    num_res = len(msas[0].sequences[0])
    num_alignments = len(int_msa)
    features = {}
    features['deletion_matrix_int'] = np.array(deletion_matrix, dtype=np.int32)
    features['msa'] = np.array(int_msa, dtype=np.int32)
    features['num_alignments'] = np.array([num_alignments] * num_res, dtype=np.int32)
    features['msa_species_identifiers'] = np.array(species_ids, dtype=np.object_)

    with open(os.path.join(msa_output_dir, msa_save_path), 'w') as fw:
        for (desc, seq) in zip(seen_desc, seen_sequences):
            fw.write(f'>{desc}\n{seq}\n')

    return features


def run_msa_tool(msa_runner, input_fasta_path: str, msa_out_path: str,
                 msa_format: str, use_precomputed_msas: bool,
                 max_sto_sequences: Optional[int] = None
                 ) -> Mapping[str, Any]:
    """Runs an MSA tool, checking if output already exists first."""
    if not use_precomputed_msas or not os.path.exists(msa_out_path):
        if msa_format == 'sto' and max_sto_sequences is not None:
            result = msa_runner.query(input_fasta_path, max_sto_sequences)[0]  # pytype: disable=wrong-arg-count
        else:
            result = msa_runner.query(input_fasta_path)[0]
        with open(msa_out_path, 'w') as f:
            f.write(result[msa_format])
    else:
        logging.warning('Reading MSA from file %s', msa_out_path)
        if msa_format == 'sto' and max_sto_sequences is not None:
            precomputed_msa = parsers.truncate_stockholm_msa(
                msa_out_path, max_sto_sequences)
            result = {'sto': precomputed_msa}
        else:
            with open(msa_out_path, 'r') as f:
                result = {msa_format: f.read()}
    return result


class DataPipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self,
                 jackhmmer_binary_path: str,
                 hhblits_binary_path: str,
                 uniref90_database_path: str,
                 mgnify_database_path: str,
                 bfd_database_path: Optional[str],
                 uniref30_database_path: Optional[str],
                 small_bfd_database_path: Optional[str],
                 template_searcher: TemplateSearcher,
                 template_featurizer: templates.TemplateHitFeaturizer,
                 use_small_bfd: bool,
                 mgnify_max_hits: int = 501,
                 uniref_max_hits: int = 10000,
                 use_precomputed_msas: bool = False):
        """Initializes the data pipeline."""
        self._use_small_bfd = use_small_bfd
        self.jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
            binary_path=jackhmmer_binary_path,
            database_path=uniref90_database_path)
        if use_small_bfd:
            self.jackhmmer_small_bfd_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=small_bfd_database_path)
        else:
            self.hhblits_bfd_uniref_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[bfd_database_path, uniref30_database_path])
        self.jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
            binary_path=jackhmmer_binary_path,
            database_path=mgnify_database_path)
        self.template_searcher = template_searcher
        self.template_featurizer = template_featurizer
        self.mgnify_max_hits = mgnify_max_hits
        self.uniref_max_hits = uniref_max_hits
        self.use_precomputed_msas = use_precomputed_msas

    def process(self, input_fasta_path: str, msa_output_dir: str, template_output_dir: str, custom_inputs) -> FeatureDict:
        """Runs alignment tools on the input sequence and creates features."""
        with open(input_fasta_path) as f:
            input_fasta_str = f.read()
        input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)
        if len(input_seqs) != 1:
            raise ValueError(
                f'More than one input sequence found in {input_fasta_path}.')
        input_sequence = input_seqs[0]
        input_description = input_descs[0]
        num_res = len(input_sequence)

        sequence_features = make_sequence_features(
            sequence=input_sequence,
            description=input_description,
            num_res=num_res)

        custom_result = None
        if custom_inputs.custom_msa is not None:
            if custom_inputs.custom_msa.find('.a3m') > 0:
                custom_msa_out_path = os.path.join(msa_output_dir, 'custom.a3m')
                os.system(f"cp {custom_inputs.custom_msa} {custom_msa_out_path}")
                with open(custom_msa_out_path, 'r') as f:
                    custom_result = parsers.parse_a3m(f.read())
            else:
                custom_msa_out_path = os.path.join(msa_output_dir, 'custom.sto')
                os.system(f"cp {custom_inputs.custom_msa} {custom_msa_out_path}")
                with open(custom_msa_out_path, 'r') as f:
                    custom_result = parsers.parse_stockholm(f.read())

        msa_features = None
        jackhmmer_uniref90_result = None
        if custom_result is None:
            uniref90_out_path = os.path.join(msa_output_dir, 'uniref90_hits.sto')
            if custom_inputs.uniref90_sto is not None:
                os.system(f"cp {custom_inputs.uniref90_sto} {uniref90_out_path}")
            jackhmmer_uniref90_result = run_msa_tool(
                msa_runner=self.jackhmmer_uniref90_runner,
                input_fasta_path=input_fasta_path,
                msa_out_path=uniref90_out_path,
                msa_format='sto',
                use_precomputed_msas=self.use_precomputed_msas,
                max_sto_sequences=self.uniref_max_hits)

            uniref90_msa = parsers.parse_stockholm(jackhmmer_uniref90_result['sto'])
            uniref90_msa = uniref90_msa.truncate(max_seqs=self.uniref_max_hits)

            mgnify_out_path = os.path.join(msa_output_dir, 'mgnify_hits.sto')
            if custom_inputs.mgnify_sto is not None:
                os.system(f"cp {custom_inputs.mgnify_sto} {mgnify_out_path}")
            jackhmmer_mgnify_result = run_msa_tool(
                msa_runner=self.jackhmmer_mgnify_runner,
                input_fasta_path=input_fasta_path,
                msa_out_path=mgnify_out_path,
                msa_format='sto',
                use_precomputed_msas=self.use_precomputed_msas,
                max_sto_sequences=self.mgnify_max_hits)

            mgnify_msa = parsers.parse_stockholm(jackhmmer_mgnify_result['sto'])
            mgnify_msa = mgnify_msa.truncate(max_seqs=self.mgnify_max_hits)

            bfd_out_path = os.path.join(msa_output_dir, 'bfd_uniref_hits.a3m')
            if custom_inputs.bfd_uniref_a3m is not None:
                os.system(f"cp {custom_inputs.bfd_uniref_a3m} {bfd_out_path}")
            elif custom_inputs.bfd_a3m is not None and custom_inputs.uniref_a3m is not None:
                combine_a3ms([custom_inputs.bfd_a3m, custom_inputs.uniref_a3m], bfd_out_path)

            hhblits_bfd_uniref_result = run_msa_tool(
                msa_runner=self.hhblits_bfd_uniref_runner,
                input_fasta_path=input_fasta_path,
                msa_out_path=bfd_out_path,
                msa_format='a3m',
                use_precomputed_msas=self.use_precomputed_msas)
            bfd_msa = parsers.parse_a3m(hhblits_bfd_uniref_result['a3m'])
            msa_features = make_msa_features((uniref90_msa, bfd_msa, mgnify_msa), msa_output_dir,
                                             msa_save_path='monomer_final.a3m')
            logging.info('Uniref90 MSA size: %d sequences.', len(uniref90_msa))
            logging.info('BFD MSA size: %d sequences.', len(bfd_msa))
            logging.info('MGnify MSA size: %d sequences.', len(mgnify_msa))
        else:
            logging.info('Custom MSA size: %d sequences.', len(custom_result))
            msa_features = make_msa_features([custom_result], msa_output_dir, msa_save_path='monomer_final.a3m')

        logging.info('Final (deduplicated) MSA size: %d sequences.',
                     msa_features['num_alignments'][0])

        templates_result_features = None
        if custom_inputs.notemplate:
            templates_result_features = mk_mock_template(input_sequence)
        elif custom_inputs.temp_struct_csv is not None:
            templates_result_features = self.template_featurizer.get_templates(query_sequence=input_sequence,
                                                                               template_pdb_dir=template_output_dir,
                                                                               hits_file=custom_inputs.temp_struct_csv).features
        else:
            if jackhmmer_uniref90_result is None:
                uniref90_out_path = os.path.join(msa_output_dir, 'uniref90_hits.sto')
                if custom_inputs.uniref90_sto is not None:
                    os.system(f"cp {custom_inputs.uniref90_sto} {uniref90_out_path}")
                jackhmmer_uniref90_result = run_msa_tool(
                    msa_runner=self.jackhmmer_uniref90_runner,
                    input_fasta_path=input_fasta_path,
                    msa_out_path=uniref90_out_path,
                    msa_format='sto',
                    use_precomputed_msas=self.use_precomputed_msas,
                    max_sto_sequences=self.uniref_max_hits)

            msa_for_templates = jackhmmer_uniref90_result['sto']
            msa_for_templates = parsers.deduplicate_stockholm_msa(msa_for_templates)
            msa_for_templates = parsers.remove_empty_columns_from_stockholm_msa(
                msa_for_templates)

            if self.template_searcher.input_format == 'sto':
                pdb_templates_result = self.template_searcher.query(msa_for_templates)
            elif self.template_searcher.input_format == 'a3m':
                uniref90_msa_as_a3m = parsers.convert_stockholm_to_a3m(msa_for_templates)
                pdb_templates_result = self.template_searcher.query(uniref90_msa_as_a3m)
            else:
                raise ValueError('Unrecognized template input format: '
                                 f'{self.template_searcher.input_format}')

            pdb_hits_out_path = os.path.join(
                msa_output_dir, f'pdb_hits.{self.template_searcher.output_format}')
            with open(pdb_hits_out_path, 'w') as f:
                f.write(pdb_templates_result)

            pdb_template_hits = self.template_searcher.get_template_hits(
                output_string=pdb_templates_result, input_sequence=input_sequence)

            templates_result = self.template_featurizer.get_templates(
                query_sequence=input_sequence,
                hits=pdb_template_hits)
            templates_result_features = templates_result.features

        if custom_inputs.notemplate:
            logging.info('Total number of templates (NB: this can include bad '
                         'templates and is later filtered to top 4): %d.',
                         len(templates_result_features['template_domain_names']))
        else:
            logging.info('Total number of templates (NB: this can include bad '
                         'templates and is later filtered to top 4): %d.',
                         templates_result_features['template_domain_names'].shape[0])

        return {**sequence_features, **msa_features, **templates_result_features}
