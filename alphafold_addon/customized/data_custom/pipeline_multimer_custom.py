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

"""Functions for building the features for the AlphaFold multimer model."""

import collections
import contextlib
import copy
import dataclasses
import json
import os
import tempfile
from typing import Mapping, MutableMapping, Sequence

from absl import logging
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data_custom import custom_params
from alphafold.data_custom import feature_processing
from alphafold.data_custom import msa_pairing
from alphafold.data_custom import parsers
from alphafold.data_custom import pipeline_custom
from alphafold.data_custom import templates
from alphafold.data_custom import templates_custom
from alphafold.data_custom.tools import jackhmmer
import numpy as np
import pickle
from alphafold.data_custom import msa_identifiers
import pandas as pd


# Internal import (7716).

def _make_chain_id_map(*,
                       sequences: Sequence[str],
                       descriptions: Sequence[str],
                       ) -> Mapping[str, templates_custom.FastaChain]:
    """Makes a mapping from PDB-format chain ID to sequence and description."""
    if len(sequences) != len(descriptions):
        raise ValueError('sequences and descriptions must have equal length. '
                         f'Got {len(sequences)} != {len(descriptions)}.')
    if len(sequences) > protein.PDB_MAX_CHAINS:
        raise ValueError('Cannot process more chains than the PDB format supports. '
                         f'Got {len(sequences)} chains.')
    chain_id_map = {}
    chain_id_seq_map = {}
    for chain_id, sequence, description in zip(
            protein.PDB_CHAIN_IDS, sequences, descriptions):
        chain_id_map[chain_id] = templates_custom.FastaChain(
            sequence=sequence, description=description)
        chain_id_seq_map[chain_id] = sequence
    return chain_id_map, chain_id_seq_map


@contextlib.contextmanager
def temp_fasta_file(fasta_str: str):
    with tempfile.NamedTemporaryFile('w', suffix='.fasta') as fasta_file:
        fasta_file.write(fasta_str)
        fasta_file.seek(0)
        yield fasta_file.name


def convert_monomer_features(
        monomer_features: pipeline_custom.FeatureDict,
        chain_id: str) -> pipeline_custom.FeatureDict:
    """Reshapes and modifies monomer features for multimer models."""
    converted = {}
    converted['auth_chain_id'] = np.asarray(chain_id, dtype=np.object_)
    unnecessary_leading_dim_feats = {
        'sequence', 'domain_name', 'num_alignments', 'seq_length'}
    for feature_name, feature in monomer_features.items():
        if feature_name in unnecessary_leading_dim_feats:
            # asarray ensures it's a np.ndarray.
            feature = np.asarray(feature[0], dtype=feature.dtype)
        elif feature_name == 'aatype':
            # The multimer model performs the one-hot operation itself.
            feature = np.argmax(feature, axis=-1).astype(np.int32)
        elif feature_name == 'template_aatype':
            feature = np.argmax(feature, axis=-1).astype(np.int32)
            new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
            feature = np.take(new_order_list, feature.astype(np.int32), axis=0)
        elif feature_name == 'template_all_atom_masks':
            feature_name = 'template_all_atom_mask'
        converted[feature_name] = feature
    return converted


def int_id_to_str_id(num: int) -> str:
    """Encodes a number as a string, using reverse spreadsheet style naming.

  Args:
    num: A positive integer.

  Returns:
    A string that encodes the positive integer using reverse spreadsheet style,
    naming e.g. 1 = A, 2 = B, ..., 27 = AA, 28 = BA, 29 = CA, ... This is the
    usual way to encode chain IDs in mmCIF files.
  """
    if num <= 0:
        raise ValueError(f'Only positive integers allowed, got {num}.')

    num = num - 1  # 1-based indexing.
    output = []
    while num >= 0:
        output.append(chr(num % 26 + ord('A')))
        num = num // 26 - 1
    return ''.join(output)


def add_assembly_features(
        all_chain_features: MutableMapping[str, pipeline_custom.FeatureDict],
) -> MutableMapping[str, pipeline_custom.FeatureDict]:
    """Add features to distinguish between chains.

  Args:
    all_chain_features: A dictionary which maps chain_id to a dictionary of
      features for each chain.

  Returns:
    all_chain_features: A dictionary which maps strings of the form
      `<seq_id>_<sym_id>` to the corresponding chain features. E.g. two
      chains from a homodimer would have keys A_1 and A_2. Two chains from a
      heterodimer would have keys A_1 and B_1.
  """
    # Group the chains by sequence
    seq_to_entity_id = {}
    grouped_chains = collections.defaultdict(list)
    for chain_id, chain_features in all_chain_features.items():
        seq = str(chain_features['sequence'])
        if seq not in seq_to_entity_id:
            seq_to_entity_id[seq] = len(seq_to_entity_id) + 1
        grouped_chains[seq_to_entity_id[seq]].append(chain_features)

    new_all_chain_features = {}
    chain_id = 1
    for entity_id, group_chain_features in grouped_chains.items():
        for sym_id, chain_features in enumerate(group_chain_features, start=1):
            new_all_chain_features[
                f'{int_id_to_str_id(entity_id)}_{sym_id}'] = chain_features
            seq_length = chain_features['seq_length']
            chain_features['asym_id'] = chain_id * np.ones(seq_length)
            chain_features['sym_id'] = sym_id * np.ones(seq_length)
            chain_features['entity_id'] = entity_id * np.ones(seq_length)
            chain_id += 1

    return new_all_chain_features


def pad_msa(np_example, min_num_seq):
    np_example = dict(np_example)
    num_seq = np_example['msa'].shape[0]
    if num_seq < min_num_seq:
        for feat in ('msa', 'deletion_matrix', 'bert_mask', 'msa_mask'):
            np_example[feat] = np.pad(
                np_example[feat], ((0, min_num_seq - num_seq), (0, 0)))
        np_example['cluster_bias_mask'] = np.pad(
            np_example['cluster_bias_mask'], ((0, min_num_seq - num_seq),))
    return np_example

FeatureDict = MutableMapping[str, np.ndarray]

def make_msa_features(msas: Sequence[parsers.Msa], msa_output_dir: str, msa_save_path: str) -> FeatureDict:
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
            # if sequence in seen_sequences:
            #     continue
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


def save_paired_msas(chain_id_map, all_chain_features, msa_output_dir):
    # all_chain_features keys and pair rows order: A_1, A_2, B_1, B_2 ....
    # chain_id_map: keys: A, B, C, D
    paired_rows = None

    paired_npy = os.path.join(msa_output_dir, 'pair_msas.npy')

    if not os.path.exists(paired_npy):
        return paired_rows

    with open(paired_npy, 'rb') as f:
        paired_rows = np.load(f)

    # Group the chains by sequence
    grouped_chains = collections.defaultdict(list)
    seq_to_entity_id = {}
    for chain_id, chain_features in all_chain_features.items():
        seq = str(chain_features['sequence'])
        if seq not in seq_to_entity_id:
            seq_to_entity_id[seq] = len(seq_to_entity_id)
        grouped_chains[seq_to_entity_id[seq]].append(chain_id)

    chain_reorder = {}
    chain_num = 0
    for entity_id, chain_ids in grouped_chains.items():
        for chain_id in chain_ids:
            chain_reorder[chain_num] = chain_id
            chain_num += 1

    for chain_num in chain_reorder:
        final_a3m = os.path.join(msa_output_dir, chain_reorder[chain_num], 'multimer_final.a3m')
        with open(final_a3m) as f:
            final_msa = parsers.parse_a3m(f.read())

        seen_desc = []
        seen_sequences = []
        for sequence_index in list(paired_rows[:, chain_num]):
            if sequence_index == -1:
                seen_desc += ['placeholder']
                seen_sequences += ['-' * len(chain_id_map[chain_reorder[chain_num]].sequence)]
            else:
                seen_desc += [final_msa.descriptions[sequence_index]]
                seen_sequences += [final_msa.sequences[sequence_index]]

        with open(os.path.join(msa_output_dir, chain_reorder[chain_num] + ".paired.a3m"), 'w') as fw:
            for (desc, seq) in zip(seen_desc, seen_sequences):
                fw.write(f'>{desc}\n{seq}\n')


class DataPipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self,
                 template_searcher=None,
                 template_featurizer=None,
                 monomer_template_featurizer=None,
                 max_uniprot_hits: int = 50000):
        """Initializes the data pipeline.

    Args:
      monomer_data_pipeline: An instance of pipeline.DataPipeline - that runs
        the data pipeline for the monomer AlphaFold system.
      jackhmmer_binary_path: Location of the jackhmmer binary.
      uniprot_database_path: Location of the unclustered uniprot sequences, that
        will be searched with jackhmmer and used for MSA pairing.
      max_uniprot_hits: The maximum number of hits to return from uniprot.
      use_precomputed_msas: Whether to use pre-existing MSAs; see run_alphafold.
    """
        self._max_uniprot_hits = max_uniprot_hits
        self.template_searcher = template_searcher
        self.template_featurizer = template_featurizer
        self.monomer_template_featurizer = monomer_template_featurizer

    def _process_single_chain(
            self,
            chain_id: str,
            sequence: str,
            description: str,
            chain_monomer_msa: str,
            chain_multimer_msa: str,
            msa_output_dir: str,
            chain_template_sto: str,
            templates_result: templates_custom.TemplateSearchResult,
            monomer_models_result: Mapping[str, templates_custom.TemplateSearchResult],
            custom_complex_msa_pairids: bool,
            is_homomer_or_monomer: bool,
            notemplate: bool) -> pipeline_custom.FeatureDict:
        """Runs the monomer pipeline on a single chain."""
        chain_fasta_str = f'>chain_{chain_id}\n{sequence}\n'
        chain_msa_output_dir = os.path.join(msa_output_dir, chain_id)
        if not os.path.exists(chain_msa_output_dir):
            os.makedirs(chain_msa_output_dir)
        with temp_fasta_file(chain_fasta_str) as chain_fasta_path:
            logging.info('Running monomer pipeline on chain %s: %s', chain_id, description)

            chain_template_features = None
            if notemplate:
                chain_template_features = templates_custom.mk_mock_template(sequence)
            else:
                chain_template_results = None
                if chain_template_sto is not None and os.path.exists(chain_template_sto):
                    # with open(chain_template_sto) as f:
                    #     msa_for_templates = f.read()
                    msa_for_templates = parsers.truncate_stockholm_msa(chain_template_sto,
                                                                       max_sequences=self._max_uniprot_hits)
                    msa_for_templates = parsers.deduplicate_stockholm_msa(msa_for_templates)
                    msa_for_templates = parsers.remove_empty_columns_from_stockholm_msa(msa_for_templates)

                    if self.template_searcher.input_format == 'sto':
                        pdb_templates_result = self.template_searcher.query(msa_for_templates)
                    elif self.template_searcher.input_format == 'a3m':
                        uniref90_msa_as_a3m = parsers.convert_stockholm_to_a3m(msa_for_templates)
                        pdb_templates_result = self.template_searcher.query(uniref90_msa_as_a3m)
                    else:
                        raise ValueError('Unrecognized template input format: '
                                         f'{self.template_searcher.input_format}')

                    pdb_hits_out_path = os.path.join(msa_output_dir, f'pdb_hits.{self.template_searcher.output_format}')
                    with open(pdb_hits_out_path, 'w') as f:
                        f.write(pdb_templates_result)

                    pdb_template_hits = self.template_searcher.get_template_hits(
                        output_string=pdb_templates_result, input_sequence=sequence)

                    chain_template_results = self.template_featurizer.get_templates(
                        query_sequence=sequence,
                        hits=pdb_template_hits)

                elif templates_result is not None:
                    chain_template_results = templates_result

                if monomer_models_result is not None:
                    if chain_template_results is None:
                        chain_template_features = monomer_models_result.features
                    else:
                        chain_template_features = chain_template_results.features
                        for template_feature_name in chain_template_features:
                            chain_template_features[template_feature_name] = []
                            for hit_feature in monomer_models_result.hits_features:
                                # print(monomer_models_result.hits_features)
                                chain_template_features[template_feature_name].append(
                                    hit_feature[template_feature_name])
                            for hit_feature in chain_template_results.hits_features:
                                chain_template_features[template_feature_name].append(
                                    hit_feature[template_feature_name])
                        for name in chain_template_features:
                            chain_template_features[name] = np.stack(chain_template_features[name], axis=0).astype(
                                templates_custom.TEMPLATE_FEATURES[name])
                else:
                    chain_template_features = chain_template_results.features

            if chain_monomer_msa.find('.sto') > 0:
                chain_monomer_msa = pipeline_custom.parsers.parse_stockholm(
                    "".join(open(chain_monomer_msa, "r").readlines()))
            else:
                chain_monomer_msa = pipeline_custom.parsers.parse_a3m("".join(open(chain_monomer_msa, "r").readlines()))

            msa_features = pipeline_custom.make_msa_features(msas=[chain_monomer_msa],
                                                             msa_output_dir=chain_msa_output_dir,
                                                             msa_save_path='monomer_final.a3m',
                                                             filter=False)
            # print(chain_template_features)
            chain_features = {
                **pipeline_custom.make_sequence_features(sequence=sequence,
                                                         description=description,
                                                         num_res=len(sequence)),
                **msa_features,
                **chain_template_features,
            }

            with open(os.path.join(chain_msa_output_dir, 'templates.pkl'), 'wb') as f:
                pickle.dump(chain_template_features, f, protocol=4)

            # We only construct the pairing features if there are 2 or more unique
            # sequences.
            if not is_homomer_or_monomer:
                all_seq_msa_features = self._all_seq_msa_features(chain_multimer_msa, chain_msa_output_dir)
                chain_features.update(all_seq_msa_features)
        return chain_features

    def _all_seq_msa_features(self, chain_multimer_msa, msa_output_dir):
        """Get MSA features for unclustered uniprot, for pairing."""
        # out_path = os.path.join(msa_output_dir, 'uniprot_hits.sto')
        # result = pipeline_custom.run_msa_tool(
        #     self._uniprot_msa_runner, input_fasta_path, out_path, 'sto',
        #     self.use_precomputed_msas)
        # msa = parsers.parse_stockholm(result['sto'])
        if chain_multimer_msa.find('.sto') > 0:
            chain_multimer_msa = pipeline_custom.parsers.parse_stockholm(
                "".join(open(chain_multimer_msa, "r").readlines()))
        else:
            chain_multimer_msa = pipeline_custom.parsers.parse_a3m("".join(open(chain_multimer_msa, "r").readlines()))
        chain_multimer_msa = chain_multimer_msa.truncate(max_seqs=self._max_uniprot_hits)
        all_seq_features = make_msa_features(msas=[chain_multimer_msa], msa_output_dir=msa_output_dir,
                                             msa_save_path='multimer_final.a3m')
        valid_feats = msa_pairing.MSA_FEATURES + (
            'msa_uniprot_accession_identifiers',
            'msa_species_identifiers',
        )
        feats = {f'{k}_all_seq': v for k, v in all_seq_features.items()
                 if k in valid_feats}
        return feats

    def process(self,
                input_fasta_path: str,
                msa_output_dir: str,
                template_output_dir: str,
                custom_inputs: custom_params.CustomizedInputs_Multimer) -> pipeline_custom.FeatureDict:
        """Runs alignment tools on the input sequences and creates features."""
        with open(input_fasta_path) as f:
            input_fasta_str = f.read()
        input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)

        chain_id_map, chain_id_seq_map = _make_chain_id_map(sequences=input_seqs,
                                                            descriptions=input_descs)
        chain_id_map_path = os.path.join(msa_output_dir, 'chain_id_map.json')
        with open(chain_id_map_path, 'w') as f:
            chain_id_map_dict = {chain_id: dataclasses.asdict(fasta_chain)
                                 for chain_id, fasta_chain in chain_id_map.items()}
            json.dump(chain_id_map_dict, f, indent=4, sort_keys=True)

        monomer_models_temp_results = {}
        if len(custom_inputs.monomer_model_paths) > 0:
            monomer_models_temp_results = self.monomer_template_featurizer.get_templates(chain_id_map,
                                                                                         template_output_dir)
            # print(monomer_models_temp_results)

        complex_template_result = None
        if os.path.exists(custom_inputs.temp_struct_csv):
            complex_template_result = self.template_featurizer.get_templates(template_output_dir=template_output_dir,
                                                                             chain_id_map=chain_id_map,
                                                                             hits_file=custom_inputs.temp_struct_csv)

        if os.path.exists(custom_inputs.temp_seq_pair_file):
            complex_template_result = templates_custom.ComplexTemplateSearchResult(monomer_results={}, errors=[],
                                                                                   warnings=[])
            temp_seq_pair_df = pd.read_csv(custom_inputs.temp_seq_pair_file)

            for monomer_count, chain_id in enumerate(chain_id_map):
                monomer_template_hits = parsers.convert_values_to_template_hit(temp_seq_pair_df, monomer_count + 1)

                hit_file = custom_inputs.template_hits_files[chain_id_map[chain_id].description]
                pdb_templates_result = open(hit_file).read()
                pdb_template_hits = self.template_searcher.get_template_hits(
                    output_string=pdb_templates_result, input_sequence=chain_id_map[chain_id].sequence)

                seen_seqs = [hit.name + hit.hit_sequence for hit in monomer_template_hits]

                for hit in pdb_template_hits:
                    #print(hit.name.split()[0] + hit.hit_sequence)
                    if hit.name.split()[0] + hit.hit_sequence not in seen_seqs:
                        hit.index = len(monomer_template_hits)
                        monomer_template_hits += [hit]
                        seen_seqs += [hit.name + hit.hit_sequence]

                print(f"{chain_id_map[chain_id].description}: {len(monomer_template_hits)}")
                complex_template_result.monomer_results[chain_id] = self.template_featurizer.get_templates(
                    query_sequence=chain_id_map[chain_id].sequence,
                    hits=monomer_template_hits,
                    sort=False,
                    multimer=True)

        if len(custom_inputs.monomer_temp_csvs) > 0:
            complex_template_result = templates_custom.ComplexTemplateSearchResult(monomer_results={}, errors=[],
                                                                                   warnings=[])
            for index, chain_id in enumerate(chain_id_map):
                complex_template_result.monomer_results[chain_id] = self.template_featurizer.get_templates(
                    query_sequence=chain_id_map[chain_id].sequence,
                    template_pdb_dir=template_output_dir,
                    hits_file=custom_inputs.monomer_temp_csvs[index])

        all_chain_features = {}
        sequence_features = {}
        is_homomer_or_monomer = len(set(input_seqs)) == 1
        for chain_id, fasta_chain in chain_id_map.items():
            # if fasta_chain.sequence in sequence_features:
            #     all_chain_features[chain_id] = copy.deepcopy(
            #         sequence_features[fasta_chain.sequence])
            #     continue

            chain_template_result = None
            chain_template_sto = None
            monomer_models_result = None
            if chain_id in monomer_models_temp_results:
                monomer_models_result = monomer_models_temp_results[chain_id]

            if complex_template_result is not None:
                chain_template_result = complex_template_result.monomer_results[chain_id]
            elif fasta_chain.description in custom_inputs.template_stos:
                chain_template_sto = custom_inputs.template_stos[fasta_chain.description]

            chain_features = self._process_single_chain(
                chain_id=chain_id,
                sequence=fasta_chain.sequence,
                description=fasta_chain.description,
                chain_monomer_msa=custom_inputs.monomer_a3ms[fasta_chain.description],
                chain_multimer_msa=custom_inputs.multimer_a3ms[fasta_chain.description],
                msa_output_dir=msa_output_dir,
                chain_template_sto=chain_template_sto,
                templates_result=chain_template_result,
                monomer_models_result=monomer_models_result,
                custom_complex_msa_pairids=os.path.exists(custom_inputs.msa_pair_file),
                is_homomer_or_monomer=is_homomer_or_monomer,
                notemplate=custom_inputs.notemplate)

            chain_features = convert_monomer_features(chain_features, chain_id=chain_id)
            all_chain_features[chain_id] = chain_features
            sequence_features[fasta_chain.sequence] = chain_features

        all_chain_features_raw = copy.deepcopy(all_chain_features)

        all_chain_features = add_assembly_features(all_chain_features)

        np_example = feature_processing.pair_and_merge(
            all_chain_features=all_chain_features,
            custom_inputs=custom_inputs,
            max_msa_hits=self._max_uniprot_hits,
            msa_output_dir=msa_output_dir
        )

        # Pad MSA to avoid zero-sized extra_msa.
        np_example = pad_msa(np_example, 512)

        # save paired monomer alignments
        save_paired_msas(chain_id_map, all_chain_features_raw, msa_output_dir)

        return np_example
