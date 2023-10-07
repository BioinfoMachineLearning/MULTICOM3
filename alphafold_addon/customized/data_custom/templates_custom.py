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

"""Functions for getting templates and calculating template features."""
import abc
import dataclasses
import datetime
import functools
import glob
import os
import re
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, List

from absl import logging
from alphafold.common import residue_constants
from alphafold.data_custom import pdb_parsing
from alphafold.data_custom import parsers
from alphafold.data_custom.tools import kalign
import numpy as np
import pandas as pd
import copy


# Internal import (7716).


class Error(Exception):
    """Base class for exceptions."""


class NoChainsError(Error):
    """An error indicating that template mmCIF didn't have any chains."""


class SequenceNotInTemplateError(Error):
    """An error indicating that template mmCIF didn't contain the sequence."""


class NoAtomDataInTemplateError(Error):
    """An error indicating that template mmCIF didn't contain atom positions."""


class TemplateAtomMaskAllZerosError(Error):
    """An error indicating that template mmCIF had all atom positions masked."""


class QueryToTemplateAlignError(Error):
    """An error indicating that the query can't be aligned to the template."""


class CaDistanceError(Error):
    """An error indicating that a CA atom distance exceeds a threshold."""


class MultipleChainsError(Error):
    """An error indicating that multiple chains were found for a given ID."""


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


TEMPLATE_FEATURES = {
    'template_aatype': np.float32,
    'template_all_atom_masks': np.float32,
    'template_all_atom_positions': np.float32,
    'template_domain_names': np.object,
    'template_sequence': np.object,
    'template_sum_probs': np.float32,
}


@dataclasses.dataclass(frozen=True)
class FastaChain:
    sequence: str
    description: str


@dataclasses.dataclass(frozen=True)
class CustomizedTemplateHit:
    query_name: str
    template_name: str
    template_chain: str
    aligned_length: int
    aln_temp: str
    tstart: int
    tend: int
    aln_query: str
    qstart: int
    qend: int
    from_predicted_structure: bool


@dataclasses.dataclass(frozen=True)
class CustomizedComplexTemplateHit:
    """Class representing a template hit."""
    template_name: str
    template_sequence: str
    monomer_hits: Mapping[str, CustomizedTemplateHit]


def _assess_template_hit(hit: CustomizedTemplateHit,
                         query_sequence: str,
                         max_subsequence_ratio: float = 0.95,
                         min_align_ratio: float = 0.1) -> bool:
    """Determines if template is valid (without parsing the template mmcif file).

  Args:
    hit: HhrHit for the template.
    hit_pdb_code: The 4 letter pdb code of the template hit. This might be
      different from the value in the actual hit since the original pdb might
      have become obsolete.
    query_sequence: Amino acid sequence of the query.
    release_dates: Dictionary mapping pdb codes to their structure release
      dates.
    release_date_cutoff: Max release date that is valid for this query.
    max_subsequence_ratio: Exclude any exact matches with this much overlap.
    min_align_ratio: Minimum overlap between the template and query.

  Returns:
    True if the hit passed the prefilter. Raises an exception otherwise.

  Raises:
    DateError: If the hit date was after the max allowed date.
    AlignRatioError: If the hit align ratio to the query was too small.
    DuplicateError: If the hit was an exact subsequence of the query.
    LengthError: If the hit was too short.
  """
    align_ratio = hit.aligned_length / len(query_sequence)

    template_sequence = hit.aln_temp.replace('-', '')
    length_ratio = float(len(template_sequence)) / len(query_sequence)

    # Check whether the template is a large subsequence or duplicate of original
    # query. This can happen due to duplicate entries in the PDB database.
    duplicate = (template_sequence in query_sequence and
                 length_ratio > max_subsequence_ratio)

    if hit.from_predicted_structure:
        duplicate = False

    if align_ratio <= min_align_ratio:
        raise AlignRatioError('Proportion of residues aligned to query too small. '
                              f'Align ratio: {align_ratio}.')

    if duplicate:
        raise DuplicateError('Template is an exact subsequence of query with large '
                             f'coverage. Length ratio: {length_ratio}.')

    if len(template_sequence) < 10:
        raise LengthError(f'Template too short. Length: {len(template_sequence)}.')

    return True


def _find_template_in_pdb(
        template_chain_id: str,
        template_sequence: str,
        pdb_object: pdb_parsing.PdbObject) -> Tuple[str, str, int]:
    """Tries to find the template chain in the given pdb file.

  This method tries the three following things in order:
    1. Tries if there is an exact match in both the chain ID and the sequence.
       If yes, the chain sequence is returned. Otherwise:
    2. Tries if there is an exact match only in the sequence.
       If yes, the chain sequence is returned. Otherwise:
    3. Tries if there is a fuzzy match (X = wildcard) in the sequence.
       If yes, the chain sequence is returned.
  If none of these succeed, a SequenceNotInTemplateError is thrown.

  Args:
    template_chain_id: The template chain ID.
    template_sequence: The template chain sequence.
    pdb_object: The PDB object to search for the template in.

  Returns:
    A tuple with:
    * The chain sequence that was found to match the template in the PDB object.
    * The ID of the chain that is being returned.
    * The offset where the template sequence starts in the chain sequence.

  Raises:
    SequenceNotInTemplateError: If no match is found after the steps described
      above.
  """
    # Try if there is an exact match in both the chain ID and the (sub)sequence.
    pdb_id = pdb_object.file_id
    chain_sequence = pdb_object.chain_to_seqres.get(template_chain_id)
    if chain_sequence and (template_sequence in chain_sequence):
        logging.info(
            'Found an exact template match %s.', pdb_id)
        mapping_offset = chain_sequence.find(template_sequence)
        return chain_sequence, template_chain_id, mapping_offset

    # Try if there is an exact match in the (sub)sequence only.
    for chain_id, chain_sequence in pdb_object.chain_to_seqres.items():
        if chain_sequence and (template_sequence in chain_sequence):
            logging.info('Found a sequence-only match %s.', pdb_id)
            mapping_offset = chain_sequence.find(template_sequence)
            return chain_sequence, chain_id, mapping_offset

    # Return a chain sequence that fuzzy matches (X = wildcard) the template.
    # Make parentheses unnamed groups (?:_) to avoid the 100 named groups limit.
    regex = ['.' if aa == 'X' else '(?:%s|X)' % aa for aa in template_sequence]
    regex = re.compile(''.join(regex))
    for chain_id, chain_sequence in pdb_object.chain_to_seqres.items():
        match = re.search(regex, chain_sequence)
        if match:
            logging.info('Found a fuzzy sequence-only match %s.', pdb_id)
            mapping_offset = match.start()
            return chain_sequence, chain_id, mapping_offset

    # No hits, raise an error.
    raise SequenceNotInTemplateError(
        'Could not find the template sequence in %s. Template sequence: %s, '
        'chain_to_seqres: %s' % (pdb_id, template_sequence,
                                 pdb_object.chain_to_seqres))


def _realign_pdb_template_to_query(
        old_template_sequence: str,
        template_chain_id: str,
        pdb_object: pdb_parsing.PdbObject,
        old_mapping: Mapping[int, int],
        kalign_binary_path: str) -> Tuple[str, Mapping[int, int]]:
    """Aligns template from the pdb_object to the query.

  In case PDB70 contains a different version of the template sequence, we need
  to perform a realignment to the actual sequence that is in the mmCIF file.
  This method performs such realignment, but returns the new sequence and
  mapping only if the sequence in the mmCIF file is 90% identical to the old
  sequence.

  Note that the old_template_sequence comes from the hit, and contains only that
  part of the chain that matches with the query while the new_template_sequence
  is the full chain.

  Args:
    old_template_sequence: The template sequence that was returned by the PDB
      template search (typically done using HHSearch).
    template_chain_id: The template chain id was returned by the PDB template
      search (typically done using HHSearch). This is used to find the right
      chain in the pdb_object chain_to_seqres mapping.
    pdb_object: A pdb_object which holds the actual template data.
    old_mapping: A mapping from the query sequence to the template sequence.
      This mapping will be used to compute the new mapping from the query
      sequence to the actual pdb_object template sequence by aligning the
      old_template_sequence and the actual template sequence.
    kalign_binary_path: The path to a kalign executable.

  Returns:
    A tuple (new_template_sequence, new_query_to_template_mapping) where:
    * new_template_sequence is the actual template sequence that was found in
      the pdb_object.
    * new_query_to_template_mapping is the new mapping from the query to the
      actual template found in the pdb_object.

  Raises:
    QueryToTemplateAlignError:
    * If there was an error thrown by the alignment tool.
    * Or if the actual template sequence differs by more than 10% from the
      old_template_sequence.
  """
    aligner = kalign.Kalign(binary_path=kalign_binary_path)
    new_template_sequence = pdb_object.chain_to_seqres.get(
        template_chain_id, '')

    # Sometimes the template chain id is unknown. But if there is only a single
    # sequence within the pdb_object, it is safe to assume it is that one.
    if not new_template_sequence:
        if len(pdb_object.chain_to_seqres) == 1:
            logging.info('Could not find %s in %s, but there is only 1 sequence, so '
                         'using that one.',
                         template_chain_id,
                         pdb_object.file_id)
            new_template_sequence = list(pdb_object.chain_to_seqres.values())[0]
        else:
            raise QueryToTemplateAlignError(
                f'Could not find chain {template_chain_id} in {pdb_object.file_id}. '
                'If there are no mmCIF parsing errors, it is possible it was not a '
                'protein chain.')

    try:
        parsed_a3m = parsers.parse_a3m(
            aligner.align([old_template_sequence, new_template_sequence]))
        old_aligned_template, new_aligned_template = parsed_a3m.sequences
    except Exception as e:
        raise QueryToTemplateAlignError(
            'Could not align old template %s to template %s (%s_%s). Error: %s' %
            (old_template_sequence, new_template_sequence, pdb_object.file_id,
             template_chain_id, str(e)))

    logging.info('Old aligned template: %s\nNew aligned template: %s',
                 old_aligned_template, new_aligned_template)

    old_to_new_template_mapping = {}
    old_template_index = -1
    new_template_index = -1
    num_same = 0
    for old_template_aa, new_template_aa in zip(
            old_aligned_template, new_aligned_template):
        if old_template_aa != '-':
            old_template_index += 1
        if new_template_aa != '-':
            new_template_index += 1
        if old_template_aa != '-' and new_template_aa != '-':
            old_to_new_template_mapping[old_template_index] = new_template_index
            if old_template_aa == new_template_aa:
                num_same += 1

    # Require at least 90 % sequence identity wrt to the shorter of the sequences.
    if float(num_same) / min(
            len(old_template_sequence), len(new_template_sequence)) < 0.9:
        raise QueryToTemplateAlignError(
            'Insufficient similarity of the sequence in the database: %s to the '
            'actual sequence in the mmCIF file %s_%s: %s. We require at least '
            '90 %% similarity wrt to the shorter of the sequences. This is not a '
            'problem unless you think this is a template that should be included.' %
            (old_template_sequence, pdb_object.file_id, template_chain_id,
             new_template_sequence))

    new_query_to_template_mapping = {}
    for query_index, old_template_index in old_mapping.items():
        new_query_to_template_mapping[query_index] = (
            old_to_new_template_mapping.get(old_template_index, -1))

    new_template_sequence = new_template_sequence.replace('-', '')

    return new_template_sequence, new_query_to_template_mapping


def _check_residue_distances(all_positions: np.ndarray,
                             all_positions_mask: np.ndarray,
                             max_ca_ca_distance: float):
    """Checks if the distance between unmasked neighbor residues is ok."""
    ca_position = residue_constants.atom_order['CA']
    prev_is_unmasked = False
    prev_calpha = None
    for i, (coords, mask) in enumerate(zip(all_positions, all_positions_mask)):
        this_is_unmasked = bool(mask[ca_position])
        if this_is_unmasked:
            this_calpha = coords[ca_position]
            if prev_is_unmasked:
                distance = np.linalg.norm(this_calpha - prev_calpha)
                if distance > max_ca_ca_distance:
                    raise CaDistanceError(
                        'The distance between residues %d and %d is %f > limit %f.' % (
                            i, i + 1, distance, max_ca_ca_distance))
            prev_calpha = this_calpha
        prev_is_unmasked = this_is_unmasked


def _get_atom_positions(
        pdb_object: pdb_parsing.PdbObject,
        auth_chain_id: str,
        max_ca_ca_distance: float) -> Tuple[np.ndarray, np.ndarray]:
    """Gets atom positions and mask from a list of Biopython Residues."""
    num_res = len(pdb_object.chain_to_seqres[auth_chain_id])

    relevant_chains = [c for c in pdb_object.structure.get_chains()
                       if c.id == auth_chain_id]

    if len(relevant_chains) != 1:
        raise MultipleChainsError(
            f'Expected exactly one chain in structure with id {auth_chain_id}.')
    chain = relevant_chains[0]

    all_positions = np.zeros([num_res, residue_constants.atom_type_num, 3])
    all_positions_mask = np.zeros([num_res, residue_constants.atom_type_num],
                                  dtype=np.int64)
    for res_index in range(num_res):
        pos = np.zeros([residue_constants.atom_type_num, 3], dtype=np.float32)
        mask = np.zeros([residue_constants.atom_type_num], dtype=np.float32)
        res_at_position = pdb_object.seqres_to_structure[auth_chain_id][res_index]
        if not res_at_position.is_missing:
            res = chain[(res_at_position.hetflag,
                         res_at_position.position.residue_number,
                         res_at_position.position.insertion_code)]
            for atom in res.get_atoms():
                atom_name = atom.get_name()
                x, y, z = atom.get_coord()
                if atom_name in residue_constants.atom_order.keys():
                    pos[residue_constants.atom_order[atom_name]] = [x, y, z]
                    mask[residue_constants.atom_order[atom_name]] = 1.0
                elif atom_name.upper() == 'SE' and res.get_resname() == 'MSE':
                    # Put the coordinates of the selenium atom in the sulphur column.
                    pos[residue_constants.atom_order['SD']] = [x, y, z]
                    mask[residue_constants.atom_order['SD']] = 1.0

            # Fix naming errors in arginine residues where NH2 is incorrectly
            # assigned to be closer to CD than NH1.
            cd = residue_constants.atom_order['CD']
            nh1 = residue_constants.atom_order['NH1']
            nh2 = residue_constants.atom_order['NH2']
            if (res.get_resname() == 'ARG' and
                    all(mask[atom_index] for atom_index in (cd, nh1, nh2)) and
                    (np.linalg.norm(pos[nh1] - pos[cd]) >
                     np.linalg.norm(pos[nh2] - pos[cd]))):
                pos[nh1], pos[nh2] = pos[nh2].copy(), pos[nh1].copy()
                mask[nh1], mask[nh2] = mask[nh2].copy(), mask[nh1].copy()

        all_positions[res_index] = pos
        all_positions_mask[res_index] = mask
    _check_residue_distances(
        all_positions, all_positions_mask, max_ca_ca_distance)
    return all_positions, all_positions_mask


def _extract_template_features(
        pdb_object: pdb_parsing.PdbObject,
        pdb_id: str,
        mapping: Mapping[int, int],
        template_sequence: str,
        query_sequence: str,
        template_chain_id: str,
        kalign_binary_path: str) -> Tuple[Dict[str, Any], Optional[str]]:
    """Parses atom positions in the target structure and aligns with the query.

  Atoms for each residue in the template structure are indexed to coincide
  with their corresponding residue in the query sequence, according to the
  alignment mapping provided.

  """
    if pdb_object is None or not pdb_object.chain_to_seqres:
        raise NoChainsError('No chains in PDB: %s_%s' % (pdb_id, template_chain_id))

    warning = None
    try:
        seqres, chain_id, mapping_offset = _find_template_in_pdb(
            template_chain_id=template_chain_id,
            template_sequence=template_sequence,
            pdb_object=pdb_object)
    except SequenceNotInTemplateError:
        # If PDB70 contains a different version of the template, we use the sequence
        # from the pdb_object.
        chain_id = template_chain_id
        warning = (
            f'The exact sequence {template_sequence} was not found in '
            f'{pdb_id}_{chain_id}. Realigning the template to the actual sequence.')
        logging.warning(warning)
        # This throws an exception if it fails to realign the hit.
        seqres, mapping = _realign_pdb_template_to_query(
            old_template_sequence=template_sequence,
            template_chain_id=template_chain_id,
            pdb_object=pdb_object,
            old_mapping=mapping,
            kalign_binary_path=kalign_binary_path)
        logging.info('Sequence in %s_%s: %s successfully realigned to %s',
                     pdb_id, chain_id, template_sequence, seqres)
        # The template sequence changed.
        template_sequence = seqres
        # No mapping offset, the query is aligned to the actual sequence.
        mapping_offset = 0

    try:
        # Essentially set to infinity - we don't want to reject templates unless
        # they're really really bad.
        all_atom_positions, all_atom_mask = _get_atom_positions(
            pdb_object, chain_id, max_ca_ca_distance=150.0)
    except (CaDistanceError, KeyError) as ex:
        raise NoAtomDataInTemplateError(
            'Could not get atom data (%s_%s): %s' % (pdb_id, chain_id, str(ex))
        ) from ex

    all_atom_positions = np.split(all_atom_positions, all_atom_positions.shape[0])
    all_atom_masks = np.split(all_atom_mask, all_atom_mask.shape[0])

    output_templates_sequence = []
    templates_all_atom_positions = []
    templates_all_atom_masks = []

    for _ in query_sequence:
        # Residues in the query_sequence that are not in the template_sequence:
        templates_all_atom_positions.append(
            np.zeros((residue_constants.atom_type_num, 3)))
        templates_all_atom_masks.append(np.zeros(residue_constants.atom_type_num))
        output_templates_sequence.append('-')

    for k, v in mapping.items():
        template_index = v + mapping_offset
        templates_all_atom_positions[k] = all_atom_positions[template_index][0]
        templates_all_atom_masks[k] = all_atom_masks[template_index][0]
        output_templates_sequence[k] = template_sequence[v]

    # Alanine (AA with the lowest number of atoms) has 5 atoms (C, CA, CB, N, O).
    if np.sum(templates_all_atom_masks) < 5:
        raise TemplateAtomMaskAllZerosError(
            'Template all atom mask was all zeros: %s_%s. Residue range: %d-%d' %
            (pdb_id, chain_id, min(mapping.values()) + mapping_offset,
             max(mapping.values()) + mapping_offset))

    output_templates_sequence = ''.join(output_templates_sequence)
    templates_aatype = residue_constants.sequence_to_onehot(
        output_templates_sequence, residue_constants.HHBLITS_AA_TO_ID)

    return (
        {
            'template_all_atom_positions': np.array(templates_all_atom_positions),
            'template_all_atom_masks': np.array(templates_all_atom_masks),
            'template_sequence': output_templates_sequence.encode(),
            'template_aatype': np.array(templates_aatype),
            'template_domain_names': f'{pdb_id.lower()}_{chain_id}'.encode(),
        },
        warning)


def _build_query_to_hit_index_mapping(
        hit_query_sequence: str,
        hit_sequence: str,
        indices_hit: Sequence[int],
        indices_query: Sequence[int],
        original_query_sequence: str) -> Mapping[int, int]:
    """Gets mapping from indices in original query sequence to indices in the hit.

  hit_query_sequence and hit_sequence are two aligned sequences containing gap
  characters. hit_query_sequence contains only the part of the original query
  sequence that matched the hit. When interpreting the indices from the .hhr, we
  need to correct for this to recover a mapping from original query sequence to
  the hit sequence.

  Args:
    hit_query_sequence: The portion of the query sequence that is in the .hhr
      hit
    hit_sequence: The portion of the hit sequence that is in the .hhr
    indices_hit: The indices for each aminoacid relative to the hit sequence
    indices_query: The indices for each aminoacid relative to the original query
      sequence
    original_query_sequence: String describing the original query sequence.

  Returns:
    Dictionary with indices in the original query sequence as keys and indices
    in the hit sequence as values.
  """
    # If the hit is empty (no aligned residues), return empty mapping
    if not hit_query_sequence:
        return {}

    # Remove gaps and find the offset of hit.query relative to original query.
    hhsearch_query_sequence = hit_query_sequence.replace('-', '')
    hit_sequence = hit_sequence.replace('-', '')
    hhsearch_query_offset = original_query_sequence.find(hhsearch_query_sequence)

    # Index of -1 used for gap characters. Subtract the min index ignoring gaps.
    min_idx = min(x for x in indices_hit if x > -1)
    fixed_indices_hit = [
        x - min_idx if x > -1 else -1 for x in indices_hit
    ]

    min_idx = min(x for x in indices_query if x > -1)
    fixed_indices_query = [x - min_idx if x > -1 else -1 for x in indices_query]

    # Zip the corrected indices, ignore case where both seqs have gap characters.
    mapping = {}
    for q_i, q_t in zip(fixed_indices_query, fixed_indices_hit):
        if q_t != -1 and q_i != -1:
            if (q_t >= len(hit_sequence) or
                    q_i + hhsearch_query_offset >= len(original_query_sequence)):
                continue
            mapping[q_i + hhsearch_query_offset] = q_t

    return mapping


@dataclasses.dataclass(frozen=True)
class SingleHitResult:
    features: Optional[Mapping[str, Any]]
    error: Optional[str]
    warning: Optional[str]


@dataclasses.dataclass(frozen=True)
class SingleComplexHitResult:
    monomer_hit_results: Mapping[str, SingleHitResult]
    error: Optional[str]
    warning: Optional[str]


@functools.lru_cache(16, typed=False)
def _read_file(path):
    with open(path, 'r') as f:
        file_data = f.read()
    return file_data


def _get_indices(sequence: str, start: int) -> List[int]:
    """Returns indices for non-gap/insert residues starting at the given index."""
    indices = []
    counter = start
    for symbol in sequence:
        # Skip gaps but add a placeholder so that the alignment is preserved.
        if symbol == '-':
            indices.append(-1)
        # Skip deleted residues, but increase the counter.
        elif symbol.islower():
            counter += 1
        # Normal aligned residue. Increase the counter and append to indices.
        else:
            indices.append(counter)
            counter += 1
    return indices


def _process_single_hit(
        query_sequence: str,
        monomer_hit: CustomizedTemplateHit,
        atom_dir: str,
        kalign_binary_path: str,
        strict_error_check: bool = False) -> SingleHitResult:
    """Tries to extract template features from a single HHSearch hit."""

    try:
        _assess_template_hit(hit=monomer_hit, query_sequence=query_sequence)
    except PrefilterError as e:
        msg = f"hit {monomer_hit.template_name} did not pass prefilter: {str(e)}"
        logging.info(msg)
        if strict_error_check and isinstance(e, (DateError, DuplicateError)):
            # In strict mode we treat some prefilter cases as errors.
            return SingleHitResult(features=None, error=msg, warning=None)

        return SingleHitResult(features=None, error=None, warning=None)

    mapping = _build_query_to_hit_index_mapping(monomer_hit.aln_query, monomer_hit.aln_temp,
                                                _get_indices(monomer_hit.aln_temp, start=monomer_hit.tstart - 1),
                                                _get_indices(monomer_hit.aln_query, start=monomer_hit.qstart - 1),
                                                query_sequence)

    # The mapping is from the query to the actual hit sequence, so we need to
    # remove gaps (which regardless have a missing confidence score).
    template_sequence = monomer_hit.aln_temp.replace('-', '')

    pdb_path = os.path.join(atom_dir, monomer_hit.template_name + '.atom')

    logging.debug('Reading PDB entry from %s. Query: %s, template: %s', pdb_path,
                  query_sequence, template_sequence)
    # Fail if we can't find the mmCIF file.
    pdb_string = _read_file(pdb_path)

    parsing_result = pdb_parsing.parse(file_id=monomer_hit.template_name,
                                       chain_id=monomer_hit.template_chain,
                                       pdb_string=pdb_string)

    try:
        features, realign_warning = _extract_template_features(
            pdb_object=parsing_result.pdb_object,
            pdb_id=monomer_hit.template_name,
            mapping=mapping,
            template_sequence=template_sequence,
            query_sequence=query_sequence,
            template_chain_id=monomer_hit.template_chain,
            kalign_binary_path=kalign_binary_path)

        features['template_sum_probs'] = [0]

        # It is possible there were some errors when parsing the other chains in the
        # mmCIF file, but the template features for the chain we want were still
        # computed. In such case the mmCIF parsing errors are not relevant.
        return SingleHitResult(
            features=features, error=None, warning=realign_warning)
    except (NoChainsError, NoAtomDataInTemplateError,
            TemplateAtomMaskAllZerosError) as e:
        # These 3 errors indicate missing mmCIF experimental data rather than a
        # problem with the template search, so turn them into warnings.
        warning = ('%s_%s: feature extracting errors: '
                   '%s, pdb parsing errors: %s'
                   % (monomer_hit.template_name, monomer_hit.template_chain,
                      str(e), parsing_result.errors))
        if strict_error_check:
            return SingleHitResult(features=None, error=warning, warning=None)
        else:
            return SingleHitResult(features=None, error=None, warning=warning)
    except Error as e:
        error = ('%s_%s: feature extracting errors: '
                 '%s, PDB parsing errors: %s'
                 % (monomer_hit.template_name, monomer_hit.template_chain,
                    str(e), parsing_result.errors))
        return SingleHitResult(features=None, error=error, warning=None)


def _process_single_complex_hit(
        complex_hit: CustomizedComplexTemplateHit,
        atom_dir: str,
        chain_id_map: str,
        kalign_binary_path: str,
        strict_error_check: bool = False) -> SingleComplexHitResult:
    errors = []
    warnings = []
    monomer_hit_results = {}

    for monomer_id in complex_hit.monomer_hits:
        monomer_result = _process_single_hit(
            query_sequence=chain_id_map[monomer_id].sequence,
            monomer_hit=complex_hit.monomer_hits[monomer_id],
            atom_dir=atom_dir,
            kalign_binary_path=kalign_binary_path,
            strict_error_check=strict_error_check)
        if monomer_result.error:
            print(monomer_result.error)
            errors.append(monomer_result.error)

        # There could be an error even if there are some results, e.g. thrown by
        # other unparsable chains in the same mmCIF file.
        if monomer_result.warning:
            print(monomer_result.warning)
            warnings.append(monomer_result.warning)

        if monomer_result.features is None:
            logging.debug('Skipped invalid hit %s, error: %s, warning: %s',
                          complex_hit.monomer_hits[monomer_id].template_name, monomer_result.error,
                          monomer_result.warning)

        monomer_hit_results[monomer_id] = monomer_result.features

        # print(monomer_result.features['template_all_atom_positions'])

    return SingleComplexHitResult(monomer_hit_results=monomer_hit_results,
                                  error='\n'.join(errors),
                                  warning='\n'.join(warnings))


@dataclasses.dataclass(frozen=True)
class TemplateSearchResult:
    features: Mapping[str, Any]
    hits_features: Sequence[Any]
    errors: Sequence[str]
    warnings: Sequence[str]


@dataclasses.dataclass(frozen=True)
class ComplexTemplateSearchResult:
    monomer_results: Mapping[str, TemplateSearchResult]
    errors: Sequence[str]
    warnings: Sequence[str]


def add_chain_info_to_atom_file(infile, chainid, outfile):
    with open(outfile, 'w') as out:
        for line in open(infile):
            if not line.startswith('ATOM'):
                out.write(line)
                continue
            out.write(line[:21] + chainid + line[22:])


class CustomizedComplexHitFeaturizer:
    """A class for turning a3m hits from hmmsearch to template features."""

    def __init__(
            self,
            input_atom_dir: str,
            # max_template_date: str,
            max_hits: int,
            kalign_binary_path: str,
            # release_dates_path: Optional[str],
            # obsolete_pdbs_path: Optional[str],
            strict_error_check: bool = False):

        self._input_atom_dir = input_atom_dir
        if not glob.glob(os.path.join(self._input_atom_dir, '*.atom')):
            logging.error('Could not find atom files in %s', self._input_atom_dir)
            raise ValueError(f'Could not find atom files in {self._input_atom_dir}')

        self._max_hits = max_hits
        self._kalign_binary_path = kalign_binary_path
        self._strict_error_check = strict_error_check

    def get_templates(
            self,
            chain_id_map: Mapping[str, FastaChain],
            template_output_dir: str,
            hits_file: str) -> TemplateSearchResult:
        """Computes the templates for given query sequence (more details above)."""

        if not os.path.exists(template_output_dir):
            os.makedirs(template_output_dir)

        pdb_hits_pd = pd.read_csv(hits_file)
        # pdb_hits_pd_sorted = pdb_hits_pd.sort_values(by=['avg_tmscore'], ascending=False, ignore_index=True)

        all_chain_template_features = {}
        for chainid in chain_id_map:
            template_features = {}
            for template_feature_name in TEMPLATE_FEATURES:
                template_features[template_feature_name] = []
            all_chain_template_features[chainid] = template_features

        already_seen = []
        errors = []
        warnings = []
        all_chain_hits_features = {}

        sequences = []
        for chainid in chain_id_map:
            sequences += [chain_id_map[chainid].sequence]
            all_chain_hits_features[chainid] = []
        logging.info('Searching for template for: %s', ','.join(sequences))

        for index, row in pdb_hits_pd.iterrows():
            # We got all the templates we wanted, stop processing hits.
            if len(already_seen) >= self._max_hits:
                break

            monomer_hits = {}
            complex_hit_query_names = []
            complex_hit_template_names = []
            complex_hit_template_sequences = []
            for i, chainid in enumerate(chain_id_map):
                hit = CustomizedTemplateHit(query_name=chainid,
                                            template_name=row[f'template{i + 1}'],
                                            template_chain=row[f'template{i + 1}'][4],
                                            aligned_length=int(row[f'aligned_length{i + 1}']),
                                            aln_temp=row[f'aln_temp{i + 1}'],
                                            tstart=int(row[f'tstart{i + 1}']),
                                            tend=int(row[f'tend{i + 1}']),
                                            aln_query=row[f'aln_query{i + 1}'],
                                            qstart=int(row[f'qstart{i + 1}']),
                                            qend=int(row[f'qend{i + 1}']),
                                            from_predicted_structure=False)
                monomer_hits[chainid] = hit
                complex_hit_query_names += [chainid]
                complex_hit_template_names += [row[f'template{i + 1}']]
                complex_hit_template_sequences += [row[f'aln_temp{i + 1}'].replace('_', '')]

                ori_atom_file = os.path.join(self._input_atom_dir, hit.template_name + '.atom')
                if not os.path.exists(ori_atom_file):
                    ori_atom_file = os.path.join(self._input_atom_dir, hit.template_name.replace('.gz', ''))

                # add chain info to the original atom file
                add_chain_info_to_atom_file(infile=ori_atom_file,
                                            chainid=hit.template_chain,
                                            outfile=os.path.join(template_output_dir, hit.template_name + '.atom'))

            complex_hit = CustomizedComplexTemplateHit(template_name='_'.join(complex_hit_template_names),
                                                       template_sequence=''.join(complex_hit_template_sequences),
                                                       monomer_hits=monomer_hits)

            logging.info('Processing complex templates: %s', complex_hit.template_name)

            complex_result = _process_single_complex_hit(
                complex_hit=complex_hit,
                atom_dir=template_output_dir,
                chain_id_map=chain_id_map,
                kalign_binary_path=self._kalign_binary_path,
                strict_error_check=self._strict_error_check)

            if complex_result.error:
                errors.append(complex_result.error)

            # There could be an error even if there are some results, e.g. thrown by
            # other unparsable chains in the same mmCIF file.
            if complex_result.warning:
                warnings.append(complex_result.warning)

            # print(chain_id_map)
            # print(complex_result.monomer_hit_results)
            valid_complex_hit = True
            for chainid in chain_id_map:
                if complex_result.monomer_hit_results[chainid] is None:
                    logging.debug('Skipped invalid hit %s, error: %s, warning: %s',
                                  complex_hit.template_name, complex_result.error, complex_result.warning)
                    valid_complex_hit = False

            if valid_complex_hit:
                already_seen_key = complex_hit.template_sequence
                if already_seen_key in already_seen:
                    continue
                # Increment the hit counter, since we got features out of this hit.
                already_seen += [already_seen_key]

                for chainid in complex_result.monomer_hit_results:
                    for k in all_chain_template_features[chainid]:
                        all_chain_template_features[chainid][k].append(complex_result.monomer_hit_results[chainid][k])
                    all_chain_hits_features[chainid] += [complex_result.monomer_hit_results[chainid]]

        # print(already_seen)
        monomer_results = {}
        if already_seen:
            for chainid in complex_result.monomer_hit_results:
                for name in all_chain_template_features[chainid]:
                    all_chain_template_features[chainid][name] = np.stack(
                        all_chain_template_features[chainid][name], axis=0).astype(TEMPLATE_FEATURES[name])
                monomer_results[chainid] = TemplateSearchResult(features=all_chain_template_features[chainid],
                                                                hits_features=all_chain_hits_features[chainid],
                                                                errors=errors, warnings=warnings)
        else:
            for chainid in chain_id_map:
                num_res = len(chain_id_map[chainid].sequence)
                # Construct a default template with all zeros.
                template_features = {
                    'template_aatype': np.zeros(
                        (1, num_res, len(residue_constants.restypes_with_x_and_gap)),
                        np.float32),
                    'template_all_atom_masks': np.zeros(
                        (1, num_res, residue_constants.atom_type_num), np.float32),
                    'template_all_atom_positions': np.zeros(
                        (1, num_res, residue_constants.atom_type_num, 3), np.float32),
                    'template_domain_names': np.array([''.encode()], dtype=np.object),
                    'template_sequence': np.array([''.encode()], dtype=np.object),
                    'template_sum_probs': np.array([0], dtype=np.float32)
                }
                all_chain_template_features[chainid] = TemplateSearchResult(features=template_features,
                                                                            hits_features=[],
                                                                            errors=errors, warnings=warnings)

        return ComplexTemplateSearchResult(
            monomer_results=monomer_results, errors=errors, warnings=warnings)


class CustomizedComplexMonomerHitFeaturizer:
    """A class for turning a3m hits from hmmsearch to template features."""

    def __init__(
            self,
            kalign_binary_path: str,
            monomer_model_paths: List[str],
            strict_error_check: bool = False,
            template_count: int = 5):

        self._kalign_binary_path = kalign_binary_path
        self._strict_error_check = strict_error_check
        self._monomer_model_paths = monomer_model_paths
        self._template_count = template_count

    def get_templates(
            self,
            chain_id_map: Mapping[str, FastaChain],
            template_output_dir: str) -> Mapping[str, TemplateSearchResult]:
        """Computes the templates for given query sequence (more details above)."""

        if not os.path.exists(template_output_dir):
            os.makedirs(template_output_dir)

        all_chain_template_features = {}
        errors = []
        warnings = []

        sequences = []
        for chainid in chain_id_map:
            sequences += [chain_id_map[chainid].sequence]
        logging.info('Searching for template for: %s', ','.join(sequences))

        if self._monomer_model_paths is not None:
            for i, chainid in enumerate(chain_id_map):
                template_features = {}
                for template_feature_name in TEMPLATE_FEATURES:
                    template_features[template_feature_name] = []

                hits_features = []
                for n in range(0, self._template_count):
                    hit = CustomizedTemplateHit(query_name=chainid,
                                                template_name=chainid + f'_{n}',
                                                template_chain=chainid,
                                                aligned_length=len(chain_id_map[chainid].sequence),
                                                aln_temp=chain_id_map[chainid].sequence,
                                                tstart=1,
                                                tend=len(chain_id_map[chainid].sequence),
                                                aln_query=chain_id_map[chainid].sequence,
                                                qstart=1,
                                                qend=len(chain_id_map[chainid].sequence),
                                                from_predicted_structure=True)

                    add_chain_info_to_atom_file(infile=os.path.join(self._monomer_model_paths[i], f'ranked_{n}.pdb'),
                                                chainid=hit.template_chain,
                                                outfile=os.path.join(template_output_dir, hit.template_name + '.atom'))

                    result = _process_single_hit(query_sequence=chain_id_map[chainid].sequence,
                                                 monomer_hit=hit,
                                                 atom_dir=template_output_dir,
                                                 kalign_binary_path=self._kalign_binary_path,
                                                 strict_error_check=self._strict_error_check)

                    if result.error:
                        errors.append(result.error)

                    # There could be an error even if there are some results, e.g. thrown by
                    # other unparsable chains in the same mmCIF file.
                    if result.warning:
                        warnings.append(result.warning)

                    if result.features is None:
                        logging.info('Skipped invalid hit %s, error: %s, warning: %s',
                                     hit.template_name, result.error, result.warning)
                    else:
                        for k in template_features:
                            template_features[k].append(result.features[k])
                        hits_features += [result.features]

                for name in template_features:
                    template_features[name] = np.stack(
                        template_features[name], axis=0).astype(TEMPLATE_FEATURES[name])

                all_chain_template_features[chainid] = TemplateSearchResult(features=template_features,
                                                                            hits_features=hits_features,
                                                                            errors=errors, warnings=warnings)

        return all_chain_template_features


class CustomizedMonomerHitFeaturizer:
    """A class for turning a3m hits from hmmsearch to template features."""

    def __init__(
            self,
            input_pdb_dir: str,
            max_hits: int,
            kalign_binary_path: str,
            strict_error_check: bool = False):

        self._kalign_binary_path = kalign_binary_path
        self._strict_error_check = strict_error_check
        self._max_hits = max_hits
        self._input_pdb_dir = input_pdb_dir

    def get_templates(self,
                      query_sequence: str,
                      template_pdb_dir: str,
                      hits_file: str) -> TemplateSearchResult:

        """Computes the templates for given query sequence (more details above)."""
        logging.info('Searching for template for: %s', query_sequence)

        template_features = {}
        for template_feature_name in TEMPLATE_FEATURES:
            template_features[template_feature_name] = []

        num_hits = 0
        errors = []
        warnings = []
        indices_map = []
        hits_features = []

        pdb_hits_pd = pd.read_csv(hits_file, sep='\t')

        for i in range(pdb_hits_pd.shape[0]):
            if num_hits >= self._max_hits:
                break

            template_chain = 'A'
            from_pdb = pdb_hits_pd.loc[i, 'target'].find('.atom.gz') > 0
            if from_pdb:
                template_chain = pdb_hits_pd.loc[i, 'target'][4]

            hit = CustomizedTemplateHit(query_name='A',
                                        template_name=pdb_hits_pd.loc[i, 'target'].split('.')[0],
                                        template_chain=template_chain,
                                        aligned_length=int(pdb_hits_pd.loc[i, 'alnlen']),
                                        aln_temp=pdb_hits_pd.loc[i, 'taln'],
                                        tstart=int(pdb_hits_pd.loc[i, 'tstart']),
                                        tend=int(pdb_hits_pd.loc[i, 'tend']),
                                        aln_query=pdb_hits_pd.loc[i, 'qaln'],
                                        qstart=int(pdb_hits_pd.loc[i, 'qstart']),
                                        qend=int(pdb_hits_pd.loc[i, 'qend']),
                                        from_predicted_structure=False)

            if from_pdb:
                # add chain info to the original atom file
                ori_atom_file = os.path.join(self._input_pdb_dir, hit.template_name + '.atom')
                if not os.path.exists(ori_atom_file):
                    ori_atom_file = os.path.join(self._input_pdb_dir, hit.template_name.replace('.gz', ''))

                add_chain_info_to_atom_file(infile=ori_atom_file,
                                            chainid=hit.template_chain,
                                            outfile=os.path.join(template_pdb_dir, hit.template_name + '.atom'))
            else:
                atom_file = os.path.join(self._input_pdb_dir, hit.template_name + ".pdb")
                if not os.path.exists(atom_file):
                    atom_file = os.path.join(self._input_pdb_dir, hit.template_name + ".atom")

                add_chain_info_to_atom_file(infile=atom_file,
                                            chainid=hit.template_chain,
                                            outfile=os.path.join(template_pdb_dir, hit.template_name + '.atom'))

            result = _process_single_hit(query_sequence=query_sequence,
                                         monomer_hit=hit,
                                         atom_dir=template_pdb_dir,
                                         kalign_binary_path=self._kalign_binary_path,
                                         strict_error_check=self._strict_error_check)

            if result.error:
                errors.append(result.error)

            # There could be an error even if there are some results, e.g. thrown by
            # other unparsable chains in the same mmCIF file.
            if result.warning:
                warnings.append(result.warning)

            if result.features is None:
                logging.info('Skipped invalid hit %s, error: %s, warning: %s',
                             hit.template_name, result.error, result.warning)
            else:
                # Increment the hit counter, since we got features out of this hit.
                indices_map.append(num_hits)
                num_hits += 1
                for k in template_features:
                    template_features[k].append(result.features[k])
                hits_features += [result.features]

        for name in template_features:
            if num_hits > 0:
                template_features[name] = np.stack(
                    template_features[name], axis=0).astype(TEMPLATE_FEATURES[name])
            else:
                # Make sure the feature has correct dtype even if empty.
                template_features[name] = np.array([], dtype=TEMPLATE_FEATURES[name])

        return TemplateSearchResult(features=template_features, hits_features=hits_features,
                                    errors=errors, warnings=warnings)

    def get_templates_alphafold(self,
                                targetname: str,
                                query_sequence: str,
                                template_pdb_dir: str) -> TemplateSearchResult:

        """Computes the templates for given query sequence (more details above)."""
        logging.info('Searching for template for: %s', query_sequence)

        template_features = {}
        for template_feature_name in TEMPLATE_FEATURES:
            template_features[template_feature_name] = []

        num_hits = 0
        errors = []
        warnings = []
        indices_map = []
        hits_features = []

        for n in range(0, self._max_hits):

            hit = CustomizedTemplateHit(query_name=targetname,
                                        template_name=targetname + f'_{n}',
                                        template_chain=targetname[4],
                                        aligned_length=len(query_sequence),
                                        aln_temp=query_sequence,
                                        tstart=1,
                                        tend=len(query_sequence),
                                        aln_query=query_sequence,
                                        qstart=1,
                                        qend=len(query_sequence),
                                        from_predicted_structure=True)

            # add chain info to the original atom file

            add_chain_info_to_atom_file(infile=os.path.join(self._input_pdb_dir, f'ranked_{n}.pdb'),
                                        chainid=hit.template_chain,
                                        outfile=os.path.join(template_pdb_dir, hit.template_name + '.atom'))

            result = _process_single_hit(query_sequence=query_sequence,
                                         monomer_hit=hit,
                                         atom_dir=template_pdb_dir,
                                         kalign_binary_path=self._kalign_binary_path,
                                         strict_error_check=self._strict_error_check)

            if result.error:
                errors.append(result.error)

            # There could be an error even if there are some results, e.g. thrown by
            # other unparsable chains in the same mmCIF file.
            if result.warning:
                warnings.append(result.warning)

            if result.features is None:
                logging.info('Skipped invalid hit %s, error: %s, warning: %s',
                             hit.template_name, result.error, result.warning)
            else:
                # Increment the hit counter, since we got features out of this hit.
                indices_map.append(num_hits)
                num_hits += 1
                for k in template_features:
                    template_features[k].append(result.features[k])
                hits_features += [result.features]

        for name in template_features:
            if num_hits > 0:
                template_features[name] = np.stack(
                    template_features[name], axis=0).astype(TEMPLATE_FEATURES[name])
            else:
                # Make sure the feature has correct dtype even if empty.
                template_features[name] = np.array([], dtype=TEMPLATE_FEATURES[name])

        return TemplateSearchResult(features=template_features, hits_features=hits_features,
                                    errors=errors, warnings=warnings)
