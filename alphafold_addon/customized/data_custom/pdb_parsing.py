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

"""Parses the mmCIF file format."""
import collections
import dataclasses
import functools
import io
from typing import Any, Mapping, Optional, Sequence, Tuple

from absl import logging
from Bio import PDB
from Bio.Data import SCOPData

# Type aliases:
ChainId = str
PdbHeader = Mapping[str, Any]
PdbStructure = PDB.Structure.Structure
SeqRes = str


@dataclasses.dataclass(frozen=True)
class Monomer:
    id: str
    num: int


# Used to map SEQRES index to a residue in the structure.
@dataclasses.dataclass(frozen=True)
class ResiduePosition:
    chain_id: str
    residue_number: int
    insertion_code: str


@dataclasses.dataclass(frozen=True)
class ResidueAtPosition:
    position: Optional[ResiduePosition]
    name: str
    is_missing: bool
    hetflag: str


@dataclasses.dataclass(frozen=True)
class PdbObject:
    """Representation of a parsed mmCIF file.

  Contains:
    file_id: A meaningful name, e.g. a pdb_id. Should be unique amongst all
      files being processed.
    header: Biopython header.
    structure: Biopython structure.
    chain_to_seqres: Dict mapping chain_id to 1 letter amino acid sequence. E.g.
      {'A': 'ABCDEFG'}
    seqres_to_structure: Dict; for each chain_id contains a mapping between
      SEQRES index and a ResidueAtPosition. e.g. {'A': {0: ResidueAtPosition,
                                                        1: ResidueAtPosition,
                                                        ...}}
    raw_string: The raw string used to construct the MmcifObject.
  """
    file_id: str
    header: PdbHeader
    structure: PdbStructure
    chain_to_seqres: Mapping[ChainId, SeqRes]
    seqres_to_structure: Mapping[ChainId, Mapping[int, ResidueAtPosition]]
    raw_string: Any


@dataclasses.dataclass(frozen=True)
class ParsingResult:
    """Returned by the parse function.

  Contains:
    mmcif_object: A MmcifObject, may be None if no chain could be successfully
      parsed.
    errors: A dict mapping (file_id, chain_id) to any exception generated.
  """
    pdb_object: Optional[PdbObject]
    errors: Mapping[Tuple[str, str], Any]


class ParseError(Exception):
    """An error indicating that an mmCIF file could not be parsed."""


@functools.lru_cache(16, typed=False)
def parse(*,
          file_id: str,
          chain_id: str,
          pdb_string: str,
          catch_all_errors: bool = True) -> ParsingResult:
    """Entry point, parses an mmcif_string.

  Args:
    file_id: A string identifier for this file. Should be unique within the
      collection of files being processed.
    mmcif_string: Contents of an mmCIF file.
    catch_all_errors: If True, all exceptions are caught and error messages are
      returned as part of the ParsingResult. If False exceptions will be allowed
      to propagate.

  Returns:
    A ParsingResult.
  """
    errors = {}
    try:
        parser = PDB.PDBParser(QUIET=True)
        handle = io.StringIO(pdb_string)
        full_structure = parser.get_structure('', handle)
        first_model_structure = _get_first_model(full_structure)

        header = parser.get_header()
        if header['resolution'] is None:
            header['resolution'] = 0.0

        chain_residue_list = []
        for residue in first_model_structure.get_residues():
            chain_residue_list.append(Monomer(id=residue.get_resname(), num=int(residue.get_full_id()[3][1])))

        if len(chain_residue_list) == 0:
            return ParsingResult(None, {(file_id, ''): 'No protein chains found in this file.'})

        valid_chains = {chain_id: chain_residue_list}

        seq_start_num = {val_chain_id: min([monomer.num for monomer in seq])
                         for val_chain_id, seq in valid_chains.items()}

        seq_to_structure_mappings = {}
        for atom in first_model_structure.get_atoms():
            hetflag = ' '
            insertion_code = ' '

            residue = atom.get_parent()
            position = ResiduePosition(chain_id=chain_id,
                                       residue_number=int(residue.get_full_id()[3][1]),
                                       insertion_code=insertion_code)

            seq_idx = int(residue.get_full_id()[3][1]) - seq_start_num[chain_id]

            current = seq_to_structure_mappings.get(chain_id, {})

            current[seq_idx] = ResidueAtPosition(position=position,
                                                 name=residue.get_resname(),
                                                 is_missing=False,
                                                 hetflag=hetflag)

            seq_to_structure_mappings[chain_id] = current

        # Add missing residue information to seq_to_structure_mappings.
        # for chain_id, seq_info in valid_chains.items():
        #     author_chain = mmcif_to_author_chain_id[chain_id]
        #     current_mapping = seq_to_structure_mappings[author_chain]
        #     for idx, monomer in enumerate(seq_info):
        #         if idx not in current_mapping:
        #             current_mapping[idx] = ResidueAtPosition(position=None,
        #                                                      name=monomer.id,
        #                                                      is_missing=True,
        #                                                      hetflag=' ')

        author_chain_to_sequence = {}
        for val_chain_id, seq_info in valid_chains.items():
            seq = []
            for monomer in seq_info:
                code = SCOPData.protein_letters_3to1.get(monomer.id, 'X')
                seq.append(code if len(code) == 1 else 'X')
            seq = ''.join(seq)
            author_chain_to_sequence[val_chain_id] = seq

        pdb_object = PdbObject(
            file_id=file_id,
            header=header,
            structure=first_model_structure,
            chain_to_seqres=author_chain_to_sequence,
            seqres_to_structure=seq_to_structure_mappings,
            raw_string=pdb_string)

        return ParsingResult(pdb_object=pdb_object, errors=errors)
    except Exception as e:  # pylint:disable=broad-except
        errors[(file_id, '')] = e
        if not catch_all_errors:
            raise
        return ParsingResult(pdb_object=None, errors=errors)


def _get_first_model(structure: PdbStructure) -> PdbStructure:
    """Returns the first model in a Biopython structure."""
    return next(structure.get_models())