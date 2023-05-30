import collections
import dataclasses
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Set, Mapping


class CustomizedInputs_Monomer:
    def __init__(self):
        self.fasta_path = None
        self.uniref_a3m = None
        self.bfd_a3m = None
        self.mgnify_sto = None
        self.uniref_bfd_a3m = None
        self.custom_msa = None
        self.uniref90_sto = None
        self.notemplate = False
        self.temp_struct_csv = None
        self.struct_atom_dir = None


class CustomizedInputs_Multimer:
    def __init__(self):
        self.fasta_path = ""
        self.monomer_a3ms = {}
        self.multimer_a3ms = {}
        self.msa_pair_file = ""
        self.template_stos = {}
        self.temp_struct_csv = ""
        self.template_hits_files = {}
        self.temp_seq_pair_file = ""
        self.monomer_model_paths = []
        self.monomer_temp_csvs = []
        self.notemplate = False
