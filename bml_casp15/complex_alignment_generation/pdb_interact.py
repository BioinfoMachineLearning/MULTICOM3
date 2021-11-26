import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *


class PDB_interact:

    def __init__(self, uniprot2pdb_mapping_file, dimers_list):

        self.mapping_file = uniprot2pdb_mapping_file
        self.dimers_list = dimers_list

    def load_data(self):

        self.uniprot2pdb_map = {}
        self.dimers_map = {}

        with open(self.mapping_file) as f:
            for line in f:
                # UniRef100_P13744        2E9QA
                uniprot_id, pdbcode = line.rstrip('\n').split()
                uniprot_id = uniprot_id.split('_')[1]
                if uniprot_id in self.uniprot2pdb_map:
                    self.uniprot2pdb_map[uniprot_id] += "," + pdbcode
                else:
                    self.uniprot2pdb_map[uniprot_id] = pdbcode

        with open(self.dimers_list) as f:
            for line in f:
                pdbcode, chain1, chain2 = line.rstrip('\n').split()
                if chain1 in self.dimers_map:
                    self.dimers_map[chain1] += "," + chain2
                else:
                    self.dimers_map[chain1] = chain2

                if chain2 in self.dimers_map:
                    self.dimers_map[chain2] += "," + chain1
                else:
                    self.dimers_map[chain2] = chain1

    def get_interactions(self, alignment1, alignment2):

        id_1 = []
        id_2 = []
        for id1 in alignment1.ids:
            for id2 in alignment2.ids:
                if id1 == id2:
                    id_1 += [id1]
                    id_2 += [id2]
                elif id1 in self.uniprot2pdb_map and id2 in self.uniprot2pdb_map:
                    pdbcodes1 = self.uniprot2pdb_map[id1]
                    pdbcodes2 = self.uniprot2pdb_map[id2]
                    interact = False
                    for pdbcode1 in pdbcodes1.split(','):
                        if interact:
                            break
                        if pdbcode1 not in self.dimers_map:
                            continue
                        for pdbcode2 in pdbcodes2.split(','):
                            if pdbcode2 in self.dimers_map[pdbcode1].split(','):
                                interact = True
                                break
                    if interact:
                        id_1 += [id1]
                        id_2 += [id2]
                        
        return pd.DataFrame(
            {"id_1": id_1, "id_2": id_2}
        )
