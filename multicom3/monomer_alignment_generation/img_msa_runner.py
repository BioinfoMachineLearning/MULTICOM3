# Modified from Alphafold2 codes

"""Library to run HHblits from Python."""

import glob
import os
import subprocess
from typing import Any, Mapping, Optional, Sequence
from absl import logging
from multicom3.tool import utils


class IMG_Msa_runner:
    """Python wrapper of the HHblits binary."""

    def __init__(self,
                 binary_path,
                 bfd_database_path,
                 img_database_path,
                 metaclust_database_path,
                 mgnify_database_path,
                 uniref90_database_path
                 ):
        self.img_database_path = img_database_path
        self.bfd_database_path = bfd_database_path
        self.metaclust_database_path = metaclust_database_path
        self.mgnify_database_path = mgnify_database_path
        self.uniref90_database_path = uniref90_database_path
        self.binary_path = binary_path

    def query(self, input_fasta_path: str, outpath: str) -> Mapping[str, Any]:
        """Queries the database using HHblits."""

        targetname = open(input_fasta_path).readlines()[0].rstrip('\n').lstrip('>')

        if not os.path.exists(outpath + '/img.running'):
            cmd = f"python {self.binary_path} " \
                  f"{input_fasta_path} " \
                  f"-hhblitsdb={self.bfd_database_path} " \
                  f"-jackhmmerdb={self.img_database_path}  " \
                  f"-hmmsearchdb={self.metaclust_database_path}:{self.mgnify_database_path}:{self.uniref90_database_path} " \
                  f"-outdir={outpath} -tmpdir={outpath}/tmp -ncpu=8  &> {outpath}/img.running &"

            print(cmd)
            os.system(cmd)
        #
        # if not os.path.exist(f'{outpath}/{targetname}.a3m'):
        #     raise Exception(f"Deepmsa failed: {cmd}")
        #
        # os.system(f"cp {outpath}/{targetname}.a3m {output_a3m_path}")

        return f"{outpath}/{targetname}.a3m"
