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

        logfile = os.path.join(outpath, 'img.running')
        if not os.path.exists(logfile):
            tmpdir = os.path.join(outpath, 'tmp')
            cmd = f"python {self.binary_path} " \
                  f"{input_fasta_path} " \
                  f"-hhblitsdb={self.bfd_database_path} " \
                  f"-jackhmmerdb={self.img_database_path}  " \
                  f"-hmmsearchdb={self.metaclust_database_path}:{self.mgnify_database_path}:{self.uniref90_database_path} " \
                  f"-outdir={outpath} -tmpdir={tmpdir} -ncpu=8  &> {logfile} &"

            print(cmd)
            os.system(cmd)

        return os.path.join(outpath, targetname+".a3m")
