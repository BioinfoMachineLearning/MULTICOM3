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

"""Library to run HHsearch from Python."""

import glob
import os
import subprocess
from typing import Sequence

from absl import logging
from bml_casp15.tool import utils
import pandas as pd
import pathlib


# Internal import (7716).

class Foldseek:
    """Python wrapper of the HHsearch binary."""

    def __init__(self,
                 *,
                 binary_path: str,
                 databases: Sequence[str],
                 maxseq: int = 1_000_000):
        """Initializes the Python HHsearch wrapper.

    Args:
      binary_path: The path to the HHsearch executable.
      databases: A sequence of HHsearch database paths. This should be the
        common prefix for the database files (i.e. up to but not including
        _hhm.ffindex etc.)
      maxseq: The maximum number of rows in an input alignment. Note that this
        parameter is only supported in HHBlits version 3.1 and higher.

    Raises:
      RuntimeError: If HHsearch binary not found within the path.
    """
        self.binary_path = binary_path
        self.databases = databases
        self.maxseq = maxseq

        for database_path in self.databases:
            if not glob.glob(database_path + '_*'):
                logging.error('Could not find HHsearch database %s', database_path)
                raise ValueError(f'Could not find HHsearch database {database_path}')

    def query(self, pdb: str, outdir: str, progressive_threshold = 20, progressive=False) -> str:
        """Queries the database using HHsearch using a given a3m."""
        input_path = os.path.join(outdir, 'query.pdb')
        os.system(f"cp {pdb} {input_path}")

        csvs = []
        result_df = pd.DataFrame(columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])
        for database in self.databases:
            database_name = pathlib.Path(database).stem
            cmd = [self.binary_path,
                   'easy-search',
                   input_path,
                   database,
                   f'{outdir}/aln.m8_{database_name}',
                   outdir + '/tmp',
                   '--format-output', 'query,target,qaln,taln,qstart,qend,tstart,tend,evalue,alnlen',
                   '--format-mode', '4']
            logging.info('Launching subprocess "%s"', ' '.join(cmd))
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with utils.timing('Foldseek query'):
                stdout, stderr = process.communicate()
                retcode = process.wait()
            if retcode:
                # Stderr is truncated to prevent proto size errors in Beam.
                raise RuntimeError(
                    'Foldseek failed:\nstdout:\n%s\n\nstderr:\n%s\n' % (
                        stdout.decode('utf-8'), stderr[:100_000].decode('utf-8')))
            csvs += [f'{outdir}/aln.m8_{database_name}']


        # search the database using tmalign mode
        if len(result_df) < progressive_threshold and progressive:
            for database in self.databases:
                database_name = pathlib.Path(database).stem
                cmd = [self.binary_path,
                       'easy-search',
                       input_path,
                       database,
                       f'{outdir}/aln.m8_{database_name}.tm',
                       outdir + '/tmp',
                       '--format-output', 'query,target,qaln,taln,qstart,qend,tstart,tend,evalue,alnlen',
                       '--format-mode', '4',
                       '--alignment-type', '1',
                       '--tmscore-threshold', '0.3']
                logging.info('Launching subprocess "%s"', ' '.join(cmd))
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                with utils.timing('Foldseek query'):
                    stdout, stderr = process.communicate()
                    retcode = process.wait()
                if retcode:
                    # Stderr is truncated to prevent proto size errors in Beam.
                    raise RuntimeError(
                        'Foldseek failed:\nstdout:\n%s\n\nstderr:\n%s\n' % (
                            stdout.decode('utf-8'), stderr[:100_000].decode('utf-8')))
                csvs += [f'{outdir}/aln.m8_{database_name}.tm']

        for csv in csvs:
            result_df = result_df.append(pd.read_csv(csv, sep='\t'))

        result_df = result_df.sort_values(by='evalue')
        result_df.to_csv(f"{outdir}/result.m8", sep='\t')
        return f"{outdir}/result.m8"
