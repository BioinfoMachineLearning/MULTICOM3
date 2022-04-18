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
                 databases: Sequence[str]):
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

        for database_path in self.databases:
            if not glob.glob(database_path + '_*'):
                logging.error('Could not find HHsearch database %s', database_path)
                raise ValueError(f'Could not find HHsearch database {database_path}')

    def query(self, pdb: str, outdir: str, progressive_threshold=1, tmscore_threshold=0.3, maxseq=2000) -> str:
        """Queries the database using HHsearch using a given a3m."""
        input_path = os.path.join(outdir, 'query.pdb')
        os.system(f"cp {pdb} {input_path}")

        result_df = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])
        evalue_df = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])
        tmscore_df = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])

        for database in self.databases:
            database_name = pathlib.Path(database).stem
            if not os.path.exists(f'{outdir}/aln.m8_{database_name}'):
                cmd = [self.binary_path,
                       'easy-search',
                       input_path,
                       database,
                       f'{outdir}/aln.m8_{database_name}',
                       outdir + '/tmp',
                       '--format-output', 'query,target,qaln,taln,qstart,qend,tstart,tend,evalue,alnlen',
                       '--format-mode', '4',
                       '--max-seqs', str(maxseq),
                       '-e', '0.001',
                       '-c', '0.5',
                       '--cov-mode', '2']
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
            evalue_df = evalue_df.append(pd.read_csv(f'{outdir}/aln.m8_{database_name}', sep='\t'))

        evalue_df = evalue_df.sort_values(by='evalue')
        evalue_df.reset_index(inplace=True, drop=True)
        evalue_df.to_csv(f"{outdir}/evalue.m8", sep='\t')

        if len(evalue_df) < progressive_threshold:
            # search the database using tmalign mode
            for database in self.databases:
                database_name = pathlib.Path(database).stem
                if not os.path.exists(f'{outdir}/aln.m8_{database_name}.tm'):
                    cmd = [self.binary_path,
                           'easy-search',
                           input_path,
                           database,
                           f'{outdir}/aln.m8_{database_name}.tm',
                           outdir + '/tmp',
                           '--format-output', 'query,target,qaln,taln,qstart,qend,tstart,tend,evalue,alnlen',
                           '--format-mode', '4',
                           '--alignment-type', '1',
                           '--tmscore-threshold', str(tmscore_threshold),
                           '--max-seqs', str(maxseq),
                           '-c', '0.5',
                           '--cov-mode', '2']
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
                tmscore_df = tmscore_df.append(pd.read_csv(f'{outdir}/aln.m8_{database_name}.tm', sep='\t'))

            tmscore_df = tmscore_df.sort_values(by='evalue', ascending=False)
            tmscore_df.reset_index(inplace=True, drop=True)
            tmscore_df.to_csv(f"{outdir}/tmscore.m8", sep='\t')

        result_df = result_df.append(evalue_df)
        result_df = result_df.append(tmscore_df)
        result_df.reset_index(inplace=True, drop=True)
        result_df.to_csv(f"{outdir}/result.m8", sep='\t')

        return {'local_alignment': evalue_df,
                'global_alignment': tmscore_df,
                'all_alignment': result_df}

    def query_with_tmalign(self, pdb: str, outdir: str, tmscore_threshold=0.3, maxseq=50000) -> str:
        """Queries the database using HHsearch using a given a3m."""
        input_path = os.path.join(outdir, 'query.pdb')
        os.system(f"cp {pdb} {input_path}")

        result_df = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])

        # search the database using tmalign mode
        for database in self.databases:
            database_name = pathlib.Path(database).stem
            if not os.path.exists(f'{outdir}/aln.m8_{database_name}.tm'):
                cmd = [self.binary_path,
                       'easy-search',
                       input_path,
                       database,
                       f'{outdir}/aln.m8_{database_name}.tm',
                       outdir + '/tmp',
                       '--format-output', 'query,target,qaln,taln,qstart,qend,tstart,tend,evalue,alnlen',
                       '--format-mode', '4',
                       '--alignment-type', '1',
                       '--tmscore-threshold', str(tmscore_threshold),
                       '--max-seqs', str(maxseq),
                       '-c', '0.5',
                       '--cov-mode', '2']
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
            result_df = result_df.append(pd.read_csv(f'{outdir}/aln.m8_{database_name}.tm', sep='\t'))

        result_df = result_df.sort_values(by='evalue', ascending=False)
        result_df.reset_index(inplace=True, drop=True)
        result_df.to_csv(f"{outdir}/tmscore.m8", sep='\t')

        return result_df

    def query_only_local(self, pdb: str, outdir: str, maxseq=2000) -> str:
        """Queries the database using HHsearch using a given a3m."""
        input_path = os.path.join(outdir, 'query.pdb')
        os.system(f"cp {pdb} {input_path}")

        result_df = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])

        for database in self.databases:
            database_name = pathlib.Path(database).stem
            if not os.path.exists(f'{outdir}/aln.m8_{database_name}'):
                cmd = [self.binary_path,
                       'easy-search',
                       input_path,
                       database,
                       f'{outdir}/aln.m8_{database_name}',
                       outdir + '/tmp',
                       '--format-output', 'query,target,qaln,taln,qstart,qend,tstart,tend,evalue,alnlen',
                       '--format-mode', '4',
                       '--max-seqs', str(maxseq),
                       '-e', '0.001',
                       '-c', '0.5',
                       '--cov-mode', '2']
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
            result_df = result_df.append(pd.read_csv(f'{outdir}/aln.m8_{database_name}', sep='\t'))

        result_df = result_df.sort_values(by='evalue')
        result_df.reset_index(inplace=True, drop=True)
        result_df.to_csv(f"{outdir}/result.m8", sep='\t')

        empty_df = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])
        return {'local_alignment': result_df,
                'global_alignment': empty_df,
                'all_alignment': result_df}