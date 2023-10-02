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
from multicom3.tool import utils
import pandas as pd
import pathlib


# Internal import (7716).

class Foldseek:
    """Python wrapper of the HHsearch binary."""

    def __init__(self,
                 *,
                 binary_path: str,
                 pdb_database: str,
                 max_template_date: datetime.datetime,
                 release_dates: Mapping[str, datetime.datetime],
                 other_databases: Sequence[str] = []):
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
        self.pdb_database = pdb_database
        self.other_databases = other_databases
        self._release_dates = release_dates
        self._max_template_date = max_template_date

        for database_path in self.other_databases:
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

        databases = [self.pdb_database] + self.other_databases
        for database in databases:
            database_name = pathlib.Path(database).stem
            outfile = os.path.join(outdir, f'aln.m8_{database_name}')
            if not os.path.exists(outfile):
                cmd = [self.binary_path,
                       'easy-search',
                       input_path,
                       database,
                       outfile,
                       os.path.join(outdir, 'tmp'),
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
            
            if database == self.pdb_database:
                keep_indices = []
                pdb_df = pd.read_csv(outfile, sep='\t')
                for i in range(len(pdb_df)):
                    target = pdb_df.loc[i, 'target']
                    if target.lower()[:4] in self._release_dates:
                        hit_release_date = datetime.datetime.strptime(self._release_dates[target.lower()[:4]], '%Y-%m-%d')
                        if hit_release_date < max_template_date:
                            keep_indices += [i]

                pdb_df_filtered = pdb_df.iloc[keep_indices]
                pdb_df_filtered.drop(pdb_df_filtered.filter(regex="Unnamed"), axis=1, inplace=True)
                pdb_df_filtered.reset_index(inplace=True, drop=True)
                pdb_df_filtered.to_csv(outfile,, sep='\t')

            evalue_df = evalue_df.append(pd.read_csv(outfile, sep='\t'))

        evalue_df = evalue_df.sort_values(by='evalue')
        evalue_df.reset_index(inplace=True, drop=True)
        evalue_df.to_csv(os.path.join(outdir, "evalue.m8"), sep='\t')

        if len(evalue_df) < progressive_threshold:
            # search the database using tmalign mode
            for database in self.databases:
                database_name = pathlib.Path(database).stem
                outfile = os.path.join(outdir, f'aln.m8_{database_name}.tm')
                if not os.path.exists(outfile):
                    cmd = [self.binary_path,
                           'easy-search',
                           input_path,
                           database,
                           outfile,
                           os.path.join(outdir, 'tmp'),
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

                if database == self.pdb_database:
                    keep_indices = []
                    pdb_df = pd.read_csv(outfile, sep='\t')
                    for i in range(len(pdb_df)):
                        target = pdb_df.loc[i, 'target']
                        if target.lower()[:4] in self._release_dates:
                            hit_release_date = datetime.datetime.strptime(self._release_dates[target.lower()[:4]], '%Y-%m-%d')
                            if hit_release_date < max_template_date:
                                keep_indices += [i]

                    pdb_df_filtered = pdb_df.iloc[keep_indices]
                    pdb_df_filtered.drop(pdb_df_filtered.filter(regex="Unnamed"), axis=1, inplace=True)
                    pdb_df_filtered.reset_index(inplace=True, drop=True)
                    pdb_df_filtered.to_csv(outfile,, sep='\t')

                tmscore_df = tmscore_df.append(pd.read_csv(outfile, sep='\t'))
            
            tmscore_df = tmscore_df.sort_values(by='evalue', ascending=False)
            tmscore_df.reset_index(inplace=True, drop=True)
            tmscore_df.to_csv(os.path.join(outdir, "tmscore.m8"), sep='\t')

        result_df = result_df.append(evalue_df)
        result_df = result_df.append(tmscore_df)
        result_df.reset_index(inplace=True, drop=True)
        result_df.to_csv(os.path.join(outdir, "result.m8"), sep='\t')

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
        databases = [self.pdb_database] + self.other_databases
        for database in databases:
            database_name = pathlib.Path(database).stem
            outfile = os.path.join(outdir, f'aln.m8_{database_name}.tm')
            if not os.path.exists(outfile):
                cmd = [self.binary_path,
                       'easy-search',
                       input_path,
                       database,
                       outfile,
                       os.path.join(outdir, 'tmp'),
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

            if database == self.pdb_database:
                keep_indices = []
                pdb_df = pd.read_csv(outfile, sep='\t')
                for i in range(len(pdb_df)):
                    target = pdb_df.loc[i, 'target']
                    if target.lower()[:4] in self._release_dates:
                        hit_release_date = datetime.datetime.strptime(self._release_dates[target.lower()[:4]], '%Y-%m-%d')
                        if hit_release_date < max_template_date:
                            keep_indices += [i]

                pdb_df_filtered = pdb_df.iloc[keep_indices]
                pdb_df_filtered.drop(pdb_df_filtered.filter(regex="Unnamed"), axis=1, inplace=True)
                pdb_df_filtered.reset_index(inplace=True, drop=True)
                pdb_df_filtered.to_csv(outfile,, sep='\t')

            result_df = result_df.append(pd.read_csv(outfile, sep='\t'))

        result_df = result_df.sort_values(by='evalue', ascending=False)
        result_df.reset_index(inplace=True, drop=True)
        result_df.to_csv(os.path.join(outdir, "tmscore.m8"), sep='\t')

        return result_df

    def query_only_local(self, pdb: str, outdir: str, maxseq=2000) -> str:
        """Queries the database using HHsearch using a given a3m."""
        input_path = os.path.join(outdir, 'query.pdb')
        os.system(f"cp {pdb} {input_path}")

        result_df = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])

        databases = [self.pdb_database] + self.other_databases
        for database in databases:
            database_name = pathlib.Path(database).stem
            outfile = os.path.join(outdir, f'aln.m8_{database_name}')
            if not os.path.exists(outfile):
                cmd = [self.binary_path,
                       'easy-search',
                       input_path,
                       database,
                       outfile,
                       os.path.join(outdir, 'tmp'),
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

            if database == self.pdb_database:
                keep_indices = []
                pdb_df = pd.read_csv(outfile, sep='\t')
                for i in range(len(pdb_df)):
                    target = pdb_df.loc[i, 'target']
                    if target.lower()[:4] in self._release_dates:
                        hit_release_date = datetime.datetime.strptime(self._release_dates[target.lower()[:4]], '%Y-%m-%d')
                        if hit_release_date < max_template_date:
                            keep_indices += [i]

                pdb_df_filtered = pdb_df.iloc[keep_indices]
                pdb_df_filtered.drop(pdb_df_filtered.filter(regex="Unnamed"), axis=1, inplace=True)
                pdb_df_filtered.reset_index(inplace=True, drop=True)
                pdb_df_filtered.to_csv(outfile,, sep='\t')

            result_df = result_df.append(pd.read_csv(outfile, sep='\t'))

        result_df = result_df.sort_values(by='evalue')
        result_df.reset_index(inplace=True, drop=True)
        result_df.to_csv(os.path.join(outdir, "result.m8"), sep='\t')

        empty_df = pd.DataFrame(
            columns=['query', 'target', 'qaln', 'taln', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'alnlen'])
        return {'local_alignment': result_df,
                'global_alignment': empty_df,
                'all_alignment': result_df}