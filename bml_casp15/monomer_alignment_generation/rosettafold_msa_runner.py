# Modified from Alphafold2 codes

"""Library to run HHblits from Python."""

import glob
import os
import subprocess
from typing import Any, Mapping, Optional, Sequence
from absl import logging
from bml_casp15.tool import utils


class RosettaFold_Msa_runner:
    """Python wrapper of the HHblits binary."""

    def __init__(self,
                 *,
                 hhblits_binary_path,
                 hhfilter_binary_path,
                 uniref30_database_path,
                 bfd_database_path,
                 n_cpu=8,
                 n_iter=4,
                 n_maxmem=3,
                 maxfilt=100000000,
                 neffmax=20,
                 cov=25,
                 e_values=['1e-30', '1e-10', '1e-6', '1e-3'],
                 maxseq=1000000,
                 realign_max=100000000):

        self.hhblits_binary_path = hhblits_binary_path
        self.hhfilter_binary_path = hhfilter_binary_path
        self.uniref30_database_path = uniref30_database_path
        self.bfd_database_path = bfd_database_path

        print(f"Using database: {uniref30_database_path} and {bfd_database_path}")
        if not glob.glob(uniref30_database_path + '_*') or not glob.glob(bfd_database_path + '_*'):
            logging.error('Could not find HHBlits database %s, %s', uniref30_database_path, bfd_database_path)
            raise ValueError(f'Could not find HHBlits database {uniref30_database_path}, {bfd_database_path}')

        self.n_cpu = n_cpu
        self.n_iter = n_iter
        self.n_maxmem = n_maxmem
        self.maxfilt = maxfilt
        self.neffmax = neffmax
        self.cov = cov
        self.e_values = e_values
        self.maxseq = maxseq
        self.realign_max = realign_max

    def query(self, input_fasta_path: str, output_a3m_path: str) -> Mapping[str, Any]:
        """Queries the database using HHblits."""

        targetname = open(input_fasta_path).readlines()[0].rstrip('\n').lstrip('>')

        outpath = os.path.dirname(os.path.abspath(output_a3m_path))

        cmd_prefix = [
            self.hhblits_binary_path,
            '-i', input_fasta_path,
            '-o', '/dev/null',
            '-mact', '0.35',
            '-maxfilt', str(self.maxfilt),
            '-neffmax', str(self.neffmax),
            '-cov', str(self.cov),
            '-cpu', str(self.n_cpu),
            '-nodiff',
            '-realign_max', str(self.realign_max),
            '-maxseq', str(self.maxseq),
            '-maxmem', str(self.n_maxmem),
            '-n', str(self.n_iter),
            '-d', str(self.uniref30_database_path),
            '-d', str(self.bfd_database_path)]

        tmp_dir = outpath + '/hhblits'

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        for evalue in self.e_values:
            if not os.path.exists(f"{tmp_dir}/t000_.{evalue}.a3m"):
                cmd = cmd_prefix + ['-oa3m', f"{tmp_dir}/t000_.{evalue}.a3m",
                                    '-e', evalue,
                                    '-v', '0']
                logging.info('Launching subprocess "%s"', ' '.join(cmd))

                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                with utils.timing('HHblits query'):
                    stdout, stderr = process.communicate()
                    retcode = process.wait()

                if retcode:
                    # Logs have a 15k character limit, so log HHblits error line by line.
                    logging.error('HHblits failed. HHblits stderr begin:')
                    for error_line in stderr.decode('utf-8').splitlines():
                        if error_line.strip():
                            logging.error(error_line.strip())
                    logging.error('HHblits stderr end')
                    raise RuntimeError('HHblits failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                        stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

            loop_dict = {"75": 2000,
                         "50": 4000}

            for coverage in loop_dict:

                cmd = [self.hhfilter_binary_path,
                       '-id', '90',
                       '-cov', coverage,
                       '-i', f'{tmp_dir}/t000_.{evalue}.a3m',
                       '-o', f'{tmp_dir}/t000_.{evalue}.id90cov{coverage}.a3m']

                logging.info('Launching subprocess "%s"', ' '.join(cmd))

                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                with utils.timing('HHfilter'):
                    stdout, stderr = process.communicate()
                    retcode = process.wait()

                if retcode:
                    # Logs have a 15k character limit, so log HHblits error line by line.
                    logging.error('HHfilter failed. HHfilter stderr begin:')
                    for error_line in stderr.decode('utf-8').splitlines():
                        if error_line.strip():
                            logging.error(error_line.strip())
                    logging.error('HHfilter stderr end')
                    raise RuntimeError('HHfilter failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                        stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

                ret_msg = os.popen(f'grep -c "^>" {tmp_dir}/t000_.{evalue}.id90cov{coverage}.a3m').read()

                num = int(ret_msg.rstrip())

                if num > loop_dict[coverage]:
                    os.system(f"cp {tmp_dir}/t000_.{evalue}.id90cov{coverage}.a3m {output_a3m_path}")
                    return dict(a3m=output_a3m_path)

        os.system(f"cp {tmp_dir}/t000_.{evalue}.id90cov{coverage}.a3m {output_a3m_path}")
        return dict(a3m=output_a3m_path)
