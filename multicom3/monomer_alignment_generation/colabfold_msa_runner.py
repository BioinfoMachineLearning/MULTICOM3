# Modified from Alphafold2 codes

"""Library to run HHblits from Python."""

import glob
import os
import subprocess
from typing import Any, Mapping, Optional, Sequence
from absl import logging
from multicom3.tool import utils


class ColabFold_Msa_runner:
    """Python wrapper of the HHblits binary."""

    def __init__(self,
                 *,
                 colabfold_search_binary_path,
                 colabfold_split_msas_binary_path,
                 mmseq_binary_path,
                 colabfold_databases):

        self.colabfold_search_binary_path = colabfold_search_binary_path
        self.colabfold_split_msas_binary_path = colabfold_split_msas_binary_path
        self.mmseq_binary_path = mmseq_binary_path
        self.colabfold_databases = colabfold_databases

        print(f"Using database: {colabfold_databases}")

        if not os.path.exists(colabfold_databases):
            logging.error('Could not find colabfold database %s', colabfold_databases)
            raise ValueError('Could not find colabfold database %s', colabfold_databases)


    def query(self, input_fasta_path: str, output_a3m_path: str) -> Mapping[str, Any]:
        """Queries the database using Colabfold."""


        targetname = open(input_fasta_path).readlines()[0].rstrip('\n').lstrip('>')

        outpath = os.path.dirname(os.path.abspath(output_a3m_path))

        tmp_dir = outpath + '/colabfold'

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        cmd = [
            'python',
            self.colabfold_search_binary_path,
            input_fasta_path,
            self.colabfold_databases,
            tmp_dir + '/search_result'
        ]

        logging.info('Colabfold subprocess "%s"', ' '.join(cmd))

        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        with utils.timing('Colabfold query'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Logs have a 15k character limit, so log HHblits error line by line.
            logging.error('HHblits failed. Colabfold stderr begin:')
            for error_line in stderr.decode('utf-8').splitlines():
                if error_line.strip():
                    logging.error(error_line.strip())
            logging.error('Colabfold stderr end')
            raise RuntimeError('Colabfold failed\nstdout:\n%s\n\nstderr:\n%s\n' % (
                stdout.decode('utf-8'), stderr[:500_000].decode('utf-8')))

        os.system(f"python {self.colabfold_split_msas_binary_path} {tmp_dir}/search_result {tmp_dir}/msas")

        if not os.path.exists(f"{tmp_dir}/msas/{targetname}.a3m"):
            raise RuntimeError(f"Cannot find the generated a3m file for {targetname}")

        # need to remove the first line
        open(output_a3m_path, 'w').writelines(open(f"{tmp_dir}/msas/{targetname}.a3m").readlines()[1:])
        # os.system(f"cp {tmp_dir}/msas/{targetname}.a3m {output_a3m_path}")
        os.system(f"rm -rf {tmp_dir}")
        return dict(a3m=output_a3m_path)
