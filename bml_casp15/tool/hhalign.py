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


# Internal import (7716).

class HHAlign:
    """Python wrapper of the HHsearch binary."""

    def __init__(self,
                 *,
                 binary_path: str):

        self.binary_path = binary_path

    @property
    def output_format(self) -> str:
        return 'hhr'

    @property
    def input_format(self) -> str:
        return 'fasta'

    def query(self, src_file: str, trg_file: str, outdir: str) -> str:
        """Queries the database using HHsearch using a given a3m."""

        if not os.path.exists(trg_file):
            logging.error('Could not find template file for HHalign %s', trg_file)
            raise ValueError(f'Could not find template file for HHalign {trg_file}')

        # filename, input_prefix = os.path.splitext(src_file)
        # input_path = os.path.join(outdir, 'query' + input_prefix)
        # os.system(f"cp {src_file} {input_path}")
        #
        # filename, input_prefix_t = os.path.splitext(trg_file)
        # input_path_t = os.path.join(outdir, 'template' + input_prefix_t)
        # os.system(f"cp {trg_file} {input_path_t}")

        hhr_path = os.path.join(outdir, 'output.hhr')

        cmd = [self.binary_path,
               '-i', src_file,
               '-t', trg_file,
               '-o', hhr_path
               ]

        logging.info('Launching subprocess "%s"', ' '.join(cmd))
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with utils.timing('HHalign query'):
            stdout, stderr = process.communicate()
            retcode = process.wait()

        if retcode:
            # Stderr is truncated to prevent proto size errors in Beam.
            raise RuntimeError(
                'HHSearch failed:\nstdout:\n%s\n\nstderr:\n%s\n' % (
                    stdout.decode('utf-8'), stderr[:100_000].decode('utf-8')))

        with open(hhr_path) as f:
            hhr = f.read()
        return hhr
