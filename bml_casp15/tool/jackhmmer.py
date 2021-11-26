# Modified from Alphafold2 codes

"""Library to run Jackhmmer from Python."""

from concurrent import futures
import glob
import os
import subprocess
from typing import Any, Callable, Mapping, Optional, Sequence
from urllib import request
from absl import logging
from bml_casp15.tool import utils

class Jackhmmer:
  """Python wrapper of the Jackhmmer binary."""

  def __init__(self,
               *,
               binary_path: str,
               database_path: str,
               n_cpu: int = 8,
               n_iter: int = 1,
               e_value: float = 0.0001,
               z_value: Optional[int] = None,
               get_tblout: bool = False,
               filter_f1: float = 0.0005,
               filter_f2: float = 0.00005,
               filter_f3: float = 0.0000005,
               incdom_e: Optional[float] = None,
               dom_e: Optional[float] = None):
    """Initializes the Python Jackhmmer wrapper.
    Args:
      binary_path: The path to the jackhmmer executable.
      database_path: The path to the jackhmmer database (FASTA format).
      n_cpu: The number of CPUs to give Jackhmmer.
      n_iter: The number of Jackhmmer iterations.
      e_value: The E-value, see Jackhmmer docs for more details.
      z_value: The Z-value, see Jackhmmer docs for more details.
      get_tblout: Whether to save tblout string.
      filter_f1: MSV and biased composition pre-filter, set to >1.0 to turn off.
      filter_f2: Viterbi pre-filter, set to >1.0 to turn off.
      filter_f3: Forward pre-filter, set to >1.0 to turn off.
      incdom_e: Domain e-value criteria for inclusion of domains in MSA/next
        round.
      dom_e: Domain e-value criteria for inclusion in tblout.
      num_streamed_chunks: Number of database chunks to stream over.
      streaming_callback: Callback function run after each chunk iteration with
        the iteration number as argument.
    """
    self.binary_path = binary_path
    self.database_path = database_path

    print(f"Using database: {self.database_path}")
    if not os.path.exists(self.database_path):
      logging.error('Could not find Jackhmmer database %s', database_path)
      raise ValueError(f'Could not find Jackhmmer database {database_path}')

    self.n_cpu = n_cpu
    self.n_iter = n_iter
    self.e_value = e_value
    self.z_value = z_value
    self.filter_f1 = filter_f1
    self.filter_f2 = filter_f2
    self.filter_f3 = filter_f3
    self.incdom_e = incdom_e
    self.dom_e = dom_e
    self.get_tblout = get_tblout

  def query(self, input_fasta_path: str, output_sto_path: str)-> Mapping[str, Any]:
    """Queries the database chunk using Jackhmmer."""

    targetname = open(input_fasta_path).readlines()[0].rstrip('\n').lstrip('>')

    # The F1/F2/F3 are the expected proportion to pass each of the filtering
    # stages (which get progressively more expensive), reducing these
    # speeds up the pipeline at the expensive of sensitivity.  They are
    # currently set very low to make querying Mgnify run in a reasonable
    # amount of time.
    cmd_flags = [
          # Don't pollute stdout with Jackhmmer output.
          '-o', '/dev/null',
          '-A', output_sto_path,
          '--noali',
          '--F1', str(self.filter_f1),
          '--F2', str(self.filter_f2),
          '--F3', str(self.filter_f3),
          '--incE', str(self.e_value),
          # Report only sequences with E-values <= x in per-sequence output.
          '-E', str(self.e_value),
          '--cpu', str(self.n_cpu),
          '-N', str(self.n_iter)
    ]
    if self.get_tblout:
        tblout_path = output_sto_path.replace(".sto", "tblout.txt")
        cmd_flags.extend(['--tblout', tblout_path])

    if self.z_value:
        cmd_flags.extend(['-Z', str(self.z_value)])

    if self.dom_e is not None:
        cmd_flags.extend(['--domE', str(self.dom_e)])

    if self.incdom_e is not None:
        cmd_flags.extend(['--incdomE', str(self.incdom_e)])

    cmd = [self.binary_path] + cmd_flags + [input_fasta_path, self.database_path]

    logging.info('Launching subprocess "%s"', ' '.join(cmd))
    print(cmd)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with utils.timing(f'Jackhmmer ({os.path.basename(self.database_path)}) query'):
        _, stderr = process.communicate()
        retcode = process.wait()

    if retcode:
        raise RuntimeError('Jackhmmer failed\nstderr:\n%s\n' % stderr.decode('utf-8'))

    # Get e-values for each target name
    tbl = ''
    if self.get_tblout:
        with open(tblout_path) as f:
            tbl = f.read()

    raw_output = dict(
        sto=output_sto_path,
        tbl=tbl,
        stderr=stderr,
        n_iter=self.n_iter,
        e_value=self.e_value)

    return raw_output
