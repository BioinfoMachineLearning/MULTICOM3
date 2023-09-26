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

"""Docker launch script for Alphafold docker image."""

import os
import pathlib
import signal
from typing import Tuple

from absl import app
from absl import flags
from absl import logging
import docker
from docker import types

flags.DEFINE_string('mode', None, 'Monomer, heteromer or homomer')
flags.DEFINE_string('fasta_path', None, 'Path to fasta')
flags.DEFINE_string('output_dir', '/tmp/multicom3',
                    'Path to a directory that will store the results.')     
flags.DEFINE_string('af_db_dir', None,
                    'Path to directory with supporting data: AlphaFold parameters and genetic '
                    'and template databases. Set to the target of download_all_databases.sh.')
flags.DEFINE_string('multicom3_db_dir', None,
                    'Path to directory with additional data for MULTICOM3')          
flags.DEFINE_string('docker_image_name', 'multicom3', 'Name of the MULTICOM3 Docker image.')
flags.DEFINE_boolean('run_img', False, 'Whether to use IMG alignment to generate models')    
flags.DEFINE_string('docker_user', f'{os.geteuid()}:{os.getegid()}',
                    'UID:GID with which to run the Docker container. The output directories '
                    'will be owned by this user:group. By default, this is the current user. '
                    'Valid options are: uid or uid:gid, non-numeric values are not recognised '
                    'by Docker unless that user has been created within the container.')
flags.DEFINE_string(
    'gpu_devices', 'all',
    'Comma separated list of devices to pass to NVIDIA_VISIBLE_DEVICES.')

FLAGS = flags.FLAGS

_ROOT_MOUNT_DIRECTORY = '/mnt/'


def _create_mount(mount_name: str, path: str) -> Tuple[types.Mount, str]:
  """Create a mount point for each file and directory used by the model."""
  path = pathlib.Path(path).absolute()
  target_path = pathlib.Path(_ROOT_MOUNT_DIRECTORY, mount_name)

  if path.is_dir():
    source_path = path
    mounted_path = target_path
  else:
    source_path = path.parent
    mounted_path = pathlib.Path(target_path, path.name)
  if not source_path.exists():
    raise ValueError(f'Failed to find source directory "{source_path}" to '
                     'mount in Docker container.')
  logging.info('Mounting %s -> %s', source_path, target_path)
  mount = types.Mount(target=str(target_path), source=str(source_path),
                      type='bind', read_only=True)
  return mount, str(mounted_path)


def main(argv):
  if len(argv) > 1:
    raise app.UsageError('Too many command-line arguments.')

  mounts = []
  command_args = [f'-m {FLAGS.mode}']
  command_args.append(f'-i {FLAGS.run_img}')

  mount, target_path = _create_mount(f'fasta_path', FLAGS.fasta_path)
  mounts.append(mount)
  command_args.append(f'-f {target_path}')

  mount, target_path = _create_mount("multicom3_db", FLAGS.multicom3_db_dir)
  mounts.append(mount)
  command_args.append(f'-d {target_path}')

  mount, target_path = _create_mount("af_db", FLAGS.af_db_dir)
  mounts.append(mount)
  command_args.append(f'-a {target_path}')

  output_target_path = os.path.join(_ROOT_MOUNT_DIRECTORY, 'output')
  mounts.append(types.Mount(output_target_path, FLAGS.output_dir, type='bind'))

  command_args.append(f'-o {output_target_path}')

  client = docker.from_env()
  device_requests = [docker.types.DeviceRequest(driver='nvidia', capabilities=[['gpu']])]
  print(command_args)
  container = client.containers.run(
      image=FLAGS.docker_image_name,
      command=command_args,
      device_requests=device_requests,
      remove=True,
      detach=True,
      mounts=mounts,
      user=FLAGS.docker_user,
      environment={
          'NVIDIA_VISIBLE_DEVICES': FLAGS.gpu_devices,
          # The following flags allow us to make predictions on proteins that
          # would typically be too long to fit into GPU memory.
          'TF_FORCE_UNIFIED_MEMORY': '1',
          'XLA_PYTHON_CLIENT_MEM_FRACTION': '4.0',
      })

  # Add signal handler to ensure CTRL+C also stops the running container.
  signal.signal(signal.SIGINT,
                lambda unused_sig, unused_frame: container.kill())

  for line in container.logs(stream=True):
    logging.info(line.strip().decode('utf-8'))


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'mode',
      'fasta_path',
      'af_db_dir',
      'multicom3_db_dir'
  ])
  app.run(main)
