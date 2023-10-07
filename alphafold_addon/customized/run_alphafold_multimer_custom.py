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

"""Full AlphaFold protein structure prediction script."""
import enum
import json
import os
import pathlib
import pickle
import random
import shutil
import sys
import time
from typing import Any, Dict, Union

from absl import app
from absl import flags
from absl import logging
from alphafold.common import protein
from alphafold.data_custom import custom_params
from alphafold.common import residue_constants
from alphafold.data_custom import pipeline
from alphafold.data_custom import pipeline_multimer
from alphafold.data_custom import pipeline_multimer_custom
from alphafold.data_custom import templates
from alphafold.data_custom import templates_custom
from alphafold.data_custom.tools import hhsearch
from alphafold.data_custom.tools import hmmsearch
from alphafold.model import config
from alphafold.model import data
from alphafold.model import model
from alphafold.relax import relax
import numpy as np
import jax.numpy as jnp
# Internal import (7716).

@enum.unique
class ModelsToRelax(enum.Enum):
  ALL = 0
  BEST = 1
  NONE = 2

logging.set_verbosity(logging.INFO)

flags.DEFINE_string('fasta_path', None, 'Path to fasta file.')
flags.DEFINE_list('monomer_a3ms', None, 'List of paths to monomer(unpaired) a3m files')
flags.DEFINE_list('multimer_a3ms', None, 'List of paths to multimer(paired) a3m files')
flags.DEFINE_string('msa_pair_file', None, 'Path of pair file')
flags.DEFINE_list('template_stos', None, 'Path of stos for template search')

flags.DEFINE_string('temp_struct_csv', None, 'Path of complex template csv file')

flags.DEFINE_list('template_hits_files', None, 'Path of complex template csv file')
flags.DEFINE_string('temp_seq_pair_file', None, 'Path of complex template csv file')

flags.DEFINE_list('monomer_temp_csvs', None, 'Paths of monomer template csv files')
flags.DEFINE_list('monomer_model_paths', None, 'Path of alphafold models for monomers')
flags.DEFINE_integer('monomer_model_count', 5, 'Number of templates used for predicted monomer structures')

flags.DEFINE_string('struct_atom_dir', None, 'Path of template file directory')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_string('env_dir', None, 'Alphafold virtual environment directory')
flags.DEFINE_string('database_dir', None, 'Alphafold databse directory')
flags.DEFINE_string('max_template_date', '9999-06-01', 'Maximum template date')
flags.DEFINE_boolean('notemplate', False, 'Whether to use template')
flags.DEFINE_integer('num_multimer_predictions_per_model', 2, 'How many '
                                                              'predictions (each with a different random seed) will be '
                                                              'generated per model. E.g. if this is 2 and there are 5 '
                                                              'models then there will be 10 predictions per input. '
                                                              'Note: this FLAG only applies if model_preset=multimer')

flags.DEFINE_integer('multimer_num_ensemble', 1, 'Number of ensemble for generating multimer models')
flags.DEFINE_integer('multimer_num_recycle', 20, 'Number of recycles for generating multimer models')

flags.DEFINE_enum('model_preset', 'multimer',
                  ['multimer'],
                  'Choose preset model configuration - the monomer model, '
                  'the monomer model with extra ensembling, monomer model with '
                  'pTM head, or multimer model')
flags.DEFINE_boolean('benchmark', False, 'Run multiple JAX model evaluations '
                     'to obtain a timing that excludes the compilation time, '
                     'which should be more indicative of the time required for '
                     'inferencing many proteins.')
flags.DEFINE_enum_class('models_to_relax', ModelsToRelax.ALL, ModelsToRelax,
                        'The models to run the final relaxation step on. '
                        'If `all`, all models are relaxed, which may be time '
                        'consuming. If `best`, only the most confident model '
                        'is relaxed. If `none`, relaxation is not run. Turning '
                        'off relaxation might result in predictions with '
                        'distracting stereochemical violations but might help '
                        'in case you are having issues with the relaxation '
                        'stage.')
flags.DEFINE_boolean('use_gpu_relax', None, 'Whether to relax on GPU. '
                     'Relax on GPU can be much faster than CPU, so it is '
                     'recommended to enable if possible. GPUs must be available'
                     ' if this setting is enabled.')


FLAGS = flags.FLAGS

MAX_TEMPLATE_HITS = 20
RELAX_MAX_ITERATIONS = 0
RELAX_ENERGY_TOLERANCE = 2.39
RELAX_STIFFNESS = 10.0
RELAX_EXCLUDE_RESIDUES = []
RELAX_MAX_OUTER_ITERATIONS = 3

def parse_fasta(fasta_file):
    fasta_string = open(fasta_file).read()
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith('>'):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append('')
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line
    return sequences, descriptions


def _reorder_chains(pdbstring):
    new_pdbstring = []
    first_chain_id = None
    for line in pdbstring.split('\n'):
        if line.startswith('ATOM') or line.startswith('TER'):
            chain_id = line[21]
            if first_chain_id is None:
                first_chain_id = chain_id
            if first_chain_id != "A":
                new_pdbstring += [line[:21] + protein.PDB_CHAIN_IDS[protein.PDB_CHAIN_IDS.find(chain_id)-1] + line[22:]]
            else:
                new_pdbstring += [line]
        else:
            new_pdbstring += [line]
    return '\n'.join(new_pdbstring)

def _check_flag(flag_name: str,
                other_flag_name: str,
                should_be_set: bool):
    if should_be_set != bool(FLAGS[flag_name].value):
        verb = 'be' if should_be_set else 'not be'
        raise ValueError(f'{flag_name} must {verb} set when running with '
                         f'"--{other_flag_name}={FLAGS[other_flag_name].value}".')

def _jnp_to_np(output: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively changes jax arrays to numpy arrays."""
    for k, v in output.items():
        if isinstance(v, dict):
            output[k] = _jnp_to_np(v)
        elif isinstance(v, jnp.ndarray):
            output[k] = np.array(v)
    return output


def predict_structure(
        fasta_path,
        fasta_name,
        output_dir,
        custom_inputs,
        data_pipeline,
        model_runners,
        amber_relaxer,
        benchmark,
        random_seed,
        models_to_relax):
    """Predicts structure using AlphaFold for the given sequence."""

    logging.info('Predicting %s', fasta_name)
    timings = {}
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    msa_output_dir = os.path.join(output_dir, 'msas')
    if not os.path.exists(msa_output_dir):
        os.makedirs(msa_output_dir)

    template_output_dir = os.path.join(output_dir, 'templates')
    if not os.path.exists(template_output_dir):
        os.makedirs(template_output_dir)

    # Get features.
    t_0 = time.time()
    feature_dict = data_pipeline.process(
        input_fasta_path=fasta_path,
        msa_output_dir=msa_output_dir,
        template_output_dir=template_output_dir,
        custom_inputs=custom_inputs)
    timings['features'] = time.time() - t_0

    # Write out features as a pickled dictionary.
    features_output_path = os.path.join(output_dir, 'features.pkl')
    with open(features_output_path, 'wb') as f:
        pickle.dump(feature_dict, f, protocol=4)

    unrelaxed_pdbs = {}
    unrelaxed_proteins = {}
    relaxed_pdbs = {}
    relax_metrics = {}
    ranking_confidences = {}

    # Run the models.
    num_models = len(model_runners)
    for model_index, (model_name, model_runner) in enumerate(model_runners.items()):
        logging.info('Running model %s on %s', model_name, fasta_name)
        t_0 = time.time()
        model_random_seed = model_index + random_seed * num_models
        processed_feature_dict = model_runner.process_features(feature_dict, random_seed=model_random_seed)
        timings[f'process_features_{model_name}'] = time.time() - t_0

        unrelaxed_pdb_path = os.path.join(output_dir, f'unrelaxed_{model_name}.pdb')
        result_output_path = os.path.join(output_dir, f'result_{model_name}.pkl')

        if not os.path.exists(unrelaxed_pdb_path) or not os.path.exists(result_output_path):
            t_0 = time.time()
            prediction_result = model_runner.predict(processed_feature_dict, random_seed=model_random_seed)
            t_diff = time.time() - t_0
            timings[f'predict_and_compile_{model_name}'] = t_diff
            logging.info('Total JAX model %s on %s predict time (includes compilation time, see --benchmark): %.1fs', model_name, fasta_name, t_diff)

            if benchmark:
                t_0 = time.time()
                model_runner.predict(processed_feature_dict, random_seed=model_random_seed)
                t_diff = time.time() - t_0
                timings[f'predict_benchmark_{model_name}'] = t_diff
                logging.info('Total JAX model %s on %s predict time (excludes compilation time): %.1fs', model_name, fasta_name, t_diff)

            plddt = prediction_result['plddt']
            ranking_confidences[model_name] = prediction_result['ranking_confidence']

            # Remove jax dependency from results.
            np_prediction_result = _jnp_to_np(dict(prediction_result))

            # Save the model outputs.
            result_output_path = os.path.join(output_dir, f'result_{model_name}.pkl')
            with open(result_output_path, 'wb') as f:
                pickle.dump(np_prediction_result, f, protocol=4)

            # Add the predicted LDDT in the b-factor column.
            # Note that higher predicted LDDT value means higher model confidence.
            plddt_b_factors = np.repeat(
                plddt[:, None], residue_constants.atom_type_num, axis=-1)
            unrelaxed_protein = protein.from_prediction(
                features=processed_feature_dict,
                result=prediction_result,
                b_factors=plddt_b_factors,
                remove_leading_feature_dimension=not model_runner.multimer_mode)

            unrelaxed_proteins[model_name] = unrelaxed_protein
            unrelaxed_pdbs[model_name] = protein.to_pdb(unrelaxed_protein)
            unrelaxed_pdb_path = os.path.join(output_dir, f'unrelaxed_{model_name}.pdb')
            with open(unrelaxed_pdb_path, 'w') as f:
                f.write(unrelaxed_pdbs[model_name])
        else:
            logging.info('%s has been generated!', model_name)

            with open(result_output_path, 'rb') as f:
                prediction_result = pickle.load(f)

            plddt = prediction_result['plddt']
            ranking_confidences[model_name] = prediction_result['ranking_confidence']

            # Add the predicted LDDT in the b-factor column.
            # Note that higher predicted LDDT value means higher model confidence.
            plddt_b_factors = np.repeat(plddt[:, None], residue_constants.atom_type_num, axis=-1)

            unrelaxed_protein = protein.from_prediction(
                features=processed_feature_dict,
                result=prediction_result,
                b_factors=plddt_b_factors,
                remove_leading_feature_dimension=not model_runner.multimer_mode)

            unrelaxed_proteins[model_name] = unrelaxed_protein
            unrelaxed_pdbs[model_name] = ''.join(open(unrelaxed_pdb_path).readlines())
            # unrelaxed_pdbs[model_name] = protein.to_pdb(unrelaxed_protein)
            # with open(unrelaxed_pdb_path, 'w') as f:
            #     f.write(unrelaxed_pdbs[model_name])

    # Rank by model confidence.
    ranked_order = [model_name for model_name, confidence in sorted(ranking_confidences.items(), key=lambda x: x[1], reverse=True)]

    # Relax predictions.
    if models_to_relax == ModelsToRelax.BEST:
        to_relax = [ranked_order[0]]
    elif models_to_relax == ModelsToRelax.ALL:
        to_relax = ranked_order
    elif models_to_relax == ModelsToRelax.NONE:
        to_relax = []

    for model_name in to_relax:
        relaxed_output_path = os.path.join(output_dir, f'relaxed_{model_name}.pdb')
        if os.path.exists(relaxed_output_path):
            relaxed_pdbs[model_name] = ''.join(open(relaxed_output_path).readlines())
            print('%s has been generated!', relaxed_output_path)
        else:
            t_0 = time.time()
            relaxed_pdb_str, _, violations = amber_relaxer.process(prot=unrelaxed_proteins[model_name])
            relax_metrics[model_name] = {
                'remaining_violations': violations,
                'remaining_violations_count': sum(violations)
            }
            timings[f'relax_{model_name}'] = time.time() - t_0

            relaxed_pdbs[model_name] = relaxed_pdb_str

            # Save the relaxed PDB.
            with open(relaxed_output_path, 'w') as f:
                f.write(relaxed_pdb_str)

    # Write out relaxed PDBs in rank order.
    for idx, model_name in enumerate(ranked_order):
        ranked_output_path = os.path.join(output_dir, f'ranked_{idx}.pdb')
        with open(ranked_output_path, 'w') as f:
            if model_name in relaxed_pdbs:
                f.write(relaxed_pdbs[model_name])
            else:
                f.write(_reorder_chains(unrelaxed_pdbs[model_name]))

    ranking_output_path = os.path.join(output_dir, 'ranking_debug.json')
    with open(ranking_output_path, 'w') as f:
        label = 'iptm+ptm' if 'iptm' in prediction_result else 'plddts'
        f.write(json.dumps({label: ranking_confidences, 'order': ranked_order}, indent=4))

    logging.info('Final timings for %s: %s', fasta_name, timings)

    timings_output_path = os.path.join(output_dir, 'timings.json')
    with open(timings_output_path, 'w') as f:
        f.write(json.dumps(timings, indent=4))
    if models_to_relax != ModelsToRelax.NONE:
        relax_metrics_path = os.path.join(output_dir, 'relax_metrics.json')
        with open(relax_metrics_path, 'w') as f:
            f.write(json.dumps(relax_metrics, indent=4))


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    custom_inputs = custom_params.CustomizedInputs_Multimer()
    custom_inputs.notemplate = FLAGS.notemplate

    template_searcher = None
    template_featurizer = None
    monomer_template_featurizer = None
    template_max_hits = MAX_TEMPLATE_HITS

    if FLAGS.monomer_model_paths is not None:
        custom_inputs.monomer_model_paths += FLAGS.monomer_model_paths
        monomer_template_featurizer = templates_custom.CustomizedComplexMonomerHitFeaturizer(
            kalign_binary_path=os.path.join(FLAGS.env_dir, 'kalign'),
            monomer_model_paths=FLAGS.monomer_model_paths,
            template_count=FLAGS.monomer_model_count)
        template_max_hits -= FLAGS.monomer_model_count

    if FLAGS.template_stos is not None:

        template_searcher = hmmsearch.Hmmsearch(
            binary_path=os.path.join(FLAGS.env_dir, 'hmmsearch'),
            hmmbuild_binary_path=os.path.join(FLAGS.env_dir, 'hmmbuild'),
            database_path=os.path.join(FLAGS.database_dir, 'pdb_seqres/pdb_seqres.txt'))

        template_featurizer = templates.HmmsearchHitFeaturizer(
            mmcif_dir=os.path.join(FLAGS.database_dir, 'pdb_mmcif/mmcif_files'),
            max_template_date=FLAGS.max_template_date,
            max_hits=template_max_hits,
            kalign_binary_path=os.path.join(FLAGS.env_dir, 'kalign'),
            release_dates_path=None,
            obsolete_pdbs_path=os.path.join(FLAGS.database_dir, 'pdb_mmcif/obsolete.dat'))

    if FLAGS.temp_seq_pair_file is not None:

        template_searcher = hhsearch.HHSearch(
            binary_path=os.path.join(FLAGS.env_dir, 'hhsearch'),
            databases=[os.path.join(FLAGS.database_dir, 'pdb70/pdb70')])

        template_featurizer = templates.HhsearchHitFeaturizer(
            mmcif_dir=os.path.join(FLAGS.database_dir, 'pdb_mmcif/mmcif_files'),
            max_template_date=FLAGS.max_template_date,
            max_hits=template_max_hits,
            kalign_binary_path=os.path.join(FLAGS.env_dir, 'kalign'),
            release_dates_path=None,
            obsolete_pdbs_path=os.path.join(FLAGS.database_dir, 'pdb_mmcif/obsolete.dat'))

        custom_inputs.temp_seq_pair_file = FLAGS.temp_seq_pair_file

    if FLAGS.temp_struct_csv is not None:
        template_featurizer = templates_custom.CustomizedComplexHitFeaturizer(
            input_atom_dir=FLAGS.struct_atom_dir,
            max_hits=template_max_hits,
            kalign_binary_path=os.path.join(FLAGS.env_dir, 'kalign'))
        custom_inputs.temp_struct_csv = FLAGS.temp_struct_csv

    if FLAGS.monomer_temp_csvs is not None:
        template_featurizer = templates_custom.CustomizedMonomerHitFeaturizer(
            input_pdb_dir=FLAGS.struct_atom_dir,
            max_hits=template_max_hits,
            kalign_binary_path=os.path.join(FLAGS.env_dir, 'kalign'))
        custom_inputs.monomer_temp_csvs = FLAGS.monomer_temp_csvs

    data_pipeline = pipeline_multimer_custom.DataPipeline(template_searcher=template_searcher,
                                                          template_featurizer=template_featurizer,
                                                          monomer_template_featurizer=monomer_template_featurizer)
    if FLAGS.msa_pair_file is not None:
        custom_inputs.msa_pair_file = FLAGS.msa_pair_file

    descs, sequences = parse_fasta(FLAGS.fasta_path)

    for i, seq in enumerate(sequences):
        seq_id = protein.PDB_CHAIN_IDS[i]
        custom_inputs.monomer_a3ms[seq_id] = FLAGS.monomer_a3ms[i]
        if custom_inputs.multimer_a3ms is not None:
            custom_inputs.multimer_a3ms[seq_id] = FLAGS.multimer_a3ms[i]
        else:
            custom_inputs.multimer_a3ms[seq_id] = FLAGS.monomer_a3ms[i]
        if FLAGS.template_stos is not None:
            custom_inputs.template_stos[seq_id] = FLAGS.template_stos[i]
        if FLAGS.template_hits_files is not None:
            custom_inputs.template_hits_files[seq_id] = FLAGS.template_hits_files[i]

    model_runners = {}
    model_names = config.MODEL_PRESETS[FLAGS.model_preset]

    num_predictions_per_model = FLAGS.num_multimer_predictions_per_model
    for model_name in model_names:
        model_config = config.model_config(model_name)
        model_config.model.num_ensemble_eval = FLAGS.multimer_num_ensemble
        model_config.model.num_recycle = FLAGS.multimer_num_recycle
        model_params = data.get_model_haiku_params(
            model_name=model_name, data_dir=FLAGS.database_dir)
        model_runner = model.RunModel(model_config, model_params)
        for i in range(num_predictions_per_model):
            model_runners[f'{model_name}_pred_{i}'] = model_runner

    logging.info('Have %d models: %s', len(model_runners),
                 list(model_runners.keys()))

    random_seed = random.randrange(sys.maxsize // len(model_runners))
    logging.info('Using random seed %d for the data pipeline', random_seed)

    # amber_relaxer = None
    amber_relaxer = relax.AmberRelaxation(
        max_iterations=RELAX_MAX_ITERATIONS,
        tolerance=RELAX_ENERGY_TOLERANCE,
        stiffness=RELAX_STIFFNESS,
        exclude_residues=RELAX_EXCLUDE_RESIDUES,
        max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS,
        use_gpu=FLAGS.use_gpu_relax)

    # print(os.path.exists(custom_inputs.temp_struct_csv))
    if os.path.exists(custom_inputs.temp_struct_csv):
        if len(custom_inputs.monomer_model_paths) > 0:
            # os.system(f"touch {FLAGS.output_dir}/structure_based_temp_with_alphafold.txt")
            logging.info('Using structure based templates with alphafold monomer models pipeline')
        else:
            # os.system(f"touch {FLAGS.output_dir}/structure_based_temp.txt")
            logging.info('Using structure based templates pipeline')

    if os.path.exists(custom_inputs.temp_seq_pair_file):
        # os.system(f"touch {FLAGS.output_dir}/sequences_based_temp.txt")
        logging.info('Using sequence based templates pipeline')

    if len(custom_inputs.template_stos) > 0:
        # os.system(f"touch {FLAGS.output_dir}/original_temp.txt")
        logging.info('Using original templates pipeline')

    fasta_name = pathlib.Path(FLAGS.fasta_path).stem
    predict_structure(
        fasta_path=FLAGS.fasta_path,
        fasta_name=fasta_name,
        output_dir=FLAGS.output_dir,
        data_pipeline=data_pipeline,
        model_runners=model_runners,
        custom_inputs=custom_inputs,
        amber_relaxer=amber_relaxer,
        benchmark=False,
        random_seed=random_seed,
        models_to_relax=FLAGS.models_to_relax)


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'fasta_path',
        'monomer_a3ms',
        'env_dir',
        'database_dir',
        'output_dir'
    ])
    app.run(main)
