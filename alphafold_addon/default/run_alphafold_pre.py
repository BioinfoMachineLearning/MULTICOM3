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
from alphafold.common import confidence
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data_custom import pipeline
from alphafold.data_custom import pipeline_multimer
from alphafold.data import templates
from alphafold.data.tools import hhsearch
from alphafold.data.tools import hmmsearch
from alphafold.model import config
from alphafold.model import data
from alphafold.model import model
from alphafold.relax import relax
import jax.numpy as jnp
import numpy as np

# Internal import (7716).

logging.set_verbosity(logging.INFO)


@enum.unique
class ModelsToRelax(enum.Enum):
  ALL = 0
  BEST = 1
  NONE = 2


flags.DEFINE_string('fasta_path', None, 'Paths to FASTA file.')
flags.DEFINE_boolean('is_prokaryote', None, 'Optional for multimer system, not used by the '
                     'single chain system. This list should contain a boolean for each fasta '
                     'specifying true where the target complex is from a prokaryote, and false '
                     'where it is not, or where the origin is unknown. These values determine '
                     'the pairing method for the MSA.')

flags.DEFINE_list('bfd_uniref_a3ms', None, 'Paths to the bfd and uniref a3ms for the subunits.')
flags.DEFINE_list('mgnify_stos', None, 'Paths to mgnify msa for the subunits.')
flags.DEFINE_list('uniref90_stos', None, 'Paths to uniref90 msa for the subunits')
flags.DEFINE_list('uniprot_stos', None, 'Paths to the uniprot msa for the subunits')

flags.DEFINE_string('output_dir', None, 'Path to a directory that will store the results.')

flags.DEFINE_string('env_dir', None, 'AlphaFold python environment directory')

flags.DEFINE_string('database_dir', None, 'AlphaFold database directory')

flags.DEFINE_string('max_template_date', None, 'Maximum template release date '
                                               'to consider. Important if folding historical test sets.')

flags.DEFINE_boolean('run_relax', True, 'Whether to run the final relaxation '
                                        'step on the predicted models. Turning relax off might '
                                        'result in predictions with distracting stereochemical '
                                        'violations but might help in case you are having issues '
                                        'with the relaxation stage.')

flags.DEFINE_boolean('use_gpu_relax', False, 'Whether to relax on GPU. '
                                             'Relax on GPU can be much faster than CPU, so it is '
                                             'recommended to enable if possible. GPUs must be available'
                                             ' if this setting is enabled.')

flags.DEFINE_integer('num_multimer_predictions_per_model', 2, 'How many '
                                                              'predictions (each with a different random seed) will be '
                                                              'generated per model. E.g. if this is 2 and there are 5 '
                                                              'models then there will be 10 predictions per input. '
                                                              'Note: this FLAG only applies if model_preset=multimer')

flags.DEFINE_integer('num_monomer_predictions_per_model', 1, 'How many '
                                                              'predictions (each with a different random seed) will be '
                                                              'generated per model. E.g. if this is 2 and there are 5 '
                                                              'models then there will be 10 predictions per input. '
                                                              'Note: this FLAG only applies if model_preset=monomer')

flags.DEFINE_enum_class('models_to_relax', ModelsToRelax.ALL, ModelsToRelax,
                        'The models to run the final relaxation step on. '
                        'If `all`, all models are relaxed, which may be time '
                        'consuming. If `best`, only the most confident model '
                        'is relaxed. If `none`, relaxation is not run. Turning '
                        'off relaxation might result in predictions with '
                        'distracting stereochemical violations but might help '
                        'in case you are having issues with the relaxation '
                        'stage.')

flags.DEFINE_integer('monomer_num_ensemble', 1, 'Number of ensemble for generating monomer models')
flags.DEFINE_integer('monomer_num_recycle', 3, 'Number of recycles for generating monomer models')
flags.DEFINE_integer('multimer_num_ensemble', 1, 'Number of ensemble for generating multimer models')
flags.DEFINE_integer('multimer_num_recycle', 20, 'Number of recycles for generating multimer models')

FLAGS = flags.FLAGS

MAX_TEMPLATE_HITS = 20
RELAX_MAX_ITERATIONS = 0
RELAX_ENERGY_TOLERANCE = 2.39
RELAX_STIFFNESS = 10.0
RELAX_EXCLUDE_RESIDUES = []
RELAX_MAX_OUTER_ITERATIONS = 3

def read_fasta(fastafile):
    seqs = {}
    current_sequence = ""
    current_id = None
    for line in open(fastafile):
        if line.startswith(">"):
            if current_id is not None:
                seqs[current_id] = current_sequence
            current_id = line.rstrip()[1:]
            current_sequence = ""
        elif not line.startswith(";"):
            current_sequence += line.rstrip()
    seqs[current_id] = current_sequence
    return seqs


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
        fasta_path: str,
        fasta_name: str,
        output_dir_base: str,
        data_pipeline: Union[pipeline.DataPipeline, pipeline_multimer.DataPipeline],
        model_runners: Dict[str, model.RunModel],
        amber_relaxer: relax.AmberRelaxation,
        benchmark: bool,
        random_seed: int,
        models_to_relax):
    """Predicts structure using AlphaFold for the given sequence."""

    logging.info('Predicting %s', fasta_name)
    timings = {}
    output_dir = output_dir_base
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    msa_output_dir = os.path.join(output_dir, 'msas')
    if not os.path.exists(msa_output_dir):
      os.makedirs(msa_output_dir)

    # Get features.
    t_0 = time.time()
    feature_dict = data_pipeline.process(
        input_fasta_path=fasta_path,
        msa_output_dir=msa_output_dir)
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
                f.write(unrelaxed_pdbs[model_name])

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

    run_multimer_system = len(read_fasta(FLAGS.fasta_path)) > 1

    if run_multimer_system:
        model_preset = 'multimer'
        template_searcher = hmmsearch.Hmmsearch(
            binary_path=FLAGS.env_dir + '/hmmsearch',
            hmmbuild_binary_path=FLAGS.env_dir + '/hmmbuild',
            database_path=FLAGS.database_dir + '/pdb_seqres/pdb_seqres.txt')

        template_featurizer = templates.HmmsearchHitFeaturizer(
            mmcif_dir=FLAGS.database_dir + '/pdb_mmcif/mmcif_files',
            max_template_date='9999-06-01',
            max_hits=MAX_TEMPLATE_HITS,
            kalign_binary_path=FLAGS.env_dir + '/kalign',
            release_dates_path=None,
            obsolete_pdbs_path=FLAGS.database_dir + '/pdb_mmcif/obsolete.dat')
    else:
        model_preset = 'monomer'
        template_searcher = hhsearch.HHSearch(
            binary_path=FLAGS.env_dir + '/hhsearch',
            databases=[FLAGS.database_dir + '/pdb70/pdb70'])

        template_featurizer = templates.HhsearchHitFeaturizer(
            mmcif_dir=FLAGS.database_dir + '/pdb_mmcif/mmcif_files',
            max_template_date='9999-06-01',
            max_hits=MAX_TEMPLATE_HITS,
            kalign_binary_path=FLAGS.env_dir + '/kalign',
            release_dates_path=None,
            obsolete_pdbs_path=FLAGS.database_dir + '/pdb_mmcif/obsolete.dat')

    monomer_data_pipeline = pipeline.DataPipeline(
        jackhmmer_binary_path=FLAGS.env_dir + '/jackhmmer',
        hhblits_binary_path=FLAGS.env_dir + '/hhblits',
        uniref90_database_path=FLAGS.database_dir + '/uniref90/uniref90.fasta',
        mgnify_database_path=FLAGS.database_dir + '/mgnify/mgy_clusters_2022_05.fa',
        bfd_database_path=FLAGS.database_dir + '/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt',
        uniref30_database_path=FLAGS.database_dir + '/uniref30/UniRef30_2021_03',
        small_bfd_database_path=FLAGS.database_dir + '/small_bfd/bfd-first_non_consensus_sequences.fasta',
        template_searcher=template_searcher,
        template_featurizer=template_featurizer,
        use_small_bfd=False,
        use_precomputed_msas=True)

    if run_multimer_system:
        num_ensemble = FLAGS.multimer_num_ensemble
        num_recycle = FLAGS.multimer_num_recycle
        num_predictions_per_model = FLAGS.num_multimer_predictions_per_model
        data_pipeline = pipeline_multimer.DataPipeline(
            monomer_data_pipeline=monomer_data_pipeline,
            jackhmmer_binary_path=FLAGS.env_dir + '/jackhmmer',
            uniprot_database_path=FLAGS.database_dir + '/uniprot/uniprot.fasta',
            use_precomputed_msas=True)
    else:
        num_ensemble = FLAGS.monomer_num_ensemble
        num_recycle = FLAGS.monomer_num_ensemble
        num_predictions_per_model = FLAGS.num_monomer_predictions_per_model
        data_pipeline = monomer_data_pipeline

    model_runners = {}
    model_names = config.MODEL_PRESETS[model_preset]
    for model_name in model_names:
        model_config = config.model_config(model_name)
        if run_multimer_system:
            model_config.model.num_ensemble_eval = num_ensemble
        else:
            model_config.data.eval.num_ensemble = num_ensemble
        model_config.model.num_recycle = num_recycle
        model_params = data.get_model_haiku_params(model_name=model_name, data_dir=FLAGS.database_dir)
        model_runner = model.RunModel(model_config, model_params)
        for i in range(num_predictions_per_model):
            model_runners[f'{model_name}_pred_{i}'] = model_runner

    logging.info('Have %d models: %s', len(model_runners),
                 list(model_runners.keys()))

    if FLAGS.run_relax:
        amber_relaxer = relax.AmberRelaxation(
            max_iterations=RELAX_MAX_ITERATIONS,
            tolerance=RELAX_ENERGY_TOLERANCE,
            stiffness=RELAX_STIFFNESS,
            exclude_residues=RELAX_EXCLUDE_RESIDUES,
            max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS,
            use_gpu=FLAGS.use_gpu_relax)
    else:
        amber_relaxer = None

    random_seed = random.randrange(sys.maxsize // len(model_runners))
    logging.info('Using random seed %d for the data pipeline', random_seed)

    fasta_name = pathlib.Path(FLAGS.fasta_path).stem
    if run_multimer_system:
        for chain_id, (i, fastaname) in zip(protein.PDB_CHAIN_IDS, enumerate(read_fasta(FLAGS.fasta_path))):
            msa_output_dir = os.path.join(FLAGS.output_dir, 'msas', chain_id)
            if not os.path.exists(msa_output_dir):
                os.makedirs(msa_output_dir)
            # copy the msas to the output directory so we don't need to run them again
            os.system(f"cp {FLAGS.mgnify_stos[i]} {msa_output_dir}/mgnify_hits.sto")
            os.system(f"cp {FLAGS.uniref90_stos[i]} {msa_output_dir}/uniref90_hits.sto")
            os.system(f"cp {FLAGS.bfd_uniref_a3ms[i]} {msa_output_dir}/bfd_uniref_hits.a3m")
            os.system(f"cp {FLAGS.uniprot_stos[i]} {msa_output_dir}/uniprot_hits.sto")
    else:
        msa_output_dir = os.path.join(FLAGS.output_dir, 'msas')
        if not os.path.exists(msa_output_dir):
            os.makedirs(msa_output_dir)
        os.system(f"cp {FLAGS.mgnify_stos[0]} {msa_output_dir}/mgnify_hits.sto")
        os.system(f"cp {FLAGS.uniref90_stos[0]} {msa_output_dir}/uniref90_hits.sto")
        os.system(f"cp {FLAGS.bfd_uniref_a3ms[0]} {msa_output_dir}/bfd_uniref_hits.a3m")


    predict_structure(
            fasta_path=FLAGS.fasta_path,
            fasta_name=fasta_name,
            output_dir_base=FLAGS.output_dir,
            data_pipeline=data_pipeline,
            model_runners=model_runners,
            amber_relaxer=amber_relaxer,
            benchmark=False,
            random_seed=random_seed,
            models_to_relax=FLAGS.models_to_relax)



if __name__ == '__main__':
  flags.mark_flags_as_required([
        'fasta_path',
        'output_dir',
        'env_dir',
        'database_dir'
  ])

  app.run(main)
