import os, sys, argparse, time
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import read_fasta, write_fasta
from bml_casp15.monomer_alignment_generation.pipeline import Monomer_alignment_generation_pipeline
from bml_casp15.complex_alignment_generation.pipeline import *
from bml_casp15.tertiary_structure_generation.pipeline import *
from bml_casp15.complex_templates_search import sequence_based_pipeline_complex_pdb, \
    sequence_based_pipeline_pdb, sequence_based_pipeline, structure_based_pipeline_v2
from bml_casp15.quaternary_structure_generation.pipeline import *
from bml_casp15.quaternary_structure_generation.iterative_search_pipeline_v0_2 import *
from bml_casp15.quaternary_structure_refinement import iterative_refine_pipeline
from absl import flags
from absl import app
import copy

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fastas')
flags.DEFINE_string('model_dir', None, 'Output directory')
flags.DEFINE_string('output_dir', None, 'Output directory')
FLAGS = flags.FLAGS

PDB_CHAIN_IDS = 'BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
PDB_MAX_CHAINS = len(PDB_CHAIN_IDS)  # := 62.


class FastaChain:
    sequence: str
    description: str


def _make_chain_id_map(sequences, descriptions):
    """Makes a mapping from PDB-format chain ID to sequence and description."""
    if len(sequences) != len(descriptions):
        raise ValueError('sequences and descriptions must have equal length. '
                         f'Got {len(sequences)} != {len(descriptions)}.')
    if len(sequences) > PDB_MAX_CHAINS:
        raise ValueError('Cannot process more chains than the PDB format supports. '
                         f'Got {len(sequences)} chains.')
    chain_id_map = {}
    chain_id_seq_map = {}
    for chain_id, sequence, description in zip(
            PDB_CHAIN_IDS, sequences, descriptions):
        chain_id_map[chain_id] = FastaChain(sequence=sequence, description=description)
        chain_id_seq_map[chain_id] = sequence
    return chain_id_map, chain_id_seq_map


def run_monomer_msa_pipeline(fasta, outdir, params):
    uniref30 = params['uniref_db_dir'] + '/' + params['uniref_db']
    uniclust30 = params['uniclust_db_dir'] + '/' + params['uniclust_db']
    uniref90_fasta = params['uniref90_fasta']
    uniprot_fasta = params['uniprot_fasta']
    smallbfd = ""  # params['smallbfd_database']
    bfd = params['bfd_database']
    mgnify = params['mgnify_database']
    hhblits_binary = params['hhblits_program']
    jackhmmer_binary = params['jackhmmer_program']
    result = None
    try:
        monomer_msa_pipeline = Monomer_alignment_generation_pipeline(jackhmmer_binary,
                                                                     hhblits_binary,
                                                                     uniref90_fasta,
                                                                     mgnify,
                                                                     smallbfd,
                                                                     bfd,
                                                                     uniref30,
                                                                     uniclust30,
                                                                     uniprot_fasta)

        result = monomer_msa_pipeline.process(fasta, outdir)
    except Exception as e:
        print(e)
    return result


def run_monomer_template_search_pipeline(targetname, sequence, a3m, outdir, params):
    template_file = None
    try:
        pipeline = monomer_sequence_based_template_search_pipeline(params)
        template_file = pipeline.search(targetname=targetname, sequence=sequence, a3m=a3m, outdir=outdir)
    except Exception as e:
        print(e)
    return template_file


def run_monomer_evaluation_pipeline(targetname, fasta_file, input_monomer_dir,
                                    chainid, input_multimer_dir, outputdir):
    makedir_if_not_exists(outputdir)
    result = None
    pipeline = Monomer_structure_evaluation_pipeline(params=params,
                                                     run_methods=["apollo", "alphafold", "enQA"],
                                                     use_gpu=True)
    try:
        result = pipeline.process(targetname=targetname, fasta_file=fasta_file,
                                  monomer_model_dir=input_monomer_dir, multimer_model_dir=input_multimer_dir,
                                  output_dir=outputdir, chainid_in_multimer=chainid)
    except Exception as e:
        print(e)
    return result


def run_concatenate_dimer_msas_pipeline(multimer, monomer_aln_dir, outputdir, params):
    chains = multimer.split('_')
    # alignment = {'outdir': f"{outputdir}/{'_'.join(chains)}"}
    alignment = {'outdir': outputdir}
    for i in range(len(chains)):
        chain = chains[i]
        chain_aln_dir = monomer_aln_dir + '/' + chain
        if os.path.exists(chain_aln_dir):
            chain_a3ms = {'name': chain,
                          'uniref30_a3m': f"{chain_aln_dir}/{chain}_uniref30.a3m",
                          'uniref90_sto': f"{chain_aln_dir}/{chain}_uniref90.sto",
                          'uniprot_sto': f"{chain_aln_dir}/{chain}_uniprot.sto",
                          'uniclust30_a3m': f"{chain_aln_dir}/{chain}_uniclust30.a3m"}
        else:
            chain_a3ms = {'name': chain}
        alignment[f"chain{i + 1}"] = chain_a3ms

    complete = True
    for name in alignment:
        if name == 'outdir':
            continue
        a3ms = alignment[name]
        if len(a3ms) == 1:
            complete = False
        for key in a3ms:
            if key.find('uni') >= 0:
                if not os.path.exists(a3ms[key]):
                    complete = False
                else:
                    contents = open(a3ms[key]).readlines()
                    if len(contents) == 0:
                        print(f"File: {a3ms[key]} is empty!")
                        complete = False
                        os.system(f"rm {a3ms[key]}")

    if complete:
        alignments = [alignment]
        print(f"Total {len(alignments)} pairs can be concatenated")
        print("Start to concatenate alignments for dimers")
        complex_alignment_concatenation_pipeline = Complex_alignment_concatenation_pipeline(params)
        alignments = complex_alignment_concatenation_pipeline.concatenate(alignments, params['hhfilter_program'])
    else:
        print("The a3ms for dimers are not complete!")


def run_complex_template_search_pipeline(multimers, monomer_aln_dir, monomer_model_dir, outdir, params):
    monomer_template_inputs = []
    for chain in multimers:
        chain_template_a3m = f"{monomer_aln_dir}/{chain}/{chain}_uniref90.sto"
        if not os.path.exists(chain1_template_a3m):
            raise Exception(f"Cannot find uniref90 alignment for {chain}")
        seq = open(f"{monomer_aln_dir}/{chain}/{chain}.fasta").readlines()[1].rstrip('\n')
        monomer_template_input = monomer_template_input(name=chain, msa_path=chain_template_a3m,
                                                        hmm_path="", seq=seq)
        monomer_template_inputs += [monomer_template_input]

    pdb_seq_dir = outdir + '/pdb_seq'
    makedir_if_not_exists(pdb_seq_dir)
    if not os.path.exists(pdb_seq_dir + '/concatenated_template_idx.csv'):
        pipeline = sequence_based_pipeline_pdb.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, pdb_seq_dir)

    complex_pdb_seq_dir = outdir + '/complex_pdb_seq'
    makedir_if_not_exists(complex_pdb_seq_dir)
    if not os.path.exists(complex_pdb_seq_dir + '/concatenated_template_idx.csv'):
        pipeline = sequence_based_pipeline_complex_pdb.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, complex_pdb_seq_dir)

    pdb70_seq_dir = outdir + '/pdb70_seq'
    makedir_if_not_exists(pdb70_seq_dir)
    if not os.path.exists(pdb70_seq_dir + '/concatenated_template_idx.csv'):
        pipeline = sequence_based_pipeline.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, pdb70_seq_dir)

    struct_temp_dir = outdir + '/struct_temp'
    makedir_if_not_exists(struct_temp_dir)
    if not os.path.exists(struct_temp_dir + '/concatenated_template_idx.csv'):
        monomer_pdbs = []
        for chain in multimers:
            monomer_pdb = f"{monomer_model_dir}/{chain}/default/ranked_0.pdb"
            if not os.path.exists(monomer_pdb):
                print(f"Cannot find teritary structure for {chain}: {monomer_pdb}")
                continue
            monomer_pdbs += [monomer_pdb]
        pipeline = structure_based_pipeline_v2.Complex_structure_based_template_search_pipeline(params)
        pipeline.search(monomer_pdbs, struct_temp_dir)


def read_qa_txt(infile):
    models = []
    scores = []
    for line in open(infile):
        line = line.rstrip('\n')
        contents = line.split()
        if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[0] == "QMODE" or \
                contents[0] == "END":
            continue
        model, score = line.split()
        models += [model]
        scores += [float(score)]
    df = pd.DataFrame({'model': models, 'score': scores})
    df = df.sort_values(by=['score'], ascending=False)
    df.reset_index(inplace=True)
    return df


def complete_result(outputdir):
    complete = True
    for i in range(0, 5):
        model = f'{outputdir}/ranked_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def parse_fasta(fasta_string):
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


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(fasta_path)

    makedir_if_not_exists(FLAGS.output_dir)

    with open(fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = _make_chain_id_map(sequences=input_seqs,
                                                        descriptions=input_descs)

    # N1_outdir = FLAGS.output_dir + '/N1_monomer_alignments_generation'
    # N1_outdir_img = FLAGS.output_dir + '/N1_monomer_alignments_generation_img'
    # N2_outdir = FLAGS.output_dir + '/N2_monomer_template_search'
    # N3_outdir = FLAGS.output_dir + '/N3_monomer_structure_generation'
    # N4_outdir = FLAGS.output_dir + '/N4_monomer_structure_evaluation'
    # img_msas = {}
    #
    # print("#################################################################################################")
    #
    # print("#################################################################################################")
    # print("1-4. Start to generate monomer models")
    #
    # makedir_if_not_exists(N1_outdir)

    # monomer_qas_res = {}
    # for chain_id in chain_id_map:
    #     monomer_id = chain_id_map[chain_id].description
    #     monomer_sequence = chain_id_map[chain_id].sequence
    #
    #     with open(f"{FLAGS.output_dir}/{monomer_id}.fasta", "w") as fw:
    #         write_fasta({monomer_id: monomer_sequence}, fw)
    #     N1_monomer_outdir = N1_outdir + '/' + monomer_id
    #     makedir_if_not_exists(N1_monomer_outdir)
    #     result = run_monomer_msa_pipeline([f"{FLAGS.output_dir}/{monomer_id}.fasta", N1_monomer_outdir, params])
    #     if result is None:
    #         raise RuntimeError(f"Program failed in step 1: monomer {monomer_id} alignment generation")
    #
    #     N1_monomer_outdir_img = N1_outdir_img + '/' + monomer_id
    #     makedir_if_not_exists(N1_monomer_outdir_img)
    #     monomer_msa_pipeline_img = Monomer_alignment_generation_pipeline_img(binary_path=params['deepmsa2_program'],
    #                                                                          bfd_database_path=params['bfd_database'],
    #                                                                          img_database_path=params['img_database'],
    #                                                                          metaclust_database_path=params[
    #                                                                              'metaclust_database'],
    #                                                                          mgnify_database_path=params[
    #                                                                              'mgnify_database'],
    #                                                                          uniref90_database_path=params[
    #                                                                              'uniref90_database'])
    #     img_msas[chain_id] = monomer_msa_pipeline_img.process(f"{FLAGS.output_dir}/{monomer_id}.fasta",
    #                                                           N1_monomer_outdir_img)
    #
    #     N2_monomer_outdir = N2_outdir + '/' + monomer_id
    #     makedir_if_not_exists(N2_monomer_outdir)
    #     template_file = run_monomer_template_search_pipeline(targetname=monomer_id, sequence=monomer_id,
    #                                                          a3m=f"{N1_monomer_outdir}/{monomer_id}_uniref90.sto",
    #                                                          outdir=N2_monomer_outdir, params=params)
    #     if template_file is None:
    #         raise RuntimeError(f"Program failed in step 2: monomer {monomer_id} template search")
    #
    #     N3_monomer_outdir = N3_outdir + '/N3_monomer_structure_generation'
    #     makedir_if_not_exists(N3_monomer_outdir)
    #     try:
    #         pipeline = Monomer_tertiary_structure_prediction_pipeline(params)
    #         pipeline.process_single(fasta_path=f"{FLAGS.output_dir}/{monomer_id}.fasta",
    #                                 alndir=N1_monomer_outdir,
    #                                 template_dir=N2_monomer_outdir,
    #                                 outdir=N3_monomer_outdir,
    #                                 run_methods=['default', 'original', 'colabfold', 'seq_temp_pdb'])
    #     except Exception as e:
    #         print(e)
    #         print(f"Program failed in step 3: monomer {monomer_id} structure generation")
    #
    #     N4_monomer_outdir = N4_outdir + '/' + monomer_id
    #     makedir_if_not_exists(N4_monomer_outdir)
    #     result = run_monomer_evaluation_pipeline(targetname=monomer_id,
    #                                              fasta_file=f"{FLAGS.output_dir}/{monomer_id}.fasta",
    #                                              input_monomer_dir=N3_monomer_outdir,
    #                                              outputdir=N4_monomer_outdir)
    #     if result is None:
    #         raise RuntimeError(f"Program failed in step 4: monomer {monomer_id} model evaluation")
    #
    #     monomer_qas_res[monomer_id] = result
    #
    # print("#################################################################################################")
    #
    # print("#################################################################################################")
    # print("5. Start to generate complex alignments")
    #
    # N5_outdir = FLAGS.output_dir + '/N5_complex_alignments_concatenation'
    # makedir_if_not_exists(N5_outdir)
    #
    # try:
    #     run_concatenate_dimer_msas_pipeline(
    #         multimer='_'.join([chain_id_map[chain_id].description for chain_id in chain_id_map]),
    #         monomer_aln_dir=N1_outdir, outputdir=N5_outdir, params=params)
    # except Exception as e:
    #     print(e)
    #     print("Program failed in step 5")
    #
    # print("#################################################################################################")
    #
    # print("#################################################################################################")
    #
    # print("6. Start to search complex templates based on monomer structures")
    #
    # N6_outdir = FLAGS.output_dir + '/N6_complex_templates_search'
    #
    # run_complex_template_search_pipeline(multimers=[chain_id_map[chain_id].description for chain_id in chain_id_map],
    #                                      monomer_aln_dir=N1_outdir,
    #                                      monomer_model_dir=N3_outdir,
    #                                      outdir=N6_outdir, params=params)
    #
    # print("#################################################################################################")
    #
    # print("#################################################################################################")
    #
    # print("7. Start to generate complex quaternary structures")
    # N7_outdir = FLAGS.output_dir + '/N7_quaternary_structure_generation'
    # makedir_if_not_exists(N7_outdir)
    #
    # try:
    #     pipeline = Quaternary_structure_prediction_pipeline(params)
    #     result = pipeline.process(fasta_path=fasta_path,
    #                               aln_dir=N1_outdir,
    #                               complex_aln_dir=N2_outdir,
    #                               template_dir=N5_outdir,
    #                               monomer_model_dir=N3_outdir,
    #                               output_dir=N7_outdir)
    # except Exception as e:
    #     print(e)
    #     print("Program failed in step 7")
    #
    # pipeline = Multimer_iterative_generation_pipeline_monomer(params)
    # for i in range(5):
    #     monomer_pdb_dirs = {}
    #     monomer_alphafold_a3ms = {}
    #     for chain_id in chain_id_map:
    #         monomer_id = chain_id_map[chain_id].description
    #         monomer_ranking = read_qa_txt(monomer_qas_res[monomer_id]['apollo'])
    #         pdb_name = monomer_ranking.loc[i, 'model']
    #         monomer_pdb_dirs[chain_id] = f"{N4_outdir}/{monomer_id}/pdb/{pdb_name}"
    #         monomer_alphafold_a3ms[chain_id] = f"{N4_outdir}/{monomer_id}/msa/{pdb_name.replace('.pdb', '.a3m')}"
    #
    #     pipeline.search_single(fasta_file=fasta_path,
    #                            chain_id_map=chain_id_map,
    #                            monomer_pdb_dirs=monomer_pdb_dirs,
    #                            monomer_alphafold_a3ms=monomer_alphafold_a3ms,
    #                            outdir=f"{N7_outdir}/iter_{i + 1}")
    #
    # print("Complex quaternary structure generation has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("8. Start to evaluate multimer models")

    N7_outdir = FLAGS.model_dir

    N8_outdir = outdir + '/N8_multimer_structure_evaluation'
    pipeline = Quaternary_structure_evaluation_pipeline(params=params)
    multimer_qa_result = pipeline.process(N7_outdir, N8_outdir)

    print("#################################################################################################")

    print("#################################################################################################")

    print("9. Start to refine multimer models based on the qa rankings")

    N9_outdir = outdir + '/N9_multimer_structure_refinement'

    makedir_if_not_exists(N9_outdir)

    ref_ranking = read_qa_txt(multimer_qa_result['alphafold'])  # apollo or average ranking or the three qas

    pipeline = iterative_refine_pipeline.Multimer_iterative_refinement_pipeline_server(params=params)
    refine_inputs = []
    for i in range(5):
        pdb_name = ref_ranking.loc[i, 'model']
        msa_paths = {}
        for chain_id in chain_id_map:
            msa_paths[chain_id] = dict(
                paired_msa=f"{N8_outdir}/msa/{chain_id_map[chain_id].description}/{pdb_name}.multimer.a3m",
                monomer_msa=f"{N8_outdir}/msa/{chain_id_map[chain_id].description}/{pdb_name}.monomer.a3m")

        refine_input = iterative_refine_pipeline.refinement_input(chain_id_map=chain_id_map,
                                                                  fasta_path=fasta_path,
                                                                  pdb_path=N8_outdir + '/pdb/' + pdb_name,
                                                                  pkl_path=N8_outdir + '/pkl/' + pdb_name.replace(
                                                                      '.pdb', '.pkl'),
                                                                  msa_paths=msa_paths)
        refine_inputs += [refine_input]
    pipeline.search(refinement_inputs=refine_inputs, outdir=N9_outdir)

    final_dir = N9_outdir + '/final'
    makedir_if_not_exists(final_dir)

    pipeline = iterative_refine_pipeline.Multimer_refinement_model_selection()
    pipeline.select_v1(indir=N9_outdir, outdir=final_dir)

    print("The refinement for the top-ranked multimer models has been finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
