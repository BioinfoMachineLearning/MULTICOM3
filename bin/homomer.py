import os, sys, argparse, time
from multiprocessing import Pool
from multicom3.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from multicom3.monomer_alignment_generation.alignment import write_fasta
from multicom3.common.protein import read_qa_txt_as_df, parse_fasta, complete_result, make_chain_id_map
from multicom3.multimer_structure_refinement import iterative_refine_pipeline_multimer
from multicom3.monomer_structure_refinement import iterative_refine_pipeline
from multicom3.common.pipeline import run_monomer_msa_pipeline, run_monomer_template_search_pipeline, \
    run_monomer_structure_generation_pipeline_v2, run_monomer_evaluation_pipeline, run_monomer_refinement_pipeline, \
    run_monomer_msas_concatenation_pipeline, run_monomer_templates_concatenation_pipeline, \
    run_multimer_structure_generation_homo_pipeline_v2, \
    run_multimer_structure_generation_pipeline_foldseek, run_multimer_structure_generation_pipeline_foldseek_old, \
    run_multimer_refinement_pipeline, run_multimer_evaluation_pipeline, run_monomer_msa_pipeline_img, \
    foldseek_iterative_monomer_input, copy_same_sequence_msas, run_multimer_structure_generation_homo_pipeline_img_v2

from absl import flags
from absl import app
import copy
import pandas as pd

flags.DEFINE_string('option_file', None, 'option file')
flags.DEFINE_string('fasta_path', None, 'Path to multimer fasta')
flags.DEFINE_string('output_dir', None, 'Output directory')
flags.DEFINE_boolean('run_img', False, 'Whether to use IMG alignment to generate models')
flags.DEFINE_boolean('run_refinement', True, 'Whether run model refinement')
FLAGS = flags.FLAGS


def main(argv):
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')

    os.environ['TF_FORCE_UNIFIED_MEMORY'] = '1'
    os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '4.0'
    
    FLAGS.fasta_path = os.path.abspath(FLAGS.fasta_path)
    FLAGS.output_dir = os.path.abspath(FLAGS.output_dir)
    
    check_file(FLAGS.option_file)

    params = read_option_file(FLAGS.option_file)

    makedir_if_not_exists(FLAGS.output_dir)

    check_dirs(params, ['hhblits_program', 'jackhmmer_program'], isdir=False)

    check_file(FLAGS.fasta_path)

    makedir_if_not_exists(FLAGS.output_dir)

    N1_outdir = os.path.join(FLAGS.output_dir, 'N1_monomer_alignments_generation')
    N1_outdir_img = os.path.join(FLAGS.output_dir, 'N1_monomer_alignments_generation_img') 
    N2_outdir = os.path.join(FLAGS.output_dir, 'N2_monomer_template_search')
    N3_outdir = os.path.join(FLAGS.output_dir, 'N3_monomer_structure_generation')
    img_msas = {}

    print("#################################################################################################")

    print("#################################################################################################")
    print("1-3. Start to generate monomer models")

    makedir_if_not_exists(N1_outdir)

    with open(FLAGS.fasta_path) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    chain_id_map, chain_id_seq_map = make_chain_id_map(sequences=input_seqs,
                                                       descriptions=input_descs)

    processed_seuqences = {}
    monomer_qas_res = {}
    for chain_id in chain_id_map:
        monomer_id = chain_id
        monomer_sequence = chain_id_map[chain_id].sequence
        monomer_fasta = os.path.join(FLAGS.output_dir, monomer_id + ".fasta")

        if monomer_sequence not in processed_seuqences:

            with open(monomer_fasta, "w") as fw:
                write_fasta({chain_id: monomer_sequence}, fw)
            N1_monomer_outdir = os.path.join(N1_outdir, monomer_id)
            makedir_if_not_exists(N1_monomer_outdir)
            result = run_monomer_msa_pipeline(monomer_fasta, N1_monomer_outdir, params)
            if result is None:
                raise RuntimeError(f"Program failed in step 1: monomer {monomer_id} alignment generation")

            if FLAGS.run_img:
                N1_monomer_outdir_img = os.path.join(N1_outdir_img, monomer_id)
                makedir_if_not_exists(N1_monomer_outdir_img)
                img_msas[chain_id] = run_monomer_msa_pipeline_img(params=params,
                                                                fasta=monomer_fasta,
                                                                outdir=N1_monomer_outdir_img)

            N2_monomer_outdir = os.path.join(N2_outdir, monomer_id)
            makedir_if_not_exists(N2_monomer_outdir)
            template_file = run_monomer_template_search_pipeline(targetname=monomer_id, sequence=monomer_id,
                                                                 a3m=os.path.join(N1_monomer_outdir, monomer_id + "_uniref90.sto"),
                                                                 outdir=N2_monomer_outdir, params=params)
            if template_file is None:
                raise RuntimeError(f"Program failed in step 2: monomer {monomer_id} template search")

            N3_monomer_outdir = os.path.join(N3_outdir, monomer_id)
            makedir_if_not_exists(N3_monomer_outdir)
            if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                                fasta_path=monomer_fasta,
                                                                alndir=N1_monomer_outdir,
                                                                templatedir=N2_monomer_outdir,
                                                                outdir=N3_monomer_outdir):
                print(f"Program failed in step 3: monomer {monomer_id} structure generation")

            processed_seuqences[monomer_sequence] = monomer_id

        else:
            with open(monomer_fasta, "w") as fw:
                write_fasta({chain_id: monomer_sequence}, fw)
            N1_monomer_outdir = os.path.join(N1_outdir, monomer_id)
            makedir_if_not_exists(N1_monomer_outdir)

            copy_same_sequence_msas(srcdir=os.path.join(N1_outdir, processed_seuqences[monomer_sequence]),
                                    trgdir=N1_monomer_outdir,
                                    srcname=processed_seuqences[monomer_sequence],
                                    trgname=monomer_id)

            N1_monomer_outdir_img = os.path.join(N1_outdir_img, monomer_id)
            makedir_if_not_exists(N1_monomer_outdir_img)

            N2_monomer_outdir = os.path.join(N2_outdir, monomer_id)
            makedir_if_not_exists(N2_monomer_outdir)

            # make a copy
            N3_monomer_outdir = os.path.join(N3_outdir, monomer_id)

            if not os.path.exists(N3_monomer_outdir):
                os.system("cp -r " + os.path.join(N3_outdir, processed_seuqences[monomer_sequence]) + " " + N3_monomer_outdir)

    print("#################################################################################################")

    print("#################################################################################################")

    print("#################################################################################################")
    print("4. Start to generate complex alignments")

    N4_outdir = os.path.join(FLAGS.output_dir, 'N4_monomer_alignments_concatenation')
    makedir_if_not_exists(N4_outdir)

    try:
        concat_methods = ['pdb_interact', 'species_interact', 'uniclust_oxmatch']
        run_monomer_msas_concatenation_pipeline(
            multimer=','.join([chain_id for chain_id in chain_id_map]),
            run_methods=concat_methods,
            monomer_aln_dir=N1_outdir, outputdir=N4_outdir, params=params, is_homomers=True)

    except Exception as e:
        print(e)
        print("Program failed in step 4")

    print("#################################################################################################")

    print("#################################################################################################")

    print("5. Start to search complex templates based on monomer structures")

    N5_outdir = os.path.join(FLAGS.output_dir, 'N5_monomer_templates_concatenation')

    run_monomer_templates_concatenation_pipeline(multimers=[chain_id for chain_id in chain_id_map],
                                         monomer_aln_dir=N1_outdir,
                                         monomer_model_dir=N3_outdir,
                                         outdir=N5_outdir, params=params)

    print("#################################################################################################")

    print("#################################################################################################")

    print("6. Start to generate complex multimer structures")
    N6_outdir = os.path.join(FLAGS.output_dir, 'N6_multimer_structure_generation')
    makedir_if_not_exists(N6_outdir)

    if not run_multimer_structure_generation_homo_pipeline_v2(params=params,
                                                                fasta_path=FLAGS.fasta_path,
                                                                chain_id_map=chain_id_map,
                                                                aln_dir=N1_outdir,
                                                                complex_aln_dir=N4_outdir,
                                                                template_dir=N5_outdir,
                                                                monomer_model_dir=N3_outdir,
                                                                output_dir=N6_outdir):
        print("Program failed in step 6")

    print("Multimer structure generation has been finished!")

    print("#################################################################################################")

    if FLAGS.run_img:

        processed_seuqences = {}
        for chain_id in chain_id_map:
            monomer_id = chain_id
            monomer_sequence = chain_id_map[chain_id].sequence
            print(monomer_id)
            if monomer_sequence not in processed_seuqences:
                N3_monomer_outdir = os.path.join(N3_outdir, monomer_id)
                default_alphafold_msa = os.path.join(N3_monomer_outdir, 'default/msas/monomer_final.a3m')
                
                while not os.path.exists(img_msas[chain_id]):
                    # sleep for 5 mins
                    print("waiting for img alignment")
                    time.sleep(300)
                print("Found img alignment, start to run monomer model generation again")
                os.system("cp " + os.path.join(N1_outdir, monomer_id, monomer_id + "_uniref90.sto") + " " + os.path.join(N1_outdir + "_img", monomer_id))
                new_contents = []
                for line_idx, line in enumerate(open(img_msas[chain_id]).readlines()):
                    if line_idx == 0:
                        new_contents.append(f">{monomer_id}\n")
                    else:
                        new_contents.append(line)
                open(img_msas[chain_id], 'w').writelines(new_contents)
                if not run_monomer_structure_generation_pipeline_v2(params=params,
                                                                    run_methods=['img', 'img+seq_template'],
                                                                    fasta_path=os.path.join(FLAGS.output_dir, monomer_id + ".fasta"),
                                                                    alndir=os.path.join(N1_outdir + '_img', monomer_id),
                                                                    templatedir=os.path.join(N2_outdir, monomer_id),
                                                                    outdir=os.path.join(N3_outdir, monomer_id)):
                    print("Program failed in step 3: monomer structure generation")
             
                processed_seuqences[monomer_sequence] = monomer_id
            
            else:
                copy_same_sequence_msas(srcdir=os.path.join(N1_outdir + '_img', processed_seuqences[monomer_sequence]),
                                        trgdir=os.path.join(N1_outdir + '_img', monomer_id),
                                        srcname=processed_seuqences[monomer_sequence],
                                        trgname=monomer_id)

                # make a copy
                N3_monomer_outdir = os.path.join(N3_outdir, monomer_id)

                if os.path.exists(N3_monomer_outdir):
                    os.system(f"rm -rf {N3_monomer_outdir}")

                makedir_if_not_exists(N3_monomer_outdir)
                os.system("cp -r " + os.path.join(N3_outdir, processed_seuqences[monomer_sequence]) + " " + N3_monomer_outdir)

            if not run_multimer_structure_generation_homo_pipeline_img_v2(params=params,
                                                                            fasta_path=FLAGS.fasta_path,
                                                                            chain_id_map=chain_id_map,
                                                                            aln_dir=N1_outdir + '_img',
                                                                            output_dir=N6_outdir):
                print("Program failed in step 6")
        print("#################################################################################################")

    print("7. Start to evaluate monomer models")

    N7_outdir = os.path.join(FLAGS.output_dir, 'N7_monomer_structure_evaluation')

    processed_seuqences = {}
    for chain_id in chain_id_map:
        monomer_id = chain_id
        monomer_sequence = chain_id_map[chain_id].sequence
        if monomer_sequence not in processed_seuqences:
            N7_monomer_outdir = os.path.join(N7_outdir, monomer_id)
            makedir_if_not_exists(N7_monomer_outdir)
            result = run_monomer_evaluation_pipeline(params=params,
                                                     targetname=monomer_id,
                                                     fasta_file=os.path.join(FLAGS.output_dir, f"{monomer_id}.fasta"),
                                                     input_monomer_dir=os.path.join(N3_outdir, monomer_id),
                                                     input_multimer_dir=N6_outdir,
                                                     outputdir=N7_monomer_outdir, generate_egnn_models=True)
            if result is None:
                raise RuntimeError(f"Program failed in step 7: monomer {monomer_id} model evaluation")
            monomer_qas_res[monomer_id] = result

            processed_seuqences[monomer_sequence] = monomer_id

        else:
            # make a copy
            N7_monomer_outdir = os.path.join(N7_outdir, monomer_id)
            os.system("cp -r " + os.path.join(N7_outdir, processed_seuqences[monomer_sequence]) + " " + N7_monomer_outdir)
            for msa in os.listdir(os.path.join(N7_monomer_outdir, 'msa')):
                os.system(f"sed -i 's/>{processed_seuqences[monomer_sequence]}/>{monomer_id}/g' " + os.path.join(N7_monomer_outdir, 'msa', msa))

    print("#################################################################################################")

    print("8. Start to run multimer iterative generation pipeline using top-ranked monomer models")

    qa_result_dir = N7_outdir

    iterative_prepare_dir = os.path.join(qa_result_dir, 'iter_prepare')
    makedir_if_not_exists(iterative_prepare_dir)

    pipeline_inputs = []
    for i in range(2):
        monomer_pdb_dirs = {}
        monomer_alphafold_a3ms = {}
        pdb_name = None

        first_monomer_id = ""
        for chain_id in chain_id_map:
            first_monomer_id = chain_id
            monomer_ranking = pd.read_csv(monomer_qas_res[first_monomer_id]['apollo_monomer'])
            pdb_name = monomer_ranking.loc[i, 'model']
            break

        current_work_dir = os.path.join(iterative_prepare_dir, str(i))
        makedir_if_not_exists(current_work_dir)

        for chain_id in chain_id_map:
            monomer_id = chain_id
            chain_pdb_dir = os.path.join(current_work_dir, monomer_id)
            makedir_if_not_exists(chain_pdb_dir)

            first_monomer_id_dir = os.path.join(qa_result_dir, first_monomer_id)
            os.system("cp " + os.path.join(first_monomer_id_dir, 'pdb', pdb_name) + " " + os.path.join(chain_pdb_dir, pdb_name))

            new_contents = []
            for idx, line in enumerate(open(os.path.join(qa_result_dir, first_monomer_id, 'msa', pdb_name.replace('.pdb', '.a3m'))).readlines()):
                if idx == 0:
                    new_contents += [f">{monomer_id}\n"]
                else:
                    new_contents += [line]

            open(os.path.join(chain_pdb_dir, pdb_name.replace('.pdb', '.a3m')), 'w').writelines(new_contents)
            monomer_pdb_dirs[chain_id] = os.path.join(chain_pdb_dir, pdb_name)
            monomer_alphafold_a3ms[chain_id] = os.path.join(chain_pdb_dir, pdb_name.replace('.pdb', '.a3m'))

        print(monomer_alphafold_a3ms)
        pipeline_inputs += [foldseek_iterative_monomer_input(monomer_pdb_dirs=monomer_pdb_dirs,
                                                             monomer_alphafold_a3ms=monomer_alphafold_a3ms)]

    if not run_multimer_structure_generation_pipeline_foldseek_old(params=params, fasta_path=FLAGS.fasta_path,
                                                                     chain_id_map=chain_id_map,
                                                                     pipeline_inputs=pipeline_inputs, outdir=N6_outdir,
                                                                     is_homomers=True):
        print("Program failed in step 6 iterative")

    if len(chain_id_map) <= 6:
        if not run_multimer_structure_generation_pipeline_foldseek(params=params, fasta_path=FLAGS.fasta_path,
                                                                     chain_id_map=chain_id_map,
                                                                     pipeline_inputs=[pipeline_inputs[0]],
                                                                     outdir=N6_outdir,
                                                                     is_homomers=True):
            print("Program failed in step 6 iterative")

    print("Multimer structure generation has been finished!")

    print("#################################################################################################")

    print("#################################################################################################")

    print("9. Start to evaluate multimer models")

    N8_outdir = os.path.join(FLAGS.output_dir, 'N8_multimer_structure_evaluation')
    multimer_qa_result = run_multimer_evaluation_pipeline(fasta_path=FLAGS.fasta_path,
                                                          params=params,
                                                          chain_id_map=chain_id_map,
                                                          indir=N6_outdir, outdir=N8_outdir, is_homomer=True)

    print("#################################################################################################")

    print("#################################################################################################")

    print("10. Start to refine multimer models based on the qa rankings")

    # if len(chain_id_map) <= 5:
    if FLAGS.run_refinement:
        
        N9_outdir = os.path.join(FLAGS.output_dir, 'N9_multimer_structure_refinement')

        makedir_if_not_exists(N9_outdir)
        ref_ranking = pd.read_csv(multimer_qa_result['pairwise_af_avg'])  # apollo or average ranking or the three qas

        refine_inputs = []
        for i in range(5):
            pdb_name = ref_ranking.loc[i, 'model']
            msa_paths = {}
            for chain_id in chain_id_map:
                msa_paths[chain_id] = dict(paired_msa=os.path.join(N8_outdir, 'msa', chain_id, pdb_name.replace('.pdb', '') + ".paired.a3m"),
                                           monomer_msa=os.path.join(N8_outdir, 'msa', chain_id, pdb_name.replace('.pdb', '') + ".monomer.a3m"))

            refine_input = iterative_refine_pipeline_multimer.refinement_input_multimer(chain_id_map=chain_id_map,
                                                                                        fasta_path=FLAGS.fasta_path,
                                                                                        pdb_path=os.path.join(N8_outdir, 'pdb', pdb_name),
                                                                                        pkl_path=os.path.join(N8_outdir, 'pkl', pdb_name.replace('.pdb', '.pkl')),
                                                                                        msa_paths=msa_paths)
            refine_inputs += [refine_input]

        final_dir = N9_outdir + '_final'
        run_multimer_refinement_pipeline(chain_id_map=chain_id_map,
                                         params=params, refinement_inputs=refine_inputs, outdir=N9_outdir,
                                         finaldir=final_dir, stoichiometry="homomer")

        print("The refinement for the top-ranked multimer models has been finished!")


if __name__ == '__main__':
    flags.mark_flags_as_required([
        'option_file',
        'fasta_path',
        'output_dir'
    ])
    app.run(main)
