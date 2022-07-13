import os, sys, argparse, time, copy
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import read_fasta, write_fasta
from bml_casp15.monomer_alignment_generation.pipeline import *
from bml_casp15.monomer_structure_generation.pipeline import *
from bml_casp15.monomer_structure_generation.pipeline_v2 import *
from bml_casp15.monomer_structure_evaluation.pipeline_sep import *
from bml_casp15.monomer_structure_evaluation.human_pipeline import *
from bml_casp15.monomer_templates_search.sequence_based_pipeline_pdb import *
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from bml_casp15.quaternary_structure_refinement import iterative_refine_pipeline_multimer
from bml_casp15.complex_alignment_generation.pipeline_v2 import *
from bml_casp15.complex_templates_search import sequence_based_pipeline_complex_pdb, \
    sequence_based_pipeline_pdb, sequence_based_pipeline, structure_based_pipeline_v2
# from bml_casp15.quaternary_structure_generation.pipeline import *
from bml_casp15.quaternary_structure_generation.pipeline_v2 import *
from bml_casp15.quaternary_structure_generation.pipeline_default import *
# from bml_casp15.quaternary_structure_generation.pipeline_homo import *
from bml_casp15.quaternary_structure_generation.pipeline_homo_v2 import *
from bml_casp15.quaternary_structure_generation.iterative_search_pipeline_v0_2 import *
from bml_casp15.quaternary_structure_generation.iterative_search_pipeline_v0_2_old import *
from bml_casp15.quaternary_structure_evaluation.pipeline import *
from bml_casp15.quaternary_structure_evaluation.human_pipeline import *
from bml_casp15.common.protein import *
import pandas as pd
import numpy as np
from bml_casp15.monomer_structure_refinement.util import cal_tmscore
from bml_casp15.monomer_structure_evaluation.alphafold_ranking import Alphafold_pkl_qa


def run_monomer_msa_pipeline(fasta, outdir, params, only_monomer=False):
    uniref30 = params['uniref_db_dir'] + '/' + params['uniref_db']
    uniref30_new = params['uniref_db_dir_new'] + '/' + params['uniref_db_new']
    uniclust30 = params['uniclust_db_dir'] + '/' + params['uniclust_db']
    uniref90_fasta = params['uniref90_fasta']

    smallbfd = ""  # params['smallbfd_database']
    bfd = params['bfd_database']
    mgnify = params['mgnify_database']

    hhblits_binary = params['hhblits_program']
    hhfilter_binary = params['hhfilter_program']
    jackhmmer_binary = params['jackhmmer_program']

    colabfold_search_binary = params['colabfold_search_program']
    colabfold_split_msas_binary = params['colabfold_split_msas_program']
    colabfold_databases = params['colabfold_databases']
    mmseq_binary = params['mmseq_program']

    if only_monomer:
        uniprot_fasta = ""
        uniprot_fasta_newest = ""
    else:
        uniprot_fasta = params['uniprot_fasta']
        uniprot_fasta_newest = params['uniprot_fasta_newest']

    result = None
    try:
        pipeline = Monomer_alignment_generation_pipeline(jackhmmer_binary_path=jackhmmer_binary,
                                                         hhblits_binary_path=hhblits_binary,
                                                         hhfilter_binary_path=hhfilter_binary,
                                                         colabfold_search_binary=colabfold_search_binary,
                                                         colabfold_split_msas_binary=colabfold_split_msas_binary,
                                                         mmseq_binary=mmseq_binary,
                                                         uniref90_database_path=uniref90_fasta,
                                                         uniref90_database_path_new=params['uniref90_fasta_newest'],
                                                         mgnify_database_path=mgnify,
                                                         small_bfd_database_path=smallbfd,
                                                         bfd_database_path=bfd,
                                                         uniref30_database_path=uniref30,
                                                         uniref30_database_path_new=uniref30_new,
                                                         uniclust30_database_path=uniclust30,
                                                         uniprot_database_path=uniprot_fasta,
                                                         uniprot_database_path_new=uniprot_fasta_newest,
                                                         colabfold_databases=colabfold_databases)
        result = pipeline.process(fasta, outdir)
    except Exception as e:
        print(e)
        return result
    return result


def copy_same_sequence_msas(srcdir, trgdir, srcname, trgname):
    for msa in os.listdir(srcdir):
        if msa.find('.a3m') > 0 or msa.find('.fasta') > 0:
            print(f"cp {srcdir}/{msa} {trgdir}/{msa}")
            os.system(f"cp {srcdir}/{msa} {trgdir}/{msa}")
            contents = open(f"{trgdir}/{msa}")
            new_contents = []
            for i, line in enumerate(contents):
                if i == 0:
                    new_contents += [f">{trgname}\n"]
                else:
                    new_contents += [line]
            fw = open(f"{trgdir}/{msa}", 'w')
            fw.writelines(new_contents)
        elif msa.find('.sto') > 0:
            os.system(f"cp {srcdir}/{msa} {trgdir}/{msa}")
            contents = []
            for line in open(f"{trgdir}/{msa}"):
                if line[:7] == "#=GF ID":
                    line = line.replace(srcname, trgname)

                tmp = line.split()
                if len(tmp) > 0 and tmp[0] == srcname:
                    line = line.replace(srcname, trgname)
                contents += [line]
            fw = open(f"{trgdir}/{msa}", 'w')
            fw.writelines(contents)
    cwd = os.getcwd()
    os.chdir(trgdir)
    print(f"rename {srcname} {trgname} *")
    os.system(f"rename {srcname} {trgname} *")
    os.chdir(cwd)


def run_monomer_msa_pipeline_img(fasta, outdir, params):
    monomer_msa_pipeline_img = Monomer_alignment_generation_pipeline_img(
        deepmsa_binary_path=params['deepmsa2_program'],
        bfd_database_path=params['bfd_database'],
        img_database_path=params['img_database'],
        metaclust_database_path=params['metaclust_database'],
        mgnify_database_path=params['mgnify_database'],
        uniref90_database_path=params['uniref90_fasta'])

    img_msa = None
    try:
        img_msa = monomer_msa_pipeline_img.process(fasta, outdir)
    except Exception as e:
        print(e)
    return img_msa


def run_monomer_template_search_pipeline(params, targetname, sequence, a3m, outdir):
    template_file = None
    try:
        pipeline = monomer_sequence_based_template_search_pipeline(params)
        template_file = pipeline.search(targetname=targetname, sequence=sequence, a3m=a3m, outdir=outdir)
    except Exception as e:
        print(e)
    return template_file


# def run_monomer_structure_generation_pipeline(params, run_methods, fasta_path, alndir, templatedir, outdir):
#     try:
#         pipeline = Monomer_structure_prediction_pipeline(params, run_methods=run_methods)
#         pipeline.process_single(fasta_path=fasta_path,
#                                 alndir=alndir,
#                                 template_dir=templatedir,
#                                 outdir=outdir)
#     except Exception as e:
#         print(e)
#         return False
#     return True


def run_monomer_structure_generation_pipeline_v2(params, fasta_path, alndir, templatedir, outdir, run_methods=None):
    try:
        pipeline = Monomer_structure_prediction_pipeline_v2(params, run_methods=run_methods)
        pipeline.process_single(fasta_path=fasta_path,
                                alndir=alndir,
                                template_dir=templatedir,
                                outdir=outdir)
    except Exception as e:
        print(e)
        return False
    return True


def select_models_monomer_only(qa_result, outputdir, params):
    if "pairwise_af_avg" not in qa_result:
        raise Exception(
            f"Cannot find pairwise ranking file for generating multicom-egnn models: {qa_result['pairwise_af_avg']}")

    selected_models = []
    pairwise_ranking_df = pd.read_csv(qa_result["pairwise_af_avg"])
    for i in range(2):
        model = pairwise_ranking_df.loc[i, 'model']
        os.system(f"cp {outputdir}/pdb/{model} {outputdir}/egnn{i + 1}.pdb")
        selected_models += [model]

    egnn_model_count = 3
    egnn_added_models = []
    top1_model = f"{outputdir}/pdb/{pairwise_ranking_df.loc[0, 'model']}"
    for i in range(2, len(pairwise_ranking_df)):
        model = pairwise_ranking_df.loc[i, 'model']
        tmscore, gdtscore = cal_tmscore(params['tmscore_program'],
                                        f"{outputdir}/pdb/{model}",
                                        top1_model,
                                        outputdir + '/tmp')
        if tmscore < 0.98:
            os.system(f"cp {outputdir}/pdb/{model} {outputdir}/egnn{egnn_model_count}.pdb")
            egnn_added_models += [model]
            selected_models += [model]
            egnn_model_count += 1
            if egnn_model_count > 5:
                break
        else:
            print(f"The tmscore between {model} and {top1_model} is larger than 0.98 ({tmscore}), skipped!")

    if egnn_model_count <= 5:
        for i in range(2, len(pairwise_ranking_df)):
            model = pairwise_ranking_df.loc[i, 'model']
            if model in egnn_added_models:
                continue
            os.system(f"cp {outputdir}/pdb/{model} {outputdir}/egnn{egnn_model_count}.pdb")
            egnn_added_models += [model]
            egnn_model_count += 1
            selected_models += [model]
            if egnn_model_count > 5:
                break
    selected_df = pd.DataFrame({'selected_models': selected_models})
    selected_df.to_csv(outputdir + '/egnn_selected.csv')

    selected_models = []
    alphafold_ranking_df = pd.read_csv(qa_result["alphafold"])
    for i in range(2):
        model = alphafold_ranking_df.loc[i, 'model']
        os.system(f"cp {outputdir}/pdb/{model} {outputdir}/deep{i + 1}.pdb")
        selected_models += [model]

    deep_model_count = 3
    deep_added_models = []
    top1_model = f"{outputdir}/pdb/{alphafold_ranking_df.loc[0, 'model']}"
    for i in range(2, len(alphafold_ranking_df)):
        model = alphafold_ranking_df.loc[i, 'model']
        tmscore, gdtscore = cal_tmscore(params['tmscore_program'],
                                        f"{outputdir}/pdb/{model}",
                                        top1_model,
                                        outputdir + '/tmp')
        if tmscore < 0.98:
            os.system(f"cp {outputdir}/pdb/{model} {outputdir}/deep{deep_model_count}.pdb")
            deep_model_count += 1
            deep_added_models += [model]
            selected_models += [model]
            if deep_model_count > 5:
                break
        else:
            print(f"The tmscore between {model} and {top1_model} is larger than 0.98 ({tmscore}), skipped!")

    if deep_model_count <= 5:
        for i in range(2, len(alphafold_ranking_df)):
            model = alphafold_ranking_df.loc[i, 'model']
            if model in deep_added_models:
                continue
            os.system(f"cp {outputdir}/pdb/{model} {outputdir}/deep{deep_model_count}.pdb")
            egnn_added_models += [model]
            deep_model_count += 1
            selected_models += [model]
            if deep_model_count > 5:
                break
    selected_df = pd.DataFrame({'selected_models': selected_models})
    selected_df.to_csv(outputdir + '/deep_selected.csv')


def select_models_with_multimer(qa_result, outputdir):
    if "pairwise_af_avg" not in qa_result:
        raise Exception(
            f"Cannot find pairwise ranking file for generating multicom-egnn models: {qa_result['pairwise_af_avg']}")

    # pdbs_from_monomer = [qa_result['pairwise_af_avg_monomer'].loc[i, 'model']
    #                      for i in range(len(qa_result['pairwise_af_avg_monomer']))]
    #
    # pdbs_from_multimer = [qa_result['pairwise_af_avg_multimer'].loc[i, 'model']
    #                       for i in range(len(qa_result['pairwise_af_avg_multimer']))]
    pairwise_af_avg_multimer_ranking = pd.read_csv(qa_result['pairwise_af_avg_multimer'])
    pairwise_af_avg_ranking = pd.read_csv(qa_result['pairwise_af_avg'])
    selected_multimer_models = [pairwise_af_avg_multimer_ranking.loc[i, 'model'] for i in range(3)]
    for i in range(len(pairwise_af_avg_ranking)):
        if pairwise_af_avg_ranking.loc[i, 'model'] in selected_multimer_models:
            continue
        else:
            selected_multimer_models += [pairwise_af_avg_ranking.loc[i, 'model']]
            if len(selected_multimer_models) >= 5:
                break

    for i in range(len(selected_multimer_models)):
        os.system(f"cp {outputdir}/pdb/{selected_multimer_models[i]} {outputdir}/egnn{i + 1}.pdb")

    selected_df = pd.DataFrame({'selected_models': selected_multimer_models})
    selected_df.to_csv(outputdir + '/egnn_selected.csv')

    alphafold_multimer_ranking = pd.read_csv(qa_result['alphafold_multimer'])
    alphafold_ranking = pd.read_csv(qa_result['alphafold'])
    selected_multimer_models = [alphafold_multimer_ranking.loc[i, 'model'] for i in range(3)]
    for i in range(len(alphafold_ranking)):
        if alphafold_ranking.loc[i, 'model'] in selected_multimer_models:
            continue
        else:
            selected_multimer_models += [alphafold_ranking.loc[i, 'model']]
            if len(selected_multimer_models) >= 5:
                break

    for i in range(len(selected_multimer_models)):
        os.system(f"cp {outputdir}/pdb/{selected_multimer_models[i]} {outputdir}/refine{i + 1}.pdb")

    selected_df = pd.DataFrame({'selected_models': selected_multimer_models})
    selected_df.to_csv(outputdir + '/refine_selected.csv')


def run_monomer_evaluation_pipeline(params, targetname, fasta_file, input_monomer_dir, outputdir, input_multimer_dir="",
                                    generate_egnn_models=False, model_count=5):
    makedir_if_not_exists(outputdir)
    qa_result = None
    pipeline = Monomer_structure_evaluation_pipeline(params=params,
                                                     use_gpu=True)
    try:
        qa_result = pipeline.process(targetname=targetname, fasta_file=fasta_file,
                                     monomer_model_dir=input_monomer_dir, multimer_model_dir=input_multimer_dir,
                                     output_dir=outputdir, model_count=model_count)
    except Exception as e:
        print(e)

    if generate_egnn_models:
        if input_multimer_dir == "" or not os.path.exists(input_multimer_dir):
            select_models_monomer_only(qa_result=qa_result, outputdir=outputdir, params=params)
        else:
            select_models_with_multimer(qa_result=qa_result, outputdir=outputdir)

    return qa_result


def select_models_monomer_only_human(qa_result, outputdir, params):
    if "pairwise_af_avg" not in qa_result:
        raise Exception(
            f"Cannot find pairwise ranking file for generating multicom models: {qa_result['pairwise_af_avg']}")

    selected_models = []
    pairwise_ranking_df = pd.read_csv(qa_result["pairwise_af_avg"])
    for i in range(4):
        model = pairwise_ranking_df.loc[i, 'model']
        os.system(f"cp {outputdir}/pdb/{model} {outputdir}/multicom{i + 1}.pdb")
        selected_models += [model]

    multicom_model_count = 5
    multicom_added_models = []
    top1_model = f"{outputdir}/pdb/{pairwise_ranking_df.loc[0, 'model']}"
    for i in range(4, len(pairwise_ranking_df)):
        model = pairwise_ranking_df.loc[i, 'model']
        tmscore, gdtscore = cal_tmscore(params['tmscore_program'],
                                        f"{outputdir}/pdb/{model}",
                                        top1_model,
                                        outputdir + '/tmp')
        if tmscore < 0.98:
            os.system(f"cp {outputdir}/pdb/{model} {outputdir}/multicom{multicom_model_count}.pdb")
            multicom_added_models += [model]
            selected_models += [model]
            multicom_model_count += 1
            if multicom_model_count > 5:
                break
        else:
            print(f"The tmscore between {model} and {top1_model} is larger than 0.98 ({tmscore}), skipped!")

    if multicom_model_count <= 5:
        for i in range(4, len(pairwise_ranking_df)):
            model = pairwise_ranking_df.loc[i, 'model']
            if model in multicom_added_models:
                continue
            os.system(f"cp {outputdir}/pdb/{model} {outputdir}/multicom{multicom_model_count}.pdb")
            multicom_added_models += [model]
            multicom_model_count += 1
            selected_models += [model]
            if multicom_model_count > 5:
                break
    selected_df = pd.DataFrame({'selected_models': selected_models})
    selected_df.to_csv(outputdir + '/multicom_selected.csv')

    selected_models = []
    casp13_human_ranking_df = read_qa_txt_as_df(qa_result["casp13"])
    casp13_human_ranking_df['model'] = casp13_human_ranking_df['model'] + '.pdb'
    for i in range(4):
        model = casp13_human_ranking_df.loc[i, 'model']
        os.system(f"cp {outputdir}/pdb/{model} {outputdir}/multicom_human{i + 1}.pdb")
        selected_models += [model]

    multicom_human_model_count = 5
    multicom_human_added_models = []
    top1_model = f"{outputdir}/pdb/{casp13_human_ranking_df.loc[0, 'model']}"
    for i in range(4, len(casp13_human_ranking_df)):
        model = casp13_human_ranking_df.loc[i, 'model']
        tmscore, gdtscore = cal_tmscore(params['tmscore_program'],
                                        f"{outputdir}/pdb/{model}",
                                        top1_model,
                                        outputdir + '/tmp')
        if tmscore < 0.98:
            os.system(f"cp {outputdir}/pdb/{model} {outputdir}/multicom_human{multicom_human_model_count}.pdb")
            multicom_human_added_models += [model]
            selected_models += [model]
            multicom_human_model_count += 1
            if multicom_human_model_count > 5:
                break
        else:
            print(f"The tmscore between {model} and {top1_model} is larger than 0.98 ({tmscore}), skipped!")

    if multicom_human_model_count <= 5:
        for i in range(4, len(casp13_human_ranking_df)):
            model = casp13_human_ranking_df.loc[i, 'model']
            if model in multicom_human_added_models:
                continue
            os.system(f"cp {outputdir}/pdb/{model} {outputdir}/multicom_human{multicom_human_model_count}.pdb")
            multicom_human_added_models += [model]
            multicom_human_model_count += 1
            selected_models += [model]
            if multicom_human_model_count > 5:
                break
    selected_df = pd.DataFrame({'selected_models': selected_models})
    selected_df.to_csv(outputdir + '/multicom_human_selected.csv')


# def select_models_with_multimer(qa_result, outputdir):
#     if "pairwise_af_avg" not in qa_result:
#         raise Exception(
#             f"Cannot find pairwise ranking file for generating multicom-egnn models: {qa_result['pairwise_af_avg']}")
#
#     # pdbs_from_monomer = [qa_result['pairwise_af_avg_monomer'].loc[i, 'model']
#     #                      for i in range(len(qa_result['pairwise_af_avg_monomer']))]
#     #
#     # pdbs_from_multimer = [qa_result['pairwise_af_avg_multimer'].loc[i, 'model']
#     #                       for i in range(len(qa_result['pairwise_af_avg_multimer']))]
#     pairwise_af_avg_multimer_ranking = pd.read_csv(qa_result['pairwise_af_avg_multimer'])
#     pairwise_af_avg_ranking = pd.read_csv(qa_result['pairwise_af_avg'])
#     selected_multimer_models = [pairwise_af_avg_multimer_ranking.loc[i, 'model'] for i in range(3)]
#     for i in range(len(pairwise_af_avg_ranking)):
#         if pairwise_af_avg_ranking.loc[i, 'model'] in selected_multimer_models:
#             continue
#         else:
#             selected_multimer_models += [pairwise_af_avg_ranking.loc[i, 'model']]
#             if len(selected_multimer_models) >= 5:
#                 break
#
#     for i in range(len(selected_multimer_models)):
#         os.system(f"cp {outputdir}/pdb/{selected_multimer_models[i]} {outputdir}/egnn{i + 1}.pdb")
#
#     selected_df = pd.DataFrame({'selected_models': selected_multimer_models})
#     selected_df.to_csv(outputdir + '/egnn_selected.csv')
#
#     alphafold_multimer_ranking = pd.read_csv(qa_result['alphafold_multimer'])
#     alphafold_ranking = pd.read_csv(qa_result['alphafold'])
#     selected_multimer_models = [alphafold_multimer_ranking.loc[i, 'model'] for i in range(3)]
#     for i in range(len(alphafold_ranking)):
#         if alphafold_ranking.loc[i, 'model'] in selected_multimer_models:
#             continue
#         else:
#             selected_multimer_models += [alphafold_ranking.loc[i, 'model']]
#             if len(selected_multimer_models) >= 5:
#                 break
#
#     for i in range(len(selected_multimer_models)):
#         os.system(f"cp {outputdir}/pdb/{selected_multimer_models[i]} {outputdir}/refine{i + 1}.pdb")
#
#     selected_df = pd.DataFrame({'selected_models': selected_multimer_models})
#     selected_df.to_csv(outputdir + '/refine_selected.csv')


def run_monomer_evaluation_pipeline_human(params, targetname, fasta_file, input_monomer_dir, outputdir,
                                          input_multimer_dir="", generate_egnn_models=False, model_count=5,
                                          dncon2_file="", dncon4_file="", distmap=""):
    makedir_if_not_exists(outputdir)
    qa_result = None
    pipeline = Monomer_structure_evaluation_human_pipeline(params=params,
                                                           use_gpu=True)
    try:
        qa_result = pipeline.process(targetname=targetname, fasta_file=fasta_file,
                                     monomer_model_dir=input_monomer_dir, multimer_model_dir=input_multimer_dir,
                                     output_dir=outputdir, model_count=model_count,
                                     dncon2_file=dncon2_file, dncon4_file=dncon4_file, distmap=distmap)
    except Exception as e:
        print(e)

    if generate_egnn_models:
        if input_multimer_dir == "" or not os.path.exists(input_multimer_dir):
            select_models_monomer_only_human(qa_result=qa_result, outputdir=outputdir, params=params)
        else:
            select_models_with_multimer(qa_result=qa_result, outputdir=outputdir)

    return qa_result


def rerun_monomer_evaluation_pipeline_human(params, targetname, fasta_file, outputdir,
                                            run_methods=None,
                                            dncon2_file="", dncon4_file="", distmap=""):
    makedir_if_not_exists(outputdir)
    result = None
    pipeline = Monomer_structure_evaluation_human_pipeline(params=params,
                                                           run_methods=run_methods,
                                                           use_gpu=True)
    try:
        result = pipeline.reprocess(targetname=targetname, fasta_file=fasta_file,
                                    output_dir=outputdir, dncon2_file=dncon2_file, dncon4_file=dncon4_file,
                                    distmap=distmap)
    except Exception as e:
        print(e)

    return result


def rerun_monomer_evaluation_pipeline(params, targetname, fasta_file, outputdir):
    makedir_if_not_exists(outputdir)
    result = None
    pipeline = Monomer_structure_evaluation_pipeline(params=params,
                                                     use_gpu=True)
    try:
        result = pipeline.reprocess(targetname=targetname, fasta_file=fasta_file,
                                    output_dir=outputdir)
    except Exception as e:
        print(e)

    return result


def run_monomer_refinement_pipeline(params, refinement_inputs, outdir, finaldir, prefix):
    pipeline = iterative_refine_pipeline.Monomer_iterative_refinement_pipeline_server(params=params)
    pipeline.search(refinement_inputs=refinement_inputs, outdir=outdir)

    makedir_if_not_exists(finaldir)

    pipeline = iterative_refine_pipeline.Monomer_refinement_model_selection(params)
    pipeline.select_v1(indir=outdir, outdir=finaldir, prefix=prefix)
    # pipeline.select_v2(indir=outdir, outdir=finaldir + '/v2', ranking_df=ranking_df)


def run_monomer_refinement_pipeline_human(params, ranking_df, refinement_inputs, outdir, finaldir, prefix, refpdb=""):
    pipeline = iterative_refine_pipeline.Monomer_iterative_refinement_pipeline_server(params=params)
    pipeline.search(refinement_inputs=refinement_inputs, outdir=outdir)

    makedir_if_not_exists(finaldir)

    pipeline = iterative_refine_pipeline.Monomer_refinement_model_selection(params)
    pipeline.select_v2(indir=outdir, outdir=finaldir, ranking_df=ranking_df, prefix=prefix, refpdb=refpdb)


def run_concatenate_dimer_msas_pipeline(multimer, run_methods, monomer_aln_dir, outputdir, params, is_homomers=False):
    chains = multimer.split(',')
    # alignment = {'outdir': f"{outputdir}/{'_'.join(chains)}"}
    alignment = {'outdir': outputdir}
    for i in range(len(chains)):
        chain = chains[i]
        chain_aln_dir = monomer_aln_dir + '/' + chain
        if os.path.exists(chain_aln_dir):
            chain_a3ms = {'name': chain,
                          'colabfold_a3m': f"{chain_aln_dir}/{chain}_colabfold.a3m",
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

        if not os.path.exists(alignment['outdir'] + '/DONE'):
            complex_alignment_concatenation_pipeline = Complex_alignment_concatenation_pipeline(params=params,
                                                                                                run_methods=run_methods)
            alignments = complex_alignment_concatenation_pipeline.concatenate(alignments, params['hhfilter_program'],
                                                                              is_homomers=is_homomers)
        else:
            print("The multimer alignments have been generated!")
    else:
        print("The a3ms for dimers are not complete!")


def run_complex_template_search_pipeline(multimers, monomer_aln_dir, monomer_model_dir, outdir, params):
    monomer_template_inputs = []
    monomer_sequences = []
    for chain in multimers:
        chain_template_a3m = f"{monomer_aln_dir}/{chain}/{chain}_uniref90.sto"
        if not os.path.exists(chain_template_a3m):
            raise Exception(f"Cannot find uniref90 alignment for {chain}")
        seq = open(f"{monomer_aln_dir}/{chain}/{chain}.fasta").readlines()[1].rstrip('\n')
        chain_template_input = sequence_based_pipeline.monomer_template_input(name=chain,
                                                                              msa_path=chain_template_a3m,
                                                                              hmm_path="", seq=seq)
        monomer_template_inputs += [chain_template_input]
        monomer_sequences += [seq]

    print("searching complex sequence based template search pipelineL RCSB_PDB")
    pdb_seq_dir = outdir + '/pdb_seq'
    makedir_if_not_exists(pdb_seq_dir)
    if not os.path.exists(pdb_seq_dir + '/sequence_templates.csv'):
        pipeline = sequence_based_pipeline_pdb.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, pdb_seq_dir)

    print("searching complex sequence based template search pipeline: Complex")
    complex_pdb_seq_dir = outdir + '/complex_pdb_seq'
    makedir_if_not_exists(complex_pdb_seq_dir)
    if not os.path.exists(complex_pdb_seq_dir + '/sequence_templates.csv'):
        pipeline = sequence_based_pipeline_complex_pdb.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, complex_pdb_seq_dir)

    print("searching complex sequence based template search pipeline: pdb70")
    pdb70_seq_dir = outdir + '/pdb70_seq'
    makedir_if_not_exists(pdb70_seq_dir)
    if not os.path.exists(pdb70_seq_dir + '/sequence_templates.csv'):
        pipeline = sequence_based_pipeline.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, pdb70_seq_dir)

    print("searching complex structure based template search pipeline")
    struct_temp_dir = outdir + '/struct_temp'
    makedir_if_not_exists(struct_temp_dir)
    if not os.path.exists(struct_temp_dir + '/structure_templates.csv'):
        monomer_pdbs = []
        for chain in multimers:
            monomer_pdb = f"{monomer_model_dir}/{chain}/default/ranked_0.pdb"
            if not os.path.exists(monomer_pdb):
                print(f"Cannot find teritary structure for {chain}: {monomer_pdb}")
                continue
            os.system(f"cp {monomer_pdb} {struct_temp_dir}/{chain}.pdb")
            monomer_pdbs += [f"{struct_temp_dir}/{chain}.pdb"]
        pipeline = structure_based_pipeline_v2.Complex_structure_based_template_search_pipeline(params)
        pipeline.search(monomer_sequences, monomer_pdbs, struct_temp_dir)


# def run_quaternary_structure_generation_pipeline(params, fasta_path, chain_id_map, aln_dir, complex_aln_dir,
#                                                  template_dir,
#                                                  monomer_model_dir, output_dir):
#     try:
#         pipeline = Quaternary_structure_prediction_pipeline(params)
#         result = pipeline.process(fasta_path=fasta_path,
#                                   chain_id_map=chain_id_map,
#                                   aln_dir=aln_dir,
#                                   complex_aln_dir=complex_aln_dir,
#                                   template_dir=template_dir,
#                                   monomer_model_dir=monomer_model_dir,
#                                   output_dir=output_dir)
#     except Exception as e:
#         print(e)
#         return False
#     return True


def run_quaternary_structure_generation_pipeline_v2(params, fasta_path, chain_id_map, aln_dir, complex_aln_dir,
                                                    template_dir,
                                                    monomer_model_dir, output_dir, run_methods=None):
    try:
        pipeline = Quaternary_structure_prediction_pipeline_v2(params, run_methods)
        result = pipeline.process(fasta_path=fasta_path,
                                  chain_id_map=chain_id_map,
                                  aln_dir=aln_dir,
                                  complex_aln_dir=complex_aln_dir,
                                  template_dir=template_dir,
                                  monomer_model_dir=monomer_model_dir,
                                  output_dir=output_dir)
    except Exception as e:
        print(e)
        return False
    return True


def run_quaternary_structure_generation_pipeline_default(params, fasta_path, chain_id_map, aln_dir, output_dir):
    try:
        pipeline = Quaternary_structure_prediction_pipeline_default(params)
        result = pipeline.process(fasta_path=fasta_path,
                                  chain_id_map=chain_id_map,
                                  aln_dir=aln_dir,
                                  output_dir=output_dir)
    except Exception as e:
        print(e)
        return False
    return True


# def run_quaternary_structure_generation_homo_pipeline(params, fasta_path, chain_id_map, aln_dir, complex_aln_dir,
#                                                       template_dir,
#                                                       monomer_model_dir, output_dir):
#     try:
#         pipeline = Quaternary_structure_prediction_homo_pipeline(params)
#         result = pipeline.process(fasta_path=fasta_path,
#                                   chain_id_map=chain_id_map,
#                                   aln_dir=aln_dir,
#                                   complex_aln_dir=complex_aln_dir,
#                                   template_dir=template_dir,
#                                   monomer_model_dir=monomer_model_dir,
#                                   output_dir=output_dir)
#     except Exception as e:
#         print(e)
#         return False
#     return True
#
#
def run_quaternary_structure_generation_homo_pipeline_img_v2(params, fasta_path, chain_id_map, aln_dir, output_dir):
    try:
        pipeline = Quaternary_structure_prediction_homo_pipeline_v2(params)
        result = pipeline.process_img(fasta_path=fasta_path,
                                      chain_id_map=chain_id_map,
                                      aln_dir=aln_dir,
                                      output_dir=output_dir)
    except Exception as e:
        print(e)
        return False
    return True


def run_quaternary_structure_generation_homo_pipeline_v2(params, fasta_path, chain_id_map, aln_dir, complex_aln_dir,
                                                         template_dir,
                                                         monomer_model_dir, output_dir, run_methods=None):
    try:
        pipeline = Quaternary_structure_prediction_homo_pipeline_v2(params, run_methods)
        result = pipeline.process(fasta_path=fasta_path,
                                  chain_id_map=chain_id_map,
                                  aln_dir=aln_dir,
                                  complex_aln_dir=complex_aln_dir,
                                  template_dir=template_dir,
                                  monomer_model_dir=monomer_model_dir,
                                  output_dir=output_dir)
    except Exception as e:
        print(e)
        return False
    return True


class foldseek_iterative_monomer_input:
    def __init__(self, monomer_pdb_dirs, monomer_alphafold_a3ms):
        self.monomer_pdb_dirs = monomer_pdb_dirs
        self.monomer_alphafold_a3ms = monomer_alphafold_a3ms


def run_quaternary_structure_generation_pipeline_foldseek(params, fasta_path, chain_id_map, pipeline_inputs, outdir,
                                                          start=0, is_homomers=False):
    pipeline = Multimer_iterative_generation_pipeline_monomer(params)
    try:
        for i, pipeline_input in enumerate(pipeline_inputs):
            if is_homomers:
                pipeline.search_single_homo(fasta_file=fasta_path,
                                            chain_id_map=chain_id_map,
                                            monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                            monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                            outdir=f"{outdir}/iter_{i + 1 + start}")
            else:
                pipeline.search_single(fasta_file=fasta_path,
                                       chain_id_map=chain_id_map,
                                       monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                       monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                       outdir=f"{outdir}/iter_{i + 1 + start}")
    except Exception as e:
        print(e)
        return False
    return True


def run_quaternary_structure_generation_pipeline_foldseek_old(params, fasta_path, chain_id_map, pipeline_inputs, outdir,
                                                              start=0, is_homomers=False):
    pipeline = Multimer_iterative_generation_pipeline_monomer_old(params)
    try:
        for i, pipeline_input in enumerate(pipeline_inputs):
            if is_homomers:
                pipeline.search_single_homo(fasta_file=fasta_path,
                                            chain_id_map=chain_id_map,
                                            monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                            monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                            outdir=f"{outdir}/iter_{i + 1 + start}_old")
            else:
                pipeline.search_single(fasta_file=fasta_path,
                                       chain_id_map=chain_id_map,
                                       monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                       monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                       outdir=f"{outdir}/iter_{i + 1 + start}_old")
    except Exception as e:
        print(e)
        return False
    return True


def extract_monomer_models_from_complex(complex_pdb, complex_pkl, chain_id_map, workdir):
    makedir_if_not_exists(workdir)

    chain_group = {}
    for chain_id in chain_id_map:
        find = False
        for chain_id_seq in chain_group:
            if chain_id_map[chain_id_seq].sequence == chain_id_map[chain_id].sequence:
                chain_group[chain_id_seq] += [chain_id]
                find = True
        if not find:
            chain_group[chain_id] = [chain_id]

    pkldir = workdir + '/pkl'
    makedir_if_not_exists(pkldir)

    alphafold_qa = Alphafold_pkl_qa(ranking_methods=['plddt_avg'])
    for chain_id in chain_group:

        chain_pdb_dict = extract_monomer_pdbs(complex_pdb=complex_pdb,
                                              sequence=chain_id_map[chain_id].sequence,
                                              output_prefix=workdir + '/chain')
        print(chain_pdb_dict)

        same_seq_chains = chain_group[chain_id]
        same_seq_pkldir = pkldir + '/' + chain_id
        makedir_if_not_exists(same_seq_pkldir)
        for same_seq_chain in same_seq_chains:
            pdbname = chain_pdb_dict[same_seq_chain]['pdbname']
            extract_pkl(src_pkl=complex_pkl,
                        residue_start=chain_pdb_dict[same_seq_chain]['chain_start'],
                        residue_end=chain_pdb_dict[same_seq_chain]['chain_end'],
                        output_pkl=same_seq_pkldir + '/' + pdbname.replace('.pdb', '.pkl'))
        alphafold_ranking = alphafold_qa.run(same_seq_pkldir)
        alphafold_ranking.to_csv(same_seq_pkldir + '_alphafold_ranking.csv')
        os.system(f"cp {workdir}/{alphafold_ranking.loc[0, 'model']} {workdir}/{chain_id}_top1.pdb")

    return chain_group


def rerun_multimer_evaluation_pipeline(params, fasta_path, chain_id_map, monomer_model_dir,
                                       outdir, stoichiometry):
    makedir_if_not_exists(outdir)
    pipeline = Quaternary_structure_evaluation_pipeline(params=params)
    multimer_qa_result = None
    try:
        multimer_qa_result = pipeline.process(fasta_path=fasta_path,
                                              chain_id_map=chain_id_map, monomer_model_dir=monomer_model_dir,
                                              output_dir=outdir, stoichiometry=stoichiometry)
    except Exception as e:
        print(e)


def run_multimer_evaluation_pipeline(params, fasta_path, chain_id_map, monomer_model_dir,
                                     indir, outdir, stoichiometry, is_homomer=False, model_count=5):
    makedir_if_not_exists(outdir)
    pipeline = Quaternary_structure_evaluation_pipeline(params=params)
    multimer_qa_result = None
    try:
        multimer_qa_result = pipeline.reprocess(fasta_path=fasta_path,
                                                chain_id_map=chain_id_map, monomer_model_dir=monomer_model_dir,
                                                model_dir=indir,
                                                output_dir=outdir, stoichiometry=stoichiometry, model_count=model_count)
    except Exception as e:
        print(e)

    alphafold_confidence_ranking = pd.read_csv(multimer_qa_result['alphafold'])
    for i in range(5):
        model_name = alphafold_confidence_ranking.loc[i, 'model']
        os.system(f"cp {outdir}/pdb/{model_name} {outdir}/qa{i + 1}.pdb")
        if not is_homomer:
            chain_group = extract_monomer_models_from_complex(complex_pdb=f"{outdir}/qa{i + 1}.pdb",
                                                              complex_pkl=f"{outdir}/pkl/{model_name.replace('.pdb', '.pkl')}",
                                                              chain_id_map=chain_id_map, workdir=f"{outdir}/qa{i + 1}")
            for chain_id in chain_group:
                chain_outdir = f"{outdir}/qa_{chain_id}"
                makedir_if_not_exists(chain_outdir)
                os.system(f"cp {outdir}/qa{i + 1}/{chain_id}_top1.pdb {chain_outdir}/qa{i + 1}.pdb")

    af_pairwise_avg_ranking = pd.read_csv(multimer_qa_result['pairwise_af_avg'])
    for i in range(5):
        model_name = af_pairwise_avg_ranking.loc[i, 'model']
        os.system(f"cp {outdir}/pdb/{model_name} {outdir}/deep{i + 1}.pdb")
        if not is_homomer:
            chain_group = extract_monomer_models_from_complex(complex_pdb=f"{outdir}/deep{i + 1}.pdb",
                                                              complex_pkl=f"{outdir}/pkl/{model_name.replace('.pdb', '.pkl')}",
                                                              chain_id_map=chain_id_map,
                                                              workdir=f"{outdir}/deep{i + 1}")
            for chain_id in chain_group:
                chain_outdir = f"{outdir}/deep_{chain_id}"
                makedir_if_not_exists(chain_outdir)
                os.system(f"cp {outdir}/deep{i + 1}/{chain_id}_top1.pdb {chain_outdir}/deep{i + 1}.pdb")

    return multimer_qa_result


def run_multimer_evaluation_pipeline_human(params, fasta_path, chain_id_map, monomer_model_dir,
                                           indir, extract_dir, outdir, stoichiometry, is_homomer=False, model_count=5):
    makedir_if_not_exists(outdir)
    pipeline = Quaternary_structure_evaluation_pipeline_human(params=params)
    multimer_qa_result = None
    try:
        multimer_qa_result = pipeline.process(fasta_path=fasta_path,
                                              chain_id_map=chain_id_map, monomer_model_dir=monomer_model_dir,
                                              model_dir=indir,
                                              extract_model_dir=extract_dir,
                                              output_dir=outdir, stoichiometry=stoichiometry, model_count=model_count)
    except Exception as e:
        print(e)

    alphafold_confidence_ranking = pd.read_csv(multimer_qa_result['alphafold'])
    for i in range(5):
        model_name = alphafold_confidence_ranking.loc[i, 'model']
        os.system(f"cp {outdir}/pdb/{model_name} {outdir}/qa{i + 1}.pdb")
        if not is_homomer:
            chain_group = extract_monomer_models_from_complex(complex_pdb=f"{outdir}/qa{i + 1}.pdb",
                                                              complex_pkl=f"{outdir}/pkl/{model_name.replace('.pdb', '.pkl')}",
                                                              chain_id_map=chain_id_map, workdir=f"{outdir}/qa{i + 1}")
            for chain_id in chain_group:
                chain_outdir = f"{outdir}/qa_{chain_id}"
                makedir_if_not_exists(chain_outdir)
                os.system(f"cp {outdir}/qa{i + 1}/{chain_id}_top1.pdb {chain_outdir}/qa{i + 1}.pdb")

    af_pairwise_avg_ranking = pd.read_csv(multimer_qa_result['pairwise_af_avg'])
    for i in range(5):
        model_name = af_pairwise_avg_ranking.loc[i, 'model']
        os.system(f"cp {outdir}/pdb/{model_name} {outdir}/deep{i + 1}.pdb")
        if not is_homomer:
            chain_group = extract_monomer_models_from_complex(complex_pdb=f"{outdir}/deep{i + 1}.pdb",
                                                              complex_pkl=f"{outdir}/pkl/{model_name.replace('.pdb', '.pkl')}",
                                                              chain_id_map=chain_id_map,
                                                              workdir=f"{outdir}/deep{i + 1}")
            for chain_id in chain_group:
                chain_outdir = f"{outdir}/deep_{chain_id}"
                makedir_if_not_exists(chain_outdir)
                os.system(f"cp {outdir}/deep{i + 1}/{chain_id}_top1.pdb {chain_outdir}/deep{i + 1}.pdb")

    return multimer_qa_result


def rerun_multimer_evaluation_pipeline_human(params, fasta_path, chain_id_map, monomer_model_dir,
                                             outdir, stoichiometry):
    makedir_if_not_exists(outdir)
    pipeline = Quaternary_structure_evaluation_pipeline_human(params=params)
    multimer_qa_result = None
    try:
        multimer_qa_result = pipeline.reprocess(fasta_path=fasta_path,
                                                chain_id_map=chain_id_map, monomer_model_dir=monomer_model_dir,
                                                output_dir=outdir, stoichiometry=stoichiometry)
    except Exception as e:
        print(e)

    return multimer_qa_result


def run_multimer_refinement_pipeline(params, chain_id_map, refinement_inputs, outdir, finaldir, stoichiometry):
    pipeline = iterative_refine_pipeline_multimer.Multimer_iterative_refinement_pipeline_server(params=params)
    pipeline.search(refinement_inputs=refinement_inputs, outdir=outdir, stoichiometry=stoichiometry)

    makedir_if_not_exists(finaldir)

    pipeline = iterative_refine_pipeline_multimer.Multimer_refinement_model_selection()
    pipeline.select_v1(indir=outdir, outdir=finaldir)

    for i in range(5):
        chain_group = extract_monomer_models_from_complex(complex_pdb=f"{finaldir}/deep{i + 1}.pdb",
                                                          complex_pkl=f"{finaldir}/deep{i + 1}.pkl",
                                                          chain_id_map=chain_id_map, workdir=f"{finaldir}/deep{i + 1}")
        for chain_id in chain_group:
            chain_outdir = f"{finaldir}/deep_{chain_id}"
            makedir_if_not_exists(chain_outdir)
            os.system(f"cp {finaldir}/deep{i + 1}/{chain_id}_top1.pdb {chain_outdir}/deep{i + 1}.pdb")
