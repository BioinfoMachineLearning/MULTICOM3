import os, sys, argparse, time, copy
from multiprocessing import Pool
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file
from bml_casp15.monomer_alignment_generation.alignment import read_fasta, write_fasta
from bml_casp15.monomer_alignment_generation.pipeline import *
from bml_casp15.monomer_structure_generation.pipeline import *
from bml_casp15.monomer_structure_evaluation.pipeline import *
from bml_casp15.monomer_templates_search.sequence_based_pipeline_pdb import *
from bml_casp15.monomer_structure_refinement import iterative_refine_pipeline
from bml_casp15.quaternary_structure_refinement import iterative_refine_pipeline_multimer
from bml_casp15.complex_alignment_generation.pipeline_v2 import *
from bml_casp15.complex_templates_search import sequence_based_pipeline_complex_pdb, \
    sequence_based_pipeline_pdb, sequence_based_pipeline, structure_based_pipeline_v2
from bml_casp15.quaternary_structure_generation.pipeline import *
from bml_casp15.quaternary_structure_generation.iterative_search_pipeline_v0_2 import *
from bml_casp15.quaternary_structure_evaluation.pipeline import *


def run_monomer_msa_pipeline(fasta, outdir, params):
    uniref30 = params['uniref_db_dir'] + '/' + params['uniref_db']
    uniclust30 = params['uniclust_db_dir'] + '/' + params['uniclust_db']
    uniref90_fasta = params['uniref90_fasta']
    uniprot_fasta = params['uniprot_fasta']
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

    result = None
    try:
        pipeline = Monomer_alignment_generation_pipeline(jackhmmer_binary_path=jackhmmer_binary,
                                                         hhblits_binary_path=hhblits_binary,
                                                         hhfilter_binary_path=hhfilter_binary,
                                                         colabfold_search_binary=colabfold_search_binary,
                                                         colabfold_split_msas_binary=colabfold_split_msas_binary,
                                                         mmseq_binary=mmseq_binary,
                                                         uniref90_database_path=uniref90_fasta,
                                                         mgnify_database_path=mgnify,
                                                         small_bfd_database_path=smallbfd,
                                                         bfd_database_path=bfd,
                                                         uniref30_database_path=uniref30,
                                                         uniclust30_database_path=uniclust30,
                                                         uniprot_database_path=uniprot_fasta,
                                                         colabfold_databases=colabfold_databases)
        result = pipeline.process(fasta, outdir)
    except Exception as e:
        print(e)
        return result
    return result


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


def run_monomer_structure_generation_pipeline(params, run_methods, fasta_path, alndir, templatedir, outdir):
    try:
        pipeline = Monomer_structure_prediction_pipeline(params, run_methods=run_methods)
        pipeline.process_single(fasta_path=fasta_path,
                                alndir=alndir,
                                template_dir=templatedir,
                                outdir=outdir)
    except Exception as e:
        print(e)
        return False
    return True


def run_monomer_evaluation_pipeline(params, targetname, fasta_file, input_monomer_dir, outputdir,
                                    chainid="", unrelaxed_chainid="", input_multimer_dir=""):
    makedir_if_not_exists(outputdir)
    result = None
    pipeline = Monomer_structure_evaluation_pipeline(params=params,
                                                     run_methods=["apollo", "alphafold", "enQA"],
                                                     use_gpu=True)
    try:
        result = pipeline.process(targetname=targetname, fasta_file=fasta_file,
                                  monomer_model_dir=input_monomer_dir, multimer_model_dir=input_multimer_dir,
                                  output_dir=outputdir, chainid_in_multimer=chainid,
                                  unrelaxed_chainid_in_multimer=unrelaxed_chainid)
    except Exception as e:
        print(e)
    return result


def rerun_monomer_evaluation_pipeline(params, targetname, fasta_file, outputdir):
    if not os.path.exists(outputdir + '/pdb'):
        raise Exception(f"cannot find the pdb directory in {outputdir}")
    result = None
    pipeline = Monomer_structure_evaluation_pipeline(params=params,
                                                     run_methods=["apollo", "alphafold", "enQA"],
                                                     use_gpu=True)
    try:
        result = pipeline.reprocess(targetname=targetname, fasta_file=fasta_file, output_dir=outputdir)
    except Exception as e:
        print(e)
    return result


def run_monomer_refinement_pipeline(params, refinement_inputs, outdir, finaldir):
    pipeline = iterative_refine_pipeline.Monomer_iterative_refinement_pipeline_server(params=params)
    pipeline.search(refinement_inputs=refinement_inputs, outdir=outdir)

    makedir_if_not_exists(finaldir)

    pipeline = iterative_refine_pipeline.Monomer_refinement_model_selection()
    pipeline.select_v1(indir=outdir, outdir=finaldir)


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
        if not os.path.exists(chain_template_a3m):
            raise Exception(f"Cannot find uniref90 alignment for {chain}")
        seq = open(f"{monomer_aln_dir}/{chain}/{chain}.fasta").readlines()[1].rstrip('\n')
        chain_template_input = sequence_based_pipeline.monomer_template_input(name=chain,
                                                                              msa_path=chain_template_a3m,
                                                                              hmm_path="", seq=seq)
        monomer_template_inputs += [chain_template_input]

    pdb_seq_dir = outdir + '/pdb_seq'
    makedir_if_not_exists(pdb_seq_dir)
    if not os.path.exists(pdb_seq_dir + '/sequence_templates.csv'):
        pipeline = sequence_based_pipeline_pdb.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, pdb_seq_dir)

    complex_pdb_seq_dir = outdir + '/complex_pdb_seq'
    makedir_if_not_exists(complex_pdb_seq_dir)
    if not os.path.exists(complex_pdb_seq_dir + '/sequence_templates.csv'):
        pipeline = sequence_based_pipeline_complex_pdb.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, complex_pdb_seq_dir)

    pdb70_seq_dir = outdir + '/pdb70_seq'
    makedir_if_not_exists(pdb70_seq_dir)
    if not os.path.exists(pdb70_seq_dir + '/sequence_templates.csv'):
        pipeline = sequence_based_pipeline.Complex_sequence_based_template_search_pipeline(params)
        pipeline.search(monomer_template_inputs, pdb70_seq_dir)

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
        pipeline.search(monomer_pdbs, struct_temp_dir)


def run_quaternary_structure_generation_pipeline(params, fasta_path, chain_id_map, aln_dir, complex_aln_dir,
                                                 template_dir,
                                                 monomer_model_dir, output_dir):
    try:
        pipeline = Quaternary_structure_prediction_pipeline(params)
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
                                                          start=0):
    pipeline = Multimer_iterative_generation_pipeline_monomer(params)
    try:
        for i, pipeline_input in enumerate(pipeline_inputs):
            pipeline.search_single(fasta_file=fasta_path,
                                   chain_id_map=chain_id_map,
                                   monomer_pdb_dirs=pipeline_input.monomer_pdb_dirs,
                                   monomer_alphafold_a3ms=pipeline_input.monomer_alphafold_a3ms,
                                   outdir=f"{outdir}/iter_{i + 1 + start}")
    except Exception as e:
        print(e)
        return False
    return True


def run_multimer_evaluation_pipeline(params, chain_id_map, indir, outdir):
    makedir_if_not_exists(outdir)
    pipeline = Quaternary_structure_evaluation_pipeline(params=params)
    multimer_qa_result = None
    try:
        multimer_qa_result = pipeline.process(chain_id_map=chain_id_map, model_dir=indir, output_dir=outdir)
    except Exception as e:
        print(e)
    return multimer_qa_result


def run_multimer_refinement_pipeline(params, refinement_inputs, outdir, finaldir):
    pipeline = iterative_refine_pipeline_multimer.Multimer_iterative_refinement_pipeline_server(params=params)
    pipeline.search(refinement_inputs=refinement_inputs, outdir=outdir)

    makedir_if_not_exists(finaldir)

    pipeline = iterative_refine_pipeline_multimer.Multimer_refinement_model_selection()
    pipeline.select_v1(indir=outdir, outdir=finaldir)
