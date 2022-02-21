import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import pathlib


def complete_result(outputdir):
    complete = True
    for i in range(0, 5):
        model = f'{outputdir}/ranked_{i}.pdb'
        if not os.path.exists(model):
            complete = False
            break
    return complete


def get_complex_alignments_by_method(monomers, concatenate_method, aln_dir):
    a3ms_path = []
    for monomer in monomers:
        if concatenate_method == 'uniclust_oxmatch_a3m':
            monomer_a3m = f"{aln_dir}/{monomer}/{monomer}_uniclust30.a3m"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniref_a3m') > 0:
            monomer_a3m = f"{aln_dir}/{monomer}/{monomer}_uniref30.a3m"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniref_sto') > 0:
            monomer_a3m = f"{aln_dir}/{monomer}/{monomer}_uniref90.sto"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniprot_sto') > 0:
            monomer_a3m = f"{aln_dir}/{monomer}/{monomer}_uniprot.sto"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
    return a3ms_path


class Quaternary_structure_prediction_pipeline:

    def __init__(self, params):  # , run_methods):

        self.params = params

        # self.run_methods = run_methods

        self.concatenate_methods = ['uniclust_oxmatch_a3m',
                                    'pdb_interact_uniref_a3m',
                                    'species_interact_uniref_a3m',
                                    'uniprot_distance_uniref_a3m',
                                    'string_interact_uniref_a3m',
                                    # 'geno_dist_uniref_a3m',
                                    'pdb_interact_uniref_sto',
                                    'species_interact_uniref_sto',
                                    'uniprot_distance_uniref_sto',
                                    'string_interact_uniref_sto',
                                    # 'geno_dist_uniref_sto',
                                    'pdb_interact_uniprot_sto',
                                    'species_interact_uniprot_sto',
                                    'uniprot_distance_uniprot_sto',
                                    'string_interact_uniprot_sto']

        self.template_methods = ['original_template_pipeline',
                                 'sequence_based_template',
                                 'structure_based_template',
                                 'structure_based_template_alphafold']

    def process(self,
                fasta_path,
                aln_dir,
                complex_aln_dir,
                structure_template_dir,
                sequence_template_dir,
                monomer_model_dir,
                output_dir):

        makedir_if_not_exists(output_dir)

        os.chdir(self.params['alphafold_program_dir'])

        monomers = []

        for line in open(fasta_path):
            line = line.rstrip('\n').strip()
            if line.startswith('>'):
                monomers += [line[1:]]

        # Customized complex alignment pipelines using original template search pipeline in alphafold
        template_stos = []
        for monomer in monomers:
            monomer_template_sto = f"{aln_dir}/{monomer}/{monomer}_uniref90.sto"
            if not os.path.exists(monomer_template_sto):
                raise Exception(f"Cannot find template stos for {monomer}: {monomer_template_sto}")
            template_stos += [monomer_template_sto]

        for concatenate_method in self.concatenate_methods:
            msa_pair_file = f"{complex_aln_dir}/{concatenate_method}/{concatenate_method}_interact.csv"
            if len(pd.read_csv(msa_pair_file)) == 0:
                continue

            outdir = f"{output_dir}/{concatenate_method}"

            if complete_result(outdir):
                continue

            makedir_if_not_exists(outdir)

            a3m_paths = get_complex_alignments_by_method(monomers=monomers,
                                                         concatenate_method=concatenate_method,
                                                         aln_dir=aln_dir)

            cmd = f"python {self.params['alphafold_multimer_program']} " \
                  f"--fasta_path {fasta_path} " \
                  f"--a3ms {','.join(a3m_paths)} " \
                  f"--msa_pair_file {msa_pair_file} " \
                  f"--template_stos {','.join(template_stos)} " \
                  f"--env_dir {self.params['alphafold_env_dir']} " \
                  f"--database_dir {self.params['alphafold_database_dir']} " \
                  f"--output_dir {outdir}"

            if complete_result(outdir):
                continue

            print(cmd)
            os.system(cmd)

        # Customized template search pipelines

        a3m_paths = []
        for monomer in monomers:
            monomer_a3m = f"{aln_dir}/{monomer}/{monomer}_uniprot.sto"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3m_paths += [monomer_a3m]

        for template_method in self.template_methods:

            outdir = f"{output_dir}/{template_method}"

            if complete_result(outdir):
                continue

            makedir_if_not_exists(outdir)

            base_cmd = f"python run_alphafold_multimer_custom_sim.py  " \
                       f"--fasta_path {fasta_path} " \
                       f"--a3ms {','.join(a3m_paths)} " \
                       f"--env_dir {self.params['alphafold_env_dir']} " \
                       f"--database_dir {self.params['alphafold_database_dir']} " \
                       f"--output_dir {outdir} "

            if template_method == "original_template_pipeline":

                base_cmd += f"--template_stos {','.join(template_stos)} "

            elif template_method == "structure_based_template":

                template_file = f"{structure_template_dir}/structure_templates.csv"

                base_cmd += f"--temp_struct_csv {template_file} "

            elif template_method == "structure_based_template_alphafold":

                template_file = f"{structure_template_dir}/structure_templates.csv"

                monomer_paths = []
                for monomer in monomers:
                    monomer_path = f"{monomer_model_dir}/{monomer}"
                    if not os.path.exists(monomer_path):
                        raise Exception(f"Cannot find monomer directory for {monomer}: {monomer_path}")
                    monomer_paths += [monomer_path]

                base_cmd += f"--temp_struct_csv {template_file} "
                base_cmd += f"--monomer_paths {','.join(monomer_paths)} "

            elif template_method == "sequence_based_template":

                template_file = f"{sequence_template_dir}/sequence_templates.csv"

                template_hits_files = []
                for monomer in monomers:
                    template_hits_file = f"{sequence_template_dir}/{monomer}/pdb_hits.hhr"
                    if not os.path.exists(template_hits_file):
                        raise Exception(f"Cannot find template hit file for {monomer}: {template_hits_file}")
                    template_hits_files += [template_hits_file]

                base_cmd += f"--temp_seq_pair_file {template_file} "
                base_cmd += f"--template_hits_files {','.join(template_hits_files)} "

            if len(base_cmd) > 0:
                print(base_cmd)
                os.system(base_cmd)

        print("The structure based template searching for dimers has finished!")