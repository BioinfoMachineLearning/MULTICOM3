import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
from bml_casp15.common.protein import complete_result, parse_fasta
import pandas as pd
from multiprocessing import Pool
import pathlib


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


class Quaternary_structure_prediction_pipeline_v2:

    def __init__(self, params, run_methods = None):

        self.params = params

        if run_methods is None:
            self.run_methods = ['default',
                                'default_mul_newest',
                                'default+sequence_based_template_pdb70',
                                'default_uniclust30',
                                'default_uniref30_22',
                                'default+structure_based_template',
                                'default+sequence_based_template_pdb',
                                'default+sequence_based_template_complex_pdb',
                                'default+alphafold_model_templates',
                                'uniclust_oxmatch_a3m',
                                'pdb_interact_uniref_a3m',
                                'species_interact_uniref_a3m',
                                'species_interact_uniref_a3m+sequence_based_template_pdb70',
                                'species_interact_uniref_a3m+structure_based_template',
                                'species_interact_uniref_a3m+sequence_based_template_pdb',
                                'species_interact_uniref_a3m+sequence_based_template_complex_pdb',
                                'species_interact_uniref_a3m+alphafold_model_templates',
                                'uniprot_distance_uniref_a3m',
                                'string_interact_uniref_a3m',
                                # 'geno_dist_uniref_a3m',
                                # 'pdb_interact_uniref_sto',
                                'species_interact_uniref_sto',
                                'uniprot_distance_uniref_sto',
                                'string_interact_uniref_sto',
                                'string_interact_uniref_sto+sequence_based_template_pdb70',
                                'string_interact_uniref_sto+structure_based_template',
                                'string_interact_uniref_sto+sequence_based_template_pdb',
                                'string_interact_uniref_sto+sequence_based_template_complex_pdb',
                                'string_interact_uniref_sto+alphafold_model_templates',
                                # 'geno_dist_uniref_sto',
                                # 'pdb_interact_uniprot_sto',
                                'species_interact_uniprot_sto',
                                'uniprot_distance_uniprot_sto',
                                'string_interact_uniprot_sto',
                                'species_colabfold_interact']
        else:
            self.run_methods = run_methods

        self.method2dir = {'default': 'default_multimer',
                           'default_mul_newest': 'default_mul_newest',
                           'default_uniclust30': 'default_uniclust30',
                           'default_uniref30_22': 'default_uniref30_22',
                           'default+structure_based_template': 'default_struct',
                           'default+sequence_based_template_pdb70': 'default_pdb70',
                           'default+sequence_based_template_pdb': 'default_pdb',
                           'default+sequence_based_template_complex_pdb': 'default_comp',
                           'default+alphafold_model_templates': 'default_af',
                           'uniclust_oxmatch_a3m': 'uniclust_oxmatch_a3m',
                           'pdb_interact_uniref_a3m': 'pdb_iter_uniref_a3m',

                           'species_interact_uniref_a3m': 'spec_iter_uniref_a3m',
                           'species_interact_uniref_a3m+structure_based_template': 'spec_struct',
                           'species_interact_uniref_a3m+sequence_based_template_pdb70': 'spec_pdb70',
                           'species_interact_uniref_a3m+sequence_based_template_pdb': 'spec_pdb',
                           'species_interact_uniref_a3m+sequence_based_template_complex_pdb': 'spec_comp',
                           'species_interact_uniref_a3m+alphafold_model_templates': 'spec_af',

                           'uniprot_distance_uniref_a3m': 'unidist_uniref_a3m',
                           'string_interact_uniref_a3m': 'str_iter_uniref_a3m',
                           'species_interact_uniref_sto': 'spec_iter_uniref_sto',
                           'uniprot_distance_uniref_sto': 'unidist_uniref_sto',

                           'string_interact_uniref_sto': 'str_iter_uniref_sto',
                           'string_interact_uniref_sto+structure_based_template': 'str_struct',
                           'string_interact_uniref_sto+sequence_based_template_pdb70': 'str_pdb70',
                           'string_interact_uniref_sto+sequence_based_template_pdb': 'str_pdb',
                           'string_interact_uniref_sto+sequence_based_template_complex_pdb': 'str_comp',
                           'string_interact_uniref_sto+alphafold_model_templates': 'str_af',

                           'species_interact_uniprot_sto': 'spec_iter_uniprot_sto',
                           'uniprot_distance_uniprot_sto': 'unidist_uniprot_sto',
                           'string_interact_uniprot_sto': 'str_iter_uniprot_sto',

                           'species_colabfold_interact': 'spec_colab_iter'}

    def process(self,
                fasta_path,
                chain_id_map,
                aln_dir,
                complex_aln_dir,
                template_dir,
                monomer_model_dir,
                output_dir):

        makedir_if_not_exists(output_dir)

        result_dirs = []

        while True:
            # run alphafold default pipeline:
            outdir = f"{output_dir}/default_multimer"
            monomers = [chain_id_map[chain_id].description for chain_id in chain_id_map]
            if not complete_result(outdir):
                os.chdir(self.params['alphafold_default_program_dir'])
                bfd_uniref_a3ms = []
                mgnify_stos = []
                uniref90_stos = []
                uniprot_stos = []
                for chain_id in chain_id_map:
                    monomer = chain_id_map[chain_id].description
                    monomer_bfd_uniref_a3m = f"{aln_dir}/{monomer}/{monomer}_uniref30_bfd.a3m"
                    if not os.path.exists(monomer_bfd_uniref_a3m):
                        raise Exception(f"Cannot find bfd and uniref30 a3m for {monomer}: {monomer_bfd_uniref_a3m}")
                    bfd_uniref_a3ms += [monomer_bfd_uniref_a3m]

                    monomer_mgnify_sto = f"{aln_dir}/{monomer}/{monomer}_mgnify.sto"
                    if not os.path.exists(monomer_mgnify_sto):
                        raise Exception(f"Cannot find mgnify sto for {monomer}: {monomer_mgnify_sto}")
                    mgnify_stos += [monomer_mgnify_sto]

                    monomer_uniref90_sto = f"{aln_dir}/{monomer}/{monomer}_uniref90.sto"
                    if not os.path.exists(monomer_uniref90_sto):
                        raise Exception(f"Cannot find uniref90 sto for {monomer}: {monomer_uniref90_sto}")
                    uniref90_stos += [monomer_uniref90_sto]

                    monomer_uniprot_sto = f"{aln_dir}/{monomer}/{monomer}_uniprot.sto"
                    if not os.path.exists(monomer_uniprot_sto):
                        raise Exception(f"Cannot find uniprot sto for {monomer}: {monomer_uniprot_sto}")
                    uniprot_stos += [monomer_uniprot_sto]

                cmd = f"python {self.params['alphafold_default_program']} " \
                      f"--fasta_path {fasta_path} " \
                      f"--bfd_uniclust_a3ms {','.join(bfd_uniref_a3ms)} " \
                      f"--mgnify_stos {','.join(mgnify_stos)} " \
                      f"--uniref90_stos {','.join(uniref90_stos)} " \
                      f"--uniprot_stos {','.join(uniprot_stos)} " \
                      f"--env_dir {self.params['alphafold_env_dir']} " \
                      f"--database_dir {self.params['alphafold_database_dir']} " \
                      f"--output_dir {outdir}"

                print(cmd)
                os.system(cmd)
            result_dirs += [outdir]

            if "default_uniref30_22" in self.run_methods:
                # run alphafold default pipeline:
                outdir = f"{output_dir}/default_uniref30_22"
                monomers = [chain_id_map[chain_id].description for chain_id in chain_id_map]
                if not complete_result(outdir):
                    os.chdir(self.params['alphafold_default_program_dir'])
                    bfd_uniclust_a3ms = []
                    mgnify_stos = []
                    uniref90_stos = []
                    uniprot_stos = []
                    for chain_id in chain_id_map:
                        monomer = chain_id_map[chain_id].description
                        monomer_bfd_uniclust_a3m = f"{aln_dir}/{monomer}/{monomer}_uniref30_22_bfd.a3m"
                        if not os.path.exists(monomer_bfd_uniclust_a3m):
                            raise Exception(f"Cannot find bfd and uniclust a3m for {monomer}: {monomer_bfd_uniclust_a3m}")
                        bfd_uniclust_a3ms += [monomer_bfd_uniclust_a3m]

                        monomer_mgnify_sto = f"{aln_dir}/{monomer}/{monomer}_mgnify.sto"
                        if not os.path.exists(monomer_mgnify_sto):
                            raise Exception(f"Cannot find mgnify sto for {monomer}: {monomer_mgnify_sto}")
                        mgnify_stos += [monomer_mgnify_sto]

                        monomer_uniref90_sto = f"{aln_dir}/{monomer}/{monomer}_uniref90.sto"
                        if not os.path.exists(monomer_uniref90_sto):
                            raise Exception(f"Cannot find uniref90 sto for {monomer}: {monomer_uniref90_sto}")
                        uniref90_stos += [monomer_uniref90_sto]

                        monomer_uniprot_sto = f"{aln_dir}/{monomer}/{monomer}_uniprot.sto"
                        if not os.path.exists(monomer_uniprot_sto):
                            raise Exception(f"Cannot find uniprot sto for {monomer}: {monomer_uniprot_sto}")
                        uniprot_stos += [monomer_uniprot_sto]

                    cmd = f"python {self.params['alphafold_default_program']} " \
                          f"--fasta_path {fasta_path} " \
                          f"--bfd_uniclust_a3ms {','.join(bfd_uniclust_a3ms)} " \
                          f"--mgnify_stos {','.join(mgnify_stos)} " \
                          f"--uniref90_stos {','.join(uniref90_stos)} " \
                          f"--uniprot_stos {','.join(uniprot_stos)} " \
                          f"--num_multimer_predictions_per_model 5 " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir']} " \
                          f"--output_dir {outdir}"

                    print(cmd)
                    os.system(cmd)

                result_dirs += [outdir]

            if "default_mul_newest" in self.run_methods:
                # run alphafold default pipeline:
                outdir = f"{output_dir}/default_mul_newest"
                monomers = [chain_id_map[chain_id].description for chain_id in chain_id_map]
                if not complete_result(outdir):
                    os.chdir(self.params['alphafold_default_program_dir'])
                    bfd_uniclust_a3ms = []
                    mgnify_stos = []
                    uniref90_stos = []
                    uniprot_stos = []
                    for chain_id in chain_id_map:
                        monomer = chain_id_map[chain_id].description
                        monomer_bfd_uniclust_a3m = f"{aln_dir}/{monomer}/{monomer}_uniref30_22_bfd.a3m"
                        if not os.path.exists(monomer_bfd_uniclust_a3m):
                            raise Exception(f"Cannot find bfd and uniclust a3m for {monomer}: {monomer_bfd_uniclust_a3m}")
                        bfd_uniclust_a3ms += [monomer_bfd_uniclust_a3m]

                        monomer_mgnify_sto = f"{aln_dir}/{monomer}/{monomer}_mgnify.sto"
                        if not os.path.exists(monomer_mgnify_sto):
                            raise Exception(f"Cannot find mgnify sto for {monomer}: {monomer_mgnify_sto}")
                        mgnify_stos += [monomer_mgnify_sto]

                        monomer_uniref90_sto = f"{aln_dir}/{monomer}/{monomer}_uniref90_new.sto"
                        if not os.path.exists(monomer_uniref90_sto):
                            raise Exception(f"Cannot find uniref90 sto for {monomer}: {monomer_uniref90_sto}")
                        uniref90_stos += [monomer_uniref90_sto]

                        monomer_uniprot_sto = f"{aln_dir}/{monomer}/{monomer}_uniprot_new.sto"
                        if not os.path.exists(monomer_uniprot_sto):
                            raise Exception(f"Cannot find uniprot sto for {monomer}: {monomer_uniprot_sto}")
                        uniprot_stos += [monomer_uniprot_sto]

                    cmd = f"python {self.params['alphafold_default_program']} " \
                          f"--fasta_path {fasta_path} " \
                          f"--bfd_uniclust_a3ms {','.join(bfd_uniclust_a3ms)} " \
                          f"--mgnify_stos {','.join(mgnify_stos)} " \
                          f"--uniref90_stos {','.join(uniref90_stos)} " \
                          f"--uniprot_stos {','.join(uniprot_stos)} " \
                          f"--num_multimer_predictions_per_model 5 " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir_newest']} " \
                          f"--output_dir {outdir}"

                    print(cmd)
                    os.system(cmd)

                result_dirs += [outdir]

            # run alphafold default pipeline:
            if "default_uniclust30" in self.run_methods:
                outdir = f"{output_dir}/default_uniclust30"
                monomers = [chain_id_map[chain_id].description for chain_id in chain_id_map]
                if not complete_result(outdir):
                    os.chdir(self.params['alphafold_default_program_dir'])
                    bfd_uniclust_a3ms = []
                    mgnify_stos = []
                    uniref90_stos = []
                    uniprot_stos = []
                    for chain_id in chain_id_map:
                        monomer = chain_id_map[chain_id].description
                        monomer_bfd_uniclust_a3m = f"{aln_dir}/{monomer}/{monomer}_uniclust30_bfd.a3m"
                        if not os.path.exists(monomer_bfd_uniclust_a3m):
                            raise Exception(f"Cannot find bfd and uniclust a3m for {monomer}: {monomer_bfd_uniclust_a3m}")
                        bfd_uniclust_a3ms += [monomer_bfd_uniclust_a3m]

                        monomer_mgnify_sto = f"{aln_dir}/{monomer}/{monomer}_mgnify.sto"
                        if not os.path.exists(monomer_mgnify_sto):
                            raise Exception(f"Cannot find mgnify sto for {monomer}: {monomer_mgnify_sto}")
                        mgnify_stos += [monomer_mgnify_sto]

                        monomer_uniref90_sto = f"{aln_dir}/{monomer}/{monomer}_uniref90.sto"
                        if not os.path.exists(monomer_uniref90_sto):
                            raise Exception(f"Cannot find uniref90 sto for {monomer}: {monomer_uniref90_sto}")
                        uniref90_stos += [monomer_uniref90_sto]

                        monomer_uniprot_sto = f"{aln_dir}/{monomer}/{monomer}_uniprot.sto"
                        if not os.path.exists(monomer_uniprot_sto):
                            raise Exception(f"Cannot find uniprot sto for {monomer}: {monomer_uniprot_sto}")
                        uniprot_stos += [monomer_uniprot_sto]

                    cmd = f"python {self.params['alphafold_default_program']} " \
                          f"--fasta_path {fasta_path} " \
                          f"--bfd_uniclust_a3ms {','.join(bfd_uniclust_a3ms)} " \
                          f"--mgnify_stos {','.join(mgnify_stos)} " \
                          f"--uniref90_stos {','.join(uniref90_stos)} " \
                          f"--uniprot_stos {','.join(uniprot_stos)} " \
                          f"--env_dir {self.params['alphafold_env_dir']} " \
                          f"--database_dir {self.params['alphafold_database_dir']} " \
                          f"--output_dir {outdir}"

                    print(cmd)
                    os.system(cmd)
                result_dirs += [outdir]
            
            # Customized complex alignment pipelines using original template search pipeline in alphafold
            default_alphafold_monomer_a3ms = []
            default_alphafold_multimer_a3ms = []
            template_stos = []
            for chain_id in chain_id_map:
                monomer = chain_id_map[chain_id].description
                monomer_template_sto = f"{aln_dir}/{monomer}/{monomer}_uniref90.sto"
                if not os.path.exists(monomer_template_sto):
                    raise Exception(f"Cannot find template stos for {monomer}: {monomer_template_sto}")
                template_stos += [monomer_template_sto]

                default_alphafold_monomer_a3m = f"{output_dir}/default_multimer/msas/{chain_id}/monomer_final.a3m"
                if not os.path.exists(default_alphafold_monomer_a3m):
                    raise Exception(
                        f"Cannot find default alphafold alignments for {monomer}: {default_alphafold_monomer_a3m}")
                default_alphafold_monomer_a3ms += [default_alphafold_monomer_a3m]

                default_alphafold_multimer_a3m = f"{output_dir}/default_multimer/msas/{monomer}.paired.a3m"
                if not os.path.exists(default_alphafold_monomer_a3m):
                    raise Exception(
                        f"Cannot find default alphafold alignments for {monomer}: {default_alphafold_multimer_a3m}")
                default_alphafold_multimer_a3ms += [default_alphafold_multimer_a3m]

            for method in self.run_methods:

                os.chdir(self.params['alphafold_program_dir'])

                if method == "default" or method == "default_uniclust30" or method == "default_uniref30_22" or method == "default_mul_newest":
                    continue
                concatenate_method = ""
                template_method = ""
                if method.find('+') > 0:
                    concatenate_method, template_method = method.split('+')
                else:
                    concatenate_method = method

                if concatenate_method == "default":
                    a3m_paths = default_alphafold_multimer_a3ms
                    msa_pair_file = f"{output_dir}/default_multimer/msas/interact.csv"
                    interact_dict = {}
                    msa_len = -1
                    for i in range(len(default_alphafold_multimer_a3ms)):
                        with open(default_alphafold_multimer_a3ms[i]) as f:
                            input_fasta_str = f.read()
                        msa_sequences, msa_descriptions = parse_fasta(input_fasta_str)
                        current_len = len(msa_descriptions)
                        if msa_len == -1:
                            msa_len = current_len
                        elif current_len != msa_len:
                            raise Exception(f"The length of each msas are not equal! {default_alphafold_multimer_a3ms}")
                        interact_dict[f'index_{i + 1}'] = [j for j in range(msa_len)]
                    interact_df = pd.DataFrame(interact_dict)
                    interact_df.to_csv(msa_pair_file)
                elif concatenate_method == "species_colabfold_interact":
                    species_msa_pair_file = f"{complex_aln_dir}/species_interact_uniref_a3m/species_interact_uniref_a3m_interact.csv"
                    msa_pair_file = f"{complex_aln_dir}/{concatenate_method}/{concatenate_method}_interact.csv"
                    if len(pd.read_csv(msa_pair_file)) <= 1 or \
                            len(pd.read_csv(msa_pair_file)) == len(pd.read_csv(species_msa_pair_file)) + 1:
                        continue
                    a3m_paths = [f"{complex_aln_dir}/{concatenate_method}/{monomer}_con.a3m" for monomer in monomers]
                else:
                    msa_pair_file = f"{complex_aln_dir}/{concatenate_method}/{concatenate_method}_interact.csv"
                    if len(pd.read_csv(msa_pair_file)) <= 1:
                        continue
                    a3m_paths = [f"{complex_aln_dir}/{concatenate_method}/{monomer}_con.a3m" for monomer in monomers]
                        # get_complex_alignments_by_method(
                        # monomers=monomers,
                        # concatenate_method=concatenate_method,
                        # aln_dir=aln_dir)

                outdir = f"{output_dir}/{self.method2dir[method]}"

                base_cmd = f"python {self.params['alphafold_multimer_program']} " \
                           f"--fasta_path {fasta_path} " \
                           f"--monomer_a3ms {','.join(default_alphafold_monomer_a3ms)} " \
                           f"--multimer_a3ms {','.join(a3m_paths)} " \
                           f"--msa_pair_file {msa_pair_file} " \
                           f"--env_dir {self.params['alphafold_env_dir']} " \
                           f"--database_dir {self.params['alphafold_database_dir']} "

                if template_method == "":
                    base_cmd += f"--template_stos {','.join(template_stos)} "

                elif template_method == "structure_based_template":
                    template_file = f"{template_dir}/struct_temp/structure_templates.csv"
                    if len(pd.read_csv(template_file)) < 4:
                        template_file = f"{template_dir}/struct_temp/structure_templates_v2.csv"
                        if len(pd.read_csv(template_file)) == 0:
                            continue
                        os.chdir(self.params['alphafold_temp_program_dir'])
                        outdir += '_v2'

                    base_cmd += f"--temp_struct_csv {template_file} "
                    base_cmd += f"--struct_atom_dir {template_dir}/struct_temp/templates "

                elif template_method == "sequence_based_template_pdb":
                    template_file = f"{template_dir}/pdb_seq/sequence_templates.csv"
                    if len(pd.read_csv(template_file)) < 4:
                        template_file = f"{template_dir}/pdb_seq/sequence_templates_v2.csv"
                        if len(pd.read_csv(template_file)) == 0:
                            continue
                        os.chdir(self.params['alphafold_temp_program_dir'])
                        outdir += '_v2'

                    base_cmd += f"--temp_struct_csv {template_file} "
                    base_cmd += f"--struct_atom_dir {template_dir}/pdb_seq/templates "

                elif template_method == "sequence_based_template_complex_pdb":
                    template_file = f"{template_dir}/complex_pdb_seq/sequence_templates.csv"
                    if len(pd.read_csv(template_file)) < 4:
                        template_file = f"{template_dir}/complex_pdb_seq/sequence_templates_v2.csv"
                        if len(pd.read_csv(template_file)) == 0:
                            continue
                        os.chdir(self.params['alphafold_temp_program_dir'])
                        outdir += '_v2'

                    base_cmd += f"--temp_struct_csv {template_file} "
                    base_cmd += f"--struct_atom_dir {template_dir}/complex_pdb_seq/templates "

                elif template_method == "sequence_based_template_pdb70":
                    template_file = f"{template_dir}/pdb70_seq/sequence_templates.csv"
                    if len(pd.read_csv(template_file)) < 4:
                        template_file = f"{template_dir}/pdb70_seq/sequence_templates_v2.csv"
                        if len(pd.read_csv(template_file)) == 0:
                            continue
                        os.chdir(self.params['alphafold_temp_program_dir'])
                        outdir += '_v2'

                    template_hits_files = []
                    for monomer in monomers:
                        template_hits_file = f"{template_dir}/pdb70_seq/{monomer}/output.hhr"
                        if not os.path.exists(template_hits_file):
                            raise Exception(f"Cannot find template hit file for {monomer}: {template_hits_file}")
                        template_hits_files += [template_hits_file]
                    base_cmd += f"--temp_seq_pair_file {template_file} "
                    base_cmd += f"--template_hits_files {','.join(template_hits_files)} "
                    
                elif template_method == "alphafold_model_templates":
                    monomer_paths = []
                    for monomer in monomers:
                        monomer_path = f"{monomer_model_dir}/{monomer}/default"
                        if not os.path.exists(monomer_path):
                            raise Exception(f"Cannot find monomer directory for {monomer}: {monomer_path}")
                        monomer_paths += [monomer_path]
                    base_cmd += f"--monomer_model_paths {','.join(monomer_paths)} "

                base_cmd += f"--output_dir {outdir} "
                if complete_result(outdir):
                    continue

                makedir_if_not_exists(outdir)

                print(base_cmd)
                os.system(base_cmd)

                result_dirs += [outdir]

            rerun = False
            for result_dir in result_dirs:
                if not complete_result(result_dir):
                    rerun = True

            if not rerun:
                break

            print("end")

        print("The quaternary structure generation for multimers has finished!")
