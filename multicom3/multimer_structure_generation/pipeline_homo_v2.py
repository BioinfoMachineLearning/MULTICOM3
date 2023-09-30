import copy
import os
import sys
import time
from multicom3.common.util import makedir_if_not_exists, check_dirs
from multicom3.common.protein import complete_result, parse_fasta
import pandas as pd
from multiprocessing import Pool
import pathlib


def get_complex_alignments_by_method(monomers, concatenate_method, aln_dir):
    a3ms_path = []
    for monomer in monomers:
        if concatenate_method == 'uniclust_oxmatch_a3m':
            monomer_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniclust30.a3m"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniref_a3m') > 0:
            monomer_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniref30.a3m"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniref_sto') > 0:
            monomer_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniref90.sto"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
        elif concatenate_method.find('_uniprot_sto') > 0:
            monomer_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniprot.sto"
            if not os.path.exists(monomer_a3m):
                raise Exception(f"Cannot find alignment for {monomer}: {monomer_a3m}")
            a3ms_path += [monomer_a3m]
    return a3ms_path


def combine_a3ms(infiles, outfile):
    targetname = None
    targetseq = None
    descriptions = []
    seqs = []
    for infile in infiles:
        for line in open(infile):
            line = line.rstrip('\n')
            if line.startswith('>'):
                descriptions += [line]
                if targetname is None:
                    targetname = line
            else:
                seqs += [line]
                if targetseq is None:
                    targetseq = line

    seen_seqs = []
    with open(outfile, 'w') as fw:
        fw.write(f"{targetname}\n{targetseq}\n")
        for (desc, seq) in zip(descriptions, seqs):
            if desc == targetname and seq == targetseq:
                continue
            if seq not in seen_seqs:
                fw.write(f"{desc}\n{seq}\n")
                seen_seqs.append(seq)


class Multimer_structure_prediction_homo_pipeline_v2:

    def __init__(self, params, run_methods=None):

        self.params = params

        if run_methods is None:
            self.run_methods = ['default',
                                'default+sequence_based_template_pdb70',
                                'default+structure_based_template',
                                'default+sequence_based_template_pdb',
                                'default+sequence_based_template_complex_pdb',
                                # 'default+alphafold_model_templates',
                                'uniclust_oxmatch_a3m',
                                'pdb_interact_uniref_a3m',
                                'species_interact_uniref_a3m',
                                'species_interact_uniref_a3m+sequence_based_template_pdb70',
                                'species_interact_uniref_a3m+structure_based_template',
                                'species_interact_uniref_a3m+sequence_based_template_pdb',
                                'species_interact_uniref_a3m+sequence_based_template_complex_pdb',
                                # 'species_interact_uniref_a3m+alphafold_model_templates',
                                'species_interact_uniref_sto',
                                'species_interact_uniprot_sto']
        else:
            self.run_methods = run_methods


        self.method2dir = {'default': 'default_multimer',
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
                           'string_interact_uniprot_sto': 'str_iter_uniprot_sto'}

    def process_img(self,
                    fasta_path,
                    chain_id_map,
                    aln_dir,
                    output_dir):

        outdir = os.path.join(output_dir, "default_img")
        if not complete_result(outdir, 5 * int(self.params['num_multimer_predictions_per_model'])):
            monomers = []
            img_default_a3ms = []
            template_stos = []
            for chain_id in chain_id_map:
                monomer = chain_id
                monomers += [monomer]
                img_a3m = os.path.join(aln_dir, monomer, f"{monomer}.a3m")
                default_a3m = os.path.join(output_dir, 'default_multimer', 'msas', chain_id, "monomer_final.a3m")
                combine_a3ms([default_a3m, img_a3m], os.path.join(aln_dir, monomer, f"{monomer}_default_img.a3m"))
                img_default_a3ms += [os.path.join(aln_dir, monomer, f"{monomer}_default_img.a3m")]
                template_stos += [os.path.join(aln_dir, monomer, f"{monomer}_uniref90.sto")]

            makedir_if_not_exists(outdir)

            base_cmd = f"python {self.params['alphafold_multimer_program']} " \
                       f"--fasta_path {fasta_path} " \
                       f"--monomer_a3ms {','.join(img_default_a3ms)} " \
                       f"--multimer_a3ms {','.join(img_default_a3ms)} " \
                       f"--template_stos {','.join(template_stos)} " \
                       f"--env_dir {self.params['alphafold_env_dir']} " \
                       f"--database_dir {self.params['alphafold_database_dir']} " \
                       f"--num_multimer_predictions_per_model {self.params['num_multimer_predictions_per_model']} " \
                       f"--multimer_num_ensemble {self.params['multimer_num_ensemble']} " \
                       f"--multimer_num_recycle {self.params['multimer_num_recycle']} " \
                       f"--output_dir {outdir} "

            print(base_cmd)
            os.system(base_cmd)

    def process(self,
                fasta_path,
                chain_id_map,
                aln_dir,
                complex_aln_dir,
                template_dir,
                monomer_model_dir,
                output_dir):

        makedir_if_not_exists(output_dir)

        # run alphafold default pipeline:
        outdir = os.path.join(output_dir, "default_multimer")
        monomers = [chain_id for chain_id in chain_id_map]
        if not complete_result(outdir, 5 * int(self.params['num_multimer_predictions_per_model'])):
            os.chdir(self.params['alphafold_program_dir'])
            bfd_uniref_a3ms = []
            mgnify_stos = []
            uniref90_stos = []
            uniprot_stos = []
            for chain_id in chain_id_map:
                monomer = chain_id
                monomer_bfd_uniref_a3m = os.path.join(aln_dir, monomer, f"{monomer}_uniref30_bfd.a3m")
                if not os.path.exists(monomer_bfd_uniref_a3m):
                    raise Exception(f"Cannot find bfd and uniclust a3m for {monomer}: {monomer_bfd_uniref_a3m}")
                bfd_uniref_a3ms += [monomer_bfd_uniref_a3m]

                monomer_mgnify_sto = os.path.join(aln_dir, monomer, f"{monomer}_mgnify.sto")
                if not os.path.exists(monomer_mgnify_sto):
                    raise Exception(f"Cannot find mgnify sto for {monomer}: {monomer_mgnify_sto}")
                mgnify_stos += [monomer_mgnify_sto]

                monomer_uniref90_sto = os.path.join(aln_dir, monomer, f"{monomer}_uniref90.sto")
                if not os.path.exists(monomer_uniref90_sto):
                    raise Exception(f"Cannot find uniref90 sto for {monomer}: {monomer_uniref90_sto}")
                uniref90_stos += [monomer_uniref90_sto]

                monomer_uniprot_sto = os.path.join(aln_dir, monomer, f"{monomer}_uniprot.sto")
                if not os.path.exists(monomer_uniprot_sto):
                    raise Exception(f"Cannot find uniprot sto for {monomer}: {monomer_uniprot_sto}")
                uniprot_stos += [monomer_uniprot_sto]

            cmd = f"python {self.params['alphafold_default_program']} " \
                  f"--fasta_path {fasta_path} " \
                  f"--bfd_uniref_a3ms {','.join(bfd_uniref_a3ms)} " \
                  f"--mgnify_stos {','.join(mgnify_stos)} " \
                  f"--uniref90_stos {','.join(uniref90_stos)} " \
                  f"--uniprot_stos {','.join(uniprot_stos)} " \
                  f"--env_dir {self.params['alphafold_env_dir']} " \
                  f"--database_dir {self.params['alphafold_database_dir']} " \
                  f"--num_multimer_predictions_per_model {self.params['num_multimer_predictions_per_model']} " \
                  f"--multimer_num_ensemble {self.params['multimer_num_ensemble']} " \
                  f"--multimer_num_recycle {self.params['multimer_num_recycle']} " \
                  f"--output_dir {outdir}"

            print(cmd)
            os.system(cmd)

        os.chdir(self.params['alphafold_program_dir'])

        # Customized complex alignment pipelines using original template search pipeline in alphafold
        default_alphafold_monomer_a3ms = []
        template_stos = []
        for chain_id in chain_id_map:
            monomer = chain_id
            monomer_template_sto = os.path.join(aln_dir, monomer, f"{monomer}_uniref90.sto")
            if not os.path.exists(monomer_template_sto):
                raise Exception(f"Cannot find template stos for {monomer}: {monomer_template_sto}")
            template_stos += [monomer_template_sto]

            default_alphafold_monomer_a3m = os.path.join(output_dir, 'default_multimer', 'msas', chain_id, "monomer_final.a3m")
            if not os.path.exists(default_alphafold_monomer_a3m):
                raise Exception(
                    f"Cannot find default alphafold alignments for {monomer}: {default_alphafold_monomer_a3m}")
            default_alphafold_monomer_a3ms += [default_alphafold_monomer_a3m]

        for method in self.run_methods:
            if method == "default" or method == "default_uniclust30" or method == "default_uniref30_22" or method == "default_mul_newest":
                continue
            concatenate_method = ""
            template_method = ""
            if method.find('+') > 0:
                concatenate_method, template_method = method.split('+')
            else:
                concatenate_method = method

            if concatenate_method == "default":
                a3m_paths = default_alphafold_monomer_a3ms
            else:
                msa_pair_file = os.path.join(complex_aln_dir, concatenate_method, concatenate_method + '_interact.csv')
                if len(pd.read_csv(msa_pair_file)) <= 1:
                    continue
                a3m_paths = [os.path.join(complex_aln_dir, concatenate_method, f"{monomer}_con.a3m") for monomer in monomers]
                print(a3m_paths)

            outdir = os.path.join(output_dir, self.method2dir[method])

            if complete_result(outdir, 5 * int(self.params['num_multimer_predictions_per_model'])):
                continue

            makedir_if_not_exists(outdir)

            base_cmd = f"python {self.params['alphafold_multimer_program']} " \
                       f"--fasta_path {fasta_path} " \
                       f"--monomer_a3ms {','.join(a3m_paths)} " \
                       f"--multimer_a3ms {','.join(a3m_paths)} " \
                       f"--env_dir {self.params['alphafold_env_dir']} " \
                       f"--database_dir {self.params['alphafold_database_dir']} " \
                       f"--num_multimer_predictions_per_model {self.params['num_multimer_predictions_per_model']} " \
                       f"--multimer_num_ensemble {self.params['multimer_num_ensemble']} " \
                       f"--multimer_num_recycle {self.params['multimer_num_recycle']} " \
                       f"--output_dir {outdir} "

            if template_method == "":
                base_cmd += f"--template_stos {','.join(template_stos)} "

            elif template_method == "structure_based_template":
                template_file = os.path.join(template_dir, "struct_temp", "structure_templates.csv")
                struct_atom_dir = os.path.join(template_dir, "struct_temp", "templates")
                if len(pd.read_csv(template_file)) == 0:
                    continue
                base_cmd += f"--temp_struct_csv {template_file} "
                base_cmd += f"--struct_atom_dir {struct_atom_dir} "

            elif template_method == "sequence_based_template_pdb":
                template_file = os.path.join(template_dir, "pdb_seq", "sequence_templates.csv")
                struct_atom_dir = os.path.join(template_dir, "pdb_seq", "templates")
                if len(pd.read_csv(template_file)) == 0:
                    continue
                base_cmd += f"--temp_struct_csv {template_file} "
                base_cmd += f"--struct_atom_dir {struct_atom_dir} "


            elif template_method == "sequence_based_template_complex_pdb":
                template_file = os.path.join(template_dir, "complex_pdb_seq", "sequence_templates.csv")
                struct_atom_dir = os.path.join(template_dir, "complex_pdb_seq", "templates")
                if len(pd.read_csv(template_file)) == 0:
                    continue
                base_cmd += f"--temp_struct_csv {template_file} "
                base_cmd += f"--struct_atom_dir {struct_atom_dir} "

            elif template_method == "sequence_based_template_pdb70":
                template_file = os.path.join(template_dir, "pdb70_seq", "sequence_templates.csv")
                if len(pd.read_csv(template_file)) == 0:
                    continue
                template_hits_files = []
                for monomer in monomers:
                    template_hits_file = os.path.join(template_dir, "pdb70_seq", monomer, "output.hhr")
                    if not os.path.exists(template_hits_file):
                        raise Exception(f"Cannot find template hit file for {monomer}: {template_hits_file}")
                    template_hits_files += [template_hits_file]
                base_cmd += f"--temp_seq_pair_file {template_file} "
                base_cmd += f"--template_hits_files {','.join(template_hits_files)} "

            elif template_method == "alphafold_model_templates":
                monomer_paths = []
                for monomer in monomers:
                    monomer_path = os.path.join(monomer_model_dir, monomer, "default")
                    if not os.path.exists(monomer_path):
                        raise Exception(f"Cannot find monomer directory for {monomer}: {monomer_path}")
                    monomer_paths += [monomer_path]
                base_cmd += f"--monomer_model_paths {','.join(monomer_paths)} "

            if complete_result(outdir, 5 * int(self.params['num_multimer_predictions_per_model'])):
                continue

            print(base_cmd)
            os.system(base_cmd)

        print("The multimer structure generation for multimers has finished!")
