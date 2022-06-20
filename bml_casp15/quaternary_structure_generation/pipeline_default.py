import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
from bml_casp15.common.protein import complete_result, parse_fasta
import pandas as pd
from multiprocessing import Pool
import pathlib


class Quaternary_structure_prediction_pipeline_default:

    def __init__(self, params):  # , run_methods):

        self.params = params

    def process(self,
                fasta_path,
                chain_id_map,
                aln_dir,
                output_dir):

        makedir_if_not_exists(output_dir)

        # run alphafold default pipeline:
        # outdir = f"{output_dir}/default_multimer"
        # monomers = [chain_id_map[chain_id].description for chain_id in chain_id_map]
        # if not complete_result(outdir):
        #     os.chdir(self.params['alphafold_default_program_dir'])
        #     bfd_uniclust_a3ms = []
        #     mgnify_stos = []
        #     uniref90_stos = []
        #     uniprot_stos = []
        #     for chain_id in chain_id_map:
        #         monomer = chain_id_map[chain_id].description
        #         monomer_bfd_uniclust_a3m = f"{aln_dir}/{monomer}/{monomer}_uniref30_bfd.a3m"
        #         if not os.path.exists(monomer_bfd_uniclust_a3m):
        #             raise Exception(f"Cannot find bfd and uniclust a3m for {monomer}: {monomer_bfd_uniclust_a3m}")
        #         bfd_uniclust_a3ms += [monomer_bfd_uniclust_a3m]
        #
        #         monomer_mgnify_sto = f"{aln_dir}/{monomer}/{monomer}_mgnify.sto"
        #         if not os.path.exists(monomer_mgnify_sto):
        #             raise Exception(f"Cannot find mgnify sto for {monomer}: {monomer_mgnify_sto}")
        #         mgnify_stos += [monomer_mgnify_sto]
        #
        #         monomer_uniref90_sto = f"{aln_dir}/{monomer}/{monomer}_uniref90.sto"
        #         if not os.path.exists(monomer_uniref90_sto):
        #             raise Exception(f"Cannot find uniref90 sto for {monomer}: {monomer_uniref90_sto}")
        #         uniref90_stos += [monomer_uniref90_sto]
        #
        #         monomer_uniprot_sto = f"{aln_dir}/{monomer}/{monomer}_uniprot.sto"
        #         if not os.path.exists(monomer_uniprot_sto):
        #             raise Exception(f"Cannot find uniprot sto for {monomer}: {monomer_uniprot_sto}")
        #         uniprot_stos += [monomer_uniprot_sto]
        #
        #     cmd = f"python {self.params['alphafold_default_program']} " \
        #           f"--fasta_path {fasta_path} " \
        #           f"--bfd_uniclust_a3ms {','.join(bfd_uniclust_a3ms)} " \
        #           f"--mgnify_stos {','.join(mgnify_stos)} " \
        #           f"--uniref90_stos {','.join(uniref90_stos)} " \
        #           f"--uniprot_stos {','.join(uniprot_stos)} " \
        #           f"--env_dir {self.params['alphafold_env_dir']} " \
        #           f"--database_dir {self.params['alphafold_database_dir']} " \
        #           f"--output_dir {outdir}"
        #
        #     print(cmd)
        #     os.system(cmd)

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
                  f"--database_dir {self.params['alphafold_database_dir']} " \
                  f"--output_dir {outdir}"

            print(cmd)
            os.system(cmd)

        print("The quaternary structure generation for multimers has finished!")
