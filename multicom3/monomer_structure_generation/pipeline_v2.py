import copy
import os
import sys
import time
from multicom3.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import pathlib
from multicom3.common.protein import complete_result


class Monomer_structure_prediction_pipeline_v2:

    def __init__(self, params, run_methods=None):

        self.params = params

        if run_methods is None:
            self.run_methods = ['default', 'default+seq_template',
                                'original', 'original+seq_template',
                                'colabfold', 'colabfold+seq_template']
        else:
            self.run_methods = run_methods

    def process_single(self, fasta_path, alndir, outdir, template_dir=None):
        
        targetname = pathlib.Path(fasta_path).stem

        cmd = ""

        makedir_if_not_exists(outdir)

        if "default" in self.run_methods:

            os.chdir(self.params['alphafold_program_dir'])

            errormsg = ""

            if not os.path.exists(alndir):
                errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

            bfd_uniref30_a3m = alndir + '/' + targetname + '_uniref30_bfd.a3m'
            if not os.path.exists(bfd_uniref30_a3m):
                errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"

            mgnify_sto = alndir + '/' + targetname + '_mgnify.sto'
            if not os.path.exists(mgnify_sto):
                errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

            uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
            if not os.path.exists(uniref90_sto):
                errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            if len(errormsg) == 0:
                if not complete_result(f"{outdir}/default", 5 * int(self.params['num_monomer_predictions_per_model'])):
                    try:
                        cmd = f"python {self.params['alphafold_default_program']} " \
                              f"--fasta_path {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--bfd_uniref_a3ms {bfd_uniref30_a3m} " \
                              f"--mgnify_stos {mgnify_sto} " \
                              f"--uniref90_stos {uniref90_sto} " \
                              f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--output_dir {outdir}/default"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)
            else:
                print(errormsg)

        if "default+seq_template" in self.run_methods:

            os.chdir(self.params['alphafold_program_dir'])

            errormsg = ""

            if not os.path.exists(alndir):
                errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

            bfd_uniref30_a3m = alndir + '/' + targetname + '_uniref30_bfd.a3m'
            if not os.path.exists(bfd_uniref30_a3m):
                errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {bfd_uniref30_a3m}\n"

            mgnify_sto = alndir + '/' + targetname + '_mgnify.sto'
            if not os.path.exists(mgnify_sto):
                errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

            uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
            if not os.path.exists(uniref90_sto):
                errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            if len(errormsg) == 0:
                if not complete_result(f"{outdir}/default_seq_temp", 5 * int(self.params['num_monomer_predictions_per_model'])):
                    try:
                        cmd = f"python {self.params['alphafold_program']} " \
                              f"--fasta_path {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--bfd_uniref_a3m {bfd_uniref30_a3m} " \
                              f"--mgnify_sto {mgnify_sto} " \
                              f"--uniref90_sto {uniref90_sto} " \
                              f"--temp_struct_csv {template_dir}/sequence_templates.csv " \
                              f"--struct_atom_dir {template_dir}/templates " \
                              f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--output_dir {outdir}/default_seq_temp"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)
            else:
                print(errormsg)

        if "original" in self.run_methods:

            os.chdir(self.params['alphafold_program_dir'])

            errormsg = ""

            if not os.path.exists(alndir):
                errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

            uniref30_a3m = alndir + '/' + targetname + '_uniref30.a3m'
            if not os.path.exists(uniref30_a3m):
                errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {uniref30_a3m}\n"

            bfd_a3m = alndir + '/' + targetname + '_bfd.a3m'
            if not os.path.exists(bfd_a3m):
                errormsg = errormsg + f"Cannot find bfd alignment for {targetname}: {bfd_a3m}\n"

            mgnify_sto = alndir + '/' + targetname + '_mgnify.sto'
            if not os.path.exists(mgnify_sto):
                errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

            uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
            if not os.path.exists(uniref90_sto):
                errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            if len(errormsg) == 0:
                if not complete_result(f"{outdir}/original", 5 * int(self.params['num_monomer_predictions_per_model'])):
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_path {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--bfd_uniref_a3m {uniref30_a3m} " \
                              f"--bfd_a3m {bfd_a3m} " \
                              f"--mgnify_sto {mgnify_sto} " \
                              f"--uniref90_sto {uniref90_sto} " \
                              f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--output_dir {outdir}/original"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)
            else:
                print(errormsg)

        if "original+seq_template" in self.run_methods:

            os.chdir(self.params['alphafold_program_dir'])

            errormsg = ""

            if not os.path.exists(alndir):
                errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {alndir}\n"

            uniref30_a3m = alndir + '/' + targetname + '_uniref30.a3m'
            if not os.path.exists(uniref30_a3m):
                errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {uniref30_a3m}\n"

            bfd_a3m = alndir + '/' + targetname + '_bfd.a3m'
            if not os.path.exists(bfd_a3m):
                errormsg = errormsg + f"Cannot find bfd alignment for {targetname}: {bfd_a3m}\n"

            mgnify_sto = alndir + '/' + targetname + '_mgnify.sto'
            if not os.path.exists(mgnify_sto):
                errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

            uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
            if not os.path.exists(uniref90_sto):
                errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            if len(errormsg) == 0:
                if not complete_result(f"{outdir}/ori_seq_temp", 5 * int(self.params['num_monomer_predictions_per_model'])):
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_path {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--bfd_uniref_a3m {uniref30_a3m} " \
                              f"--bfd_a3m {bfd_a3m} " \
                              f"--mgnify_sto {mgnify_sto} " \
                              f"--uniref90_sto {uniref90_sto} " \
                              f"--temp_struct_csv {template_dir}/sequence_templates.csv " \
                              f"--struct_atom_dir {template_dir}/templates " \
                              f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--output_dir {outdir}/ori_seq_temp"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)
            else:
                print(errormsg)

        if "colabfold" in self.run_methods:
            os.chdir(self.params['alphafold_program_dir'])
            errormsg = ""
            uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
            if not os.path.exists(uniref90_sto):
                errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            colabfold_a3m = alndir + '/' + targetname + '_colabfold.a3m'
            if not os.path.exists(colabfold_a3m):
                errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {colabfold_a3m}\n"

            if len(errormsg) == 0:
                if not complete_result(f"{outdir}/colabfold", 5 * int(self.params['num_monomer_predictions_per_model'])):
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_path {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--custom_msa {colabfold_a3m} " \
                              f"--uniref90_sto {uniref90_sto} " \
                              f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--output_dir {outdir}/colabfold"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)
            else:
                print(errormsg)

        if "colabfold+seq_template" in self.run_methods:
            os.chdir(self.params['alphafold_program_dir'])
            errormsg = ""
            uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
            if not os.path.exists(uniref90_sto):
                errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            colabfold_a3m = alndir + '/' + targetname + '_colabfold.a3m'
            if not os.path.exists(colabfold_a3m):
                errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {colabfold_a3m}\n"

            if len(errormsg) == 0:
                if not complete_result(f"{outdir}/colab_seq_temp", 5 * int(self.params['num_monomer_predictions_per_model'])):
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_path {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--custom_msa {colabfold_a3m} " \
                              f"--uniref90_sto {uniref90_sto} " \
                              f"--temp_struct_csv {template_dir}/sequence_templates.csv " \
                              f"--struct_atom_dir {template_dir}/templates " \
                              f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--output_dir {outdir}/colab_seq_temp"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)
            else:
                print(errormsg)

        if "img" in self.run_methods:
            os.chdir(self.params['alphafold_program_dir'])
            errormsg = ""
            uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
            if not os.path.exists(uniref90_sto):
                errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            img_a3m = alndir + '/' + targetname + '.a3m'
            if not os.path.exists(img_a3m):
                errormsg = errormsg + f"Cannot find img alignment for {targetname}: {img_a3m}\n"

            if len(errormsg) == 0:
                if not complete_result(f"{outdir}/img", 5 * int(self.params['num_monomer_predictions_per_model'])):
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_path {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--custom_msa {img_a3m} " \
                              f"--uniref90_sto {uniref90_sto} " \
                              f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--output_dir {outdir}/img"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)
            else:
                print(errormsg)

        if "img+seq_template" in self.run_methods:
            os.chdir(self.params['alphafold_program_dir'])
            errormsg = ""
            uniref90_sto = alndir + '/' + targetname + '_uniref90.sto'
            if not os.path.exists(uniref90_sto):
                errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            img_a3m = alndir + '/' + targetname + '.a3m'
            if not os.path.exists(img_a3m):
                errormsg = errormsg + f"Cannot find img alignment for {targetname}: {img_a3m}\n"

            if len(errormsg) == 0:
                if not complete_result(f"{outdir}/img_seq_temp", 5 * int(self.params['num_monomer_predictions_per_model'])):
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_path {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--custom_msa {img_a3m} " \
                              f"--uniref90_sto {uniref90_sto} " \
                              f"--temp_struct_csv {template_dir}/sequence_templates.csv " \
                              f"--struct_atom_dir {template_dir}/templates " \
                              f"--monomer_num_ensemble {self.params['monomer_num_ensemble']} " \
                              f"--monomer_num_recycle {self.params['monomer_num_recycle']} " \
                              f"--num_monomer_predictions_per_model {self.params['num_monomer_predictions_per_model']} " \
                              f"--output_dir {outdir}/img_seq_temp"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)
            else:
                print(errormsg)


    def process(self, monomers, alndir, outdir, templatedir=None):
        outdir = os.path.abspath(outdir) + "/"
        for fasta_path in monomers:
            fasta_name = pathlib.Path(fasta_path).stem
            monomer_aln_dir = alndir + '/' + fasta_name
            monomer_outdir = outdir + '/' + fasta_name
            monomer_template_dir = ""
            if templatedir is not None:
                monomer_template_dir = templatedir + '/' + fasta_name
            self.process_single(fasta_path=fasta_path, alndir=monomer_aln_dir, outdir=monomer_outdir,
                                template_dir=monomer_template_dir)

        print("The tertiary structure generation for monomers has finished!")
