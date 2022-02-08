import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import pathlib




class Monomer_tertiary_structure_prediction_pipeline:

    def __init__(self, params, run_methods):

        self.params = params

        self.run_methods = run_methods

    def process(self, monomers, alndir, outdir):

        outdir = os.path.abspath(outdir) + "/"

        os.chdir(self.params['alphafold_program_dir'])

        for fasta_path in monomers:

            fasta_name = pathlib.Path(fasta_path).stem

            cmd = ""

            targetname = fasta_name

            monomer_aln_dir = alndir + '/' + targetname

            if "original" in self.run_methods:

                errormsg = ""

                if not os.path.exists(monomer_aln_dir):
                    errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {monomer_aln_dir}\n"

                uniclust30_a3m = monomer_aln_dir + '/' + targetname + '_uniclust30.a3m'
                if not os.path.exists(uniclust30_a3m):
                    errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {uniclust30_a3m}\n"

                bfd_a3m = monomer_aln_dir + '/' + targetname + '_bfd.a3m'
                if not os.path.exists(bfd_a3m):
                    errormsg = errormsg + f"Cannot find bfd alignment for {targetname}: {bfd_a3m}\n"

                mgnify_sto = monomer_aln_dir + '/' + targetname + '_mgnify.sto'
                if not os.path.exists(mgnify_sto):
                    errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

                uniref90_sto = monomer_aln_dir + '/' + targetname + '_uniref90.sto'
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                if len(errormsg) == 0:
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_paths {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--uniclust_a3ms {uniclust30_a3m} " \
                              f"--bfd_a3ms {bfd_a3m} " \
                              f"--mgnify_stos {mgnify_sto} " \
                              f"--uniref90_stos {uniref90_sto} " \
                              f"--output_dir {outdir}/original"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)

            if "rosettafold" in self.run_methods:
                errormsg = ""
                uniref90_sto = monomer_aln_dir + '/' + targetname + '_uniref90.sto'
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                rosettafold_a3m = monomer_aln_dir + '/' + targetname + '_rosettafold.a3m'
                if not os.path.exists(rosettafold_a3m):
                    errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {rosettafold_a3m}\n"

                if len(errormsg) == 0:
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_paths {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--custom_msas {rosettafold_a3m} " \
                              f"--uniref90_stos {uniref90_sto} " \
                              f"--output_dir {outdir}/rosettafold"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)

            if "colabfold" in self.run_methods:
                errormsg = ""
                uniref90_sto = monomer_aln_dir + '/' + targetname + '_uniref90.sto'
                if not os.path.exists(uniref90_sto):
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                colabfold_a3m = monomer_aln_dir + '/' + targetname + '_colabfold.a3m'
                if not os.path.exists(colabfold_a3m):
                    errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {colabfold_a3m}\n"

                if len(errormsg) == 0:
                    try:
                        cmd = f"python {self.params['alphafold_program']} --fasta_paths {fasta_path} " \
                              f"--env_dir {self.params['alphafold_env_dir']} " \
                              f"--database_dir {self.params['alphafold_database_dir']} " \
                              f"--custom_msas {colabfold_a3m} " \
                              f"--uniref90_stos {uniref90_sto} " \
                              f"--output_dir {outdir}/colabfold"
                        print(cmd)
                        os.system(cmd)
                    except Exception as e:
                        print(e)

        print("The tertiary structure generation for monomers has finished!")
