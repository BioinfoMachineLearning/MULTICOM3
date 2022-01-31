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

        fasta_paths = monomers

        fasta_names = [pathlib.Path(fasta_path).stem for fasta_path in fasta_paths]

        uniclust_a3ms = []

        bfd_a3ms = []

        mgnify_stos = []

        uniref90_stos = []

        custom_msas = []

        errormsg = ""

        for fasta_name in fasta_names:

            targetname = fasta_name

            monomer_aln_dir = alndir + '/' + targetname

            if "original" in self.run_methods:

                if not os.path.exists(monomer_aln_dir):
                    errormsg = errormsg + f"Cannot find alignment directory for {targetname}: {monomer_aln_dir}\n"
                    continue

                uniclust30_a3m = monomer_aln_dir + '/' + targetname + '_uniclust30.a3m'
                if os.path.exists(uniclust30_a3m):
                    uniclust_a3ms += [uniclust30_a3m]
                else:
                    errormsg = errormsg + f"Cannot find uniclust30 alignment for {targetname}: {uniclust30_a3m}\n"

                bfd_a3m = monomer_aln_dir + '/' + targetname + '_bfd.a3m'
                if os.path.exists(bfd_a3m):
                    bfd_a3ms += [bfd_a3m]
                else:
                    errormsg = errormsg + f"Cannot find bfd alignment for {targetname}: {bfd_a3m}\n"

                mgnify_sto = monomer_aln_dir + '/' + targetname + '_mgnify.sto'
                if os.path.exists(mgnify_sto):
                    mgnify_stos += [mgnify_sto]
                else:
                    errormsg = errormsg + f"Cannot find mgnify alignment for {targetname}: {mgnify_sto}\n"

                uniref90_sto = monomer_aln_dir + '/' + targetname + '_uniref90.sto'
                if os.path.exists(uniref90_sto):
                    uniref90_stos += [uniref90_sto]
                else:
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

            if "rosettafold" in self.run_methods:

                uniref90_sto = monomer_aln_dir + '/' + targetname + '_uniref90.sto'
                if os.path.exists(uniref90_sto):
                    uniref90_stos += [uniref90_sto]
                else:
                    errormsg = errormsg + f"Cannot find uniref90 alignment for {targetname}: {uniref90_sto}\n"

                rosettafold_a3m = monomer_aln_dir + '/' + targetname + '_rosettafold.a3m'
                if os.path.exists(rosettafold_a3m):
                    custom_msas += [rosettafold_a3m]
                else:
                    errormsg = errormsg + f"Cannot find rosettafold alignment for {targetname}: {rosettafold_a3m}\n"

        if len(errormsg) > 0:
            #raise ValueError(errormsg)
            print(errormsg)

        # cmd = f"sh {alphafold_program} {','.join(fasta_paths)} {self.params['alphafold_env_dir']} " \
        #       f"{self.params['alphafold_database_dir']} {','.join(uniclust_a3ms)} {','.join(bfd_a3ms)} " \
        #       f"{','.join(mgnify_stos)} {','.join(uniref90_stos)} {outdir}"

        os.chdir(self.params['alphafold_program_dir'])

        if "original" in self.run_methods:
            cmd = f"python {self.params['alphafold_program']} --fasta_paths {','.join(fasta_paths)} " \
                  f"--env_dir {self.params['alphafold_env_dir']} " \
                  f"--database_dir {self.params['alphafold_database_dir']} " \
                  f"--uniclust_a3ms {','.join(uniclust_a3ms)} " \
                  f"--bfd_a3ms {','.join(bfd_a3ms)} " \
                  f"--mgnify_stos {','.join(mgnify_stos)} " \
                  f"--uniref90_stos {','.join(uniref90_stos)} " \
                  f"--output_dir {outdir}/original"
            print(cmd)
            os.system(cmd)

        if "rosettafold" in self.run_methods:
            cmd = f"python {self.params['alphafold_program']} --fasta_paths {','.join(fasta_paths)} " \
                  f"--env_dir {self.params['alphafold_env_dir']} " \
                  f"--database_dir {self.params['alphafold_database_dir']} " \
                  f"--custom_msas {','.join(custom_msas)} " \
                  f"--uniref90_stos {','.join(uniref90_stos)} " \
                  f"--output_dir {outdir}/rosettafold"
            print(cmd)
            os.system(cmd)

        print("The structure based template searching for dimers has finished!")
