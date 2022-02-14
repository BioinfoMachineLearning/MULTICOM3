import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from bml_casp15.common.util import makedir_if_not_exists
import glob

class En_qa:

    def __init__(self, program_python, program_path, program_script, use_gpu):
        self.program_python = program_python
        self.program_path = program_path
        self.program_script = program_script
        self.use_gpu = use_gpu

    def run(self, input_dir, alphafold_prediction_dir, outputdir):

        input_dir = os.path.abspath(input_dir)
        alphafold_prediction_dir = os.path.abspath(alphafold_prediction_dir)
        outputdir = os.path.abspath(outputdir)

        makedir_if_not_exists(outputdir)

        os.system(f"cp -r {alphafold_prediction_dir} {outputdir}/alphafold_prediction")

        if not glob.glob(f"{outputdir}/alphafold_prediction/relaxed_model_*.pdb"):
            if not glob.glob(f"{outputdir}/alphafold_prediction/unrelaxed_model_*.pdb"):
                raise Exception(f"Cannot find relaxed models in {alphafold_prediction_dir}")
            else:
                os.chdir(f"{outputdir}/alphafold_prediction")
                for i in range(1, 6):
                    os.system(f"ln -s unrelaxed_model_{i}.pdb relaxed_model_{i}.pdb")

        os.chdir(self.program_path)

        print("Try to run the ensemble model first")
        cmd = f"{self.program_python} {self.program_script} --input {input_dir} --output {outputdir} " \
              f"--method ensemble --alphafold_prediction {outputdir}/alphafold_prediction"
        if not self.use_gpu:
            cmd += " --cpu"
            
        print(cmd)
        try:
            os.system(cmd)
        except Exception as e:
            print(e)

        if not os.path.exists(f"{outputdir}/result.txt"):
            print("Try to run the EGNN_esto9 model")
            cmd = f"{self.program_python} {self.program_script} --input {input_dir} --output {outputdir} " \
                  f"--method EGNN_esto9 --alphafold_prediction {outputdir}/alphafold_prediction"
            if not self.use_gpu:
                cmd += " --cpu"
            print(cmd)
            try:
                os.system(cmd)
            except Exception as e:
                print(e)

        models = []
        scores = []
        if os.path.exists(f"{outputdir}/result.txt"):
            for line in open(f"{outputdir}/result.txt"):
                line = line.rstrip('\n')
                contents = line.split()
                if contents[0] == "PFRMAT" or contents[0] == "TARGET" or contents[0] == "MODEL" or contents[
                    0] == "QMODE":
                    continue
                model, score = line.split()
                models += [model]
                scores += [score]
        else:
            pdbs = os.listdir(input_dir)
            for i in range(len(pdbs)):
                models += [pdbs[i]]
                scores += [0.0]

        return pd.DataFrame({'model': models, 'score': scores})
