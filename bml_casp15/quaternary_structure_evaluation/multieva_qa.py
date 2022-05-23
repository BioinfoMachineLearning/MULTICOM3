import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from bml_casp15.common.util import makedir_if_not_exists
import glob
import pathlib

class MultiEva_qa:

    def __init__(self, multieva_program):
        self.multieva_program = multieva_program

    def run(self, chain_id_map, fasta_path, input_dir, stoichiometry,
            alphafold_prediction_dir, outputdir):

        makedir_if_not_exists(outputdir)

        cmd = f"python {self.multieva_program} {fasta_path} {input_dir} {stoichiometry} {alphafold_prediction_dir} 20 {outputdir}"
            
        print(cmd)
        try:
            os.system(cmd)
        except Exception as e:
            print(e)

        result_csv = f"{outputdir}/{pathlib.Path(fasta_path).stem}.csv"
        print(result_csv)
        if not os.path.exists(result_csv):
            return None

        df = pd.read_csv(result_csv)
        if 'MMalign score' in df:
            df = df.sort_values(by=['MMalign score'], ascending=False)
        else:
            df = df.sort_values(by=['Final_Rank'])
        df.reset_index(inplace=True,drop=True)
        df.to_csv(result_csv)
        return result_csv