import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from bml_casp15.common.util import makedir_if_not_exists
import glob


class DPROQ:

    def __init__(self, dproq_program):
        self.dproq_program = dproq_program

    def run(self, indir, outdir):

        indir = os.path.abspath(indir)
        outdir = os.path.abspath(outdir)
        makedir_if_not_exists(outdir)

        workdir = outdir + '/work'
        resultdir = outdir + '/result'

        if not os.path.exists(f"{outdir}/DOCKQ_ranking.csv") or not os.path.exists(
                f"{outdir}/Escore_ranking.csv"):
            print("Start to run DRPOQ QA on multimer models")
            cmd = f"sh {self.dproq_program} {indir} {workdir} {resultdir}"
            print(cmd)
            try:
                os.system(cmd)
                os.system(f"cp {resultdir}/DOCKQ_ranking.csv {outdir}")
                os.system(f"cp {resultdir}/Escore_ranking.csv {outdir}")
            except Exception as e:
                print(e)
                return None, None

        return pd.read_csv(f"{outdir}/DOCKQ_ranking.csv"), pd.read_csv(f"{outdir}/Escore_ranking.csv")
