import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool


class Complex_sequence_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

        self.template_searcher = hmmsearch.Hmmsearch(
            binary_path=params['hmmsearch_program'],
            hmmbuild_binary_path=params['hmmbuild_program'],
            database_path=params['dimer_database'])

    def search(self, monomers, outdir, is_homodimer=False):
        return None
