import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
from bml_casp15.common.util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from bml_casp15.monomer_alignment_generation.alignment import *
from bml_casp15.tool import hhblits
from bml_casp15.tool import jackhmmer


def run_msa_tool(inparams):
    msa_runner, input_fasta_path, msa_out_path, msa_key = inparams
    """Runs an MSA tool, checking if output already exists first."""
    if not os.path.exists(msa_out_path) or len(open(msa_out_path).readlines()) == 0:
        result = msa_runner.query(input_fasta_path, msa_out_path)
    return msa_key, msa_out_path


class Monomer_alignment_generation_pipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(self,
                 jackhmmer_binary_path,
                 hhblits_binary_path,
                 uniref90_database_path,
                 mgnify_database_path,
                 small_bfd_database_path,
                 bfd_database_path,
                 uniref30_database_path,
                 uniclust30_database_path,
                 uniprot_database_path,
                 mgnify_max_hits: int = 501,
                 uniref_max_hits: int = 10000,
                 use_precomputed_msas: bool = False):
        """Initializes the data pipeline."""

        # alignment generation pipeline from alphafold

        self.jackhmmer_uniref90_runner = None
        self.hhblits_bfd_runner = None
        self.hhblits_uniref_runner = None
        self.jackhmmer_mgnify_runner = None
        self.hhblits_uniclust_runner = None
        self.jackhmmer_uniprot_runner = None
        self.hhblits_uniclust_folddock_runner = None
        self.jackhmmer_small_bfd_runner = None

        if len(uniref90_database_path) > 0:
            self.jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=uniref90_database_path,
                get_tblout=True)

        if len(bfd_database_path) > 0:
            self.hhblits_bfd_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[bfd_database_path])

        if len(uniref30_database_path) > 0:
            self.hhblits_uniref_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[uniref30_database_path])

        if len(mgnify_database_path) > 0:
            self.jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=mgnify_database_path,
                get_tblout=True)

        if len(uniclust30_database_path) > 0:
            self.hhblits_uniclust_runner = hhblits.HHBlits(
                binary_path=hhblits_binary_path,
                databases=[uniclust30_database_path])

        if len(uniprot_database_path) > 0:
            self.jackhmmer_uniprot_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=uniprot_database_path,
                get_tblout=True)

        # if os.path.exists(uniclust30_database_path):
        #     self.hhblits_uniclust_folddock_runner = hhblits.HHBlits(
        #         binary_path=hhblits_binary_path,
        #         databases=[uniclust30_database_path],
        #         all_seqs=True)

        if len(small_bfd_database_path) > 0:
            self.jackhmmer_small_bfd_runner = jackhmmer.Jackhmmer(
                binary_path=jackhmmer_binary_path,
                database_path=small_bfd_database_path,
                get_tblout=True)


    def process(self, input_fasta_path, msa_output_dir):
        """Runs alignment tools on the input sequence and creates features."""

        input_id = open(input_fasta_path).readlines()[0].rstrip('\n').lstrip('>')

        msa_process_list = []

        if self.jackhmmer_uniref90_runner is not None:
            uniref90_out_path = os.path.join(msa_output_dir, f'{input_id}_uniref90.sto')
            msa_process_list.append([self.jackhmmer_uniref90_runner, input_fasta_path, uniref90_out_path, 'uniref90_sto'])

        if self.jackhmmer_mgnify_runner is not None:
            mgnify_out_path = os.path.join(msa_output_dir, f'{input_id}_mgnify.sto')
            msa_process_list.append([self.jackhmmer_mgnify_runner, input_fasta_path, mgnify_out_path, 'mgnify_sto'])

        if self.jackhmmer_small_bfd_runner is not None:
            small_bfd_out_path = os.path.join(msa_output_dir, f'{input_id}_smallbfd.sto')
            msa_process_list.append([self.jackhmmer_small_bfd_runner, input_fasta_path, small_bfd_out_path, 'smallbfd_sto'])

        if self.hhblits_bfd_runner is not None:
            bfd_out_path = os.path.join(msa_output_dir, f'{input_id}_bfd.a3m')
            msa_process_list.append([self.hhblits_bfd_runner, input_fasta_path, bfd_out_path, 'bfd_a3m'])

        if self.hhblits_uniref_runner is not None:
            uniref30_out_path = os.path.join(msa_output_dir, f'{input_id}_uniref30.a3m')
            msa_process_list.append([self.hhblits_uniref_runner, input_fasta_path, uniref30_out_path, 'uniref30_a3m'])

        if self.hhblits_uniclust_runner is not None:
            uniclust30_out_path = os.path.join(msa_output_dir, f'{input_id}_uniclust30.a3m')
            msa_process_list.append([self.hhblits_uniclust_runner, input_fasta_path, uniclust30_out_path, 'uniclust30_a3m'])

        if self.hhblits_uniclust_folddock_runner is not None:
            uniclust30_all_out_path = os.path.join(msa_output_dir, f'{input_id}_uniclust30_all.a3m')
            msa_process_list.append([self.hhblits_uniclust_folddock_runner, input_fasta_path, uniclust30_all_out_path, 'uniclust30_all_a3m'])

        if self.jackhmmer_uniprot_runner is not None:
            uniprot_out_path = os.path.join(msa_output_dir, f'{input_id}_uniprot.sto')
            msa_process_list.append([self.jackhmmer_uniprot_runner, input_fasta_path, uniprot_out_path, 'uniprot_sto'])

        # pool = Pool(processes=len(msa_process_list))
        # results = pool.map(run_msa_tool, msa_process_list)
        # pool.close()
        # pool.join()
        #
        # result_dict = {}
        # for result in results:
        #     msa_key, msa_out_path = result
        #     if os.path.exists(msa_out_path):
        #         result_dict[msa_key] = msa_out_path

        result_dict = {}
        for msa_process_params in msa_process_list:
            msa_key, msa_out_path = run_msa_tool(msa_process_params)
            if os.path.exists(msa_out_path):
                result_dict[msa_key] = msa_out_path

        return result_dict
