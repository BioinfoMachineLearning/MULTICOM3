import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses
from bml_casp15.tool.foldseek import *
import pickle
import numpy as np

def _convert_taln_seq_to_a3m(query_non_gaps, aln):
    for is_query_res_non_gap, sequence_res in zip(query_non_gaps, aln):
        if is_query_res_non_gap:
            yield sequence_res


class Monomer_iterative_generation_pipeline:

    def __init__(self, params):

        self.params = params

        self.max_iteration = 5

    def search_templates(self, inpdb, outdir):
        makedir_if_not_exists(outdir)
        foldseek_program = self.params['foldseek_program']
        foldseek_pdb_database = self.params['foldseek_pdb_database']
        foldseek_af_database = self.params['foldseek_af_database']
        foldseek_runner = Foldseek(binary_path=foldseek_program, databases=[foldseek_pdb_database, foldseek_af_database])
        return foldseek_runner.query(pdb=inpdb, outdir=outdir)

    def generate_msa_from_templates(self, fasta_file, template_file, outfile):
        targetname = None
        seq = None
        for line in open(fasta_file):
            line = line.rstrip('\n')
            if line.startswith('>'):
                targetname = line[1:]
            else:
                seq = line

        templates = pd.read_csv(template_file, sep='\t')

        alignments = {targetname: seq}
        seen_seq = []
        for i in range(len(templates)):
            target = templates.loc[i, 'target']
            qaln = templates.loc[i, 'qaln']
            qstart = int(templates.loc[i, 'qstart'])
            qend = int(templates.loc[i, 'qend'])
            taln = templates.loc[i, 'taln']
            tstart = templates.loc[i, 'tstart']
            tend = templates.loc[i, 'tend']

            query_non_gaps = [res != '-' for res in qaln]
            out_sequence = ''.join(_convert_taln_seq_to_a3m(query_non_gaps, taln))

            aln_full = ['-'] * len(seq)
            aln_full[qstart - 1:qend] = out_sequence
            taln_full_seq = ''.join(aln_full)
            if taln_full_seq in seen_seq:
                continue
            alignments[target] = taln_full_seq
            seen_seq += [taln_full_seq]

        fasta_chunks = (f">{k}\n{alignments[k]}" for k in alignments)

        with open(outfile + '.work', 'w') as fw:
            fw.write('\n'.join(fasta_chunks) + '\n')

        hhfilter_program = self.params['hhfilter_program']

        cmd = f"{hhfilter_program} -diff 1024 -i {outfile}.work -o {outfile}"

        os.system(cmd)


    def search(self, fasta_file, input_pdb_dir, outdir):

        input_pdb_dir = os.path.abspath(input_pdb_dir)

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        for i in range(0, 5):
            model_outdir = f"{outdir}/ranked_{i}"
            makedir_if_not_exists(model_outdir)

            os.system(f"cp {input_pdb_dir}/ranked_{i}.pdb {model_outdir}/start.pdb")
            os.system(f"cp {input_pdb_dir}/result_model_{i + 1}.pkl {model_outdir}/start.pkl")

            start_avg_lddt = 0
            with open(f"{model_outdir}/start.pkl", 'rb') as f:
                prediction_result = pickle.load(f)
                start_avg_lddt = np.mean(prediction_result['plddt'])

            start_pdb = f"{model_outdir}/start.pdb"
            for num_iteration in range(self.max_iteration):
                workdir = f"{model_outdir}/iteration{num_iteration + 1}"
                makedir_if_not_exists(workdir)
                foldseek_res = self.search_templates(inpdb=start_pdb, outdir=workdir + '/foldseek')
                os.system(f"cp {foldseek_res} {workdir}/structure_templates.csv")
                self.generate_msa_from_templates(fasta_file, template_file=f"{workdir}/structure_templates.csv", outfile=f"{workdir}/template.a3m")