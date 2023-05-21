import re
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from multicom3.monomer_alignment_generation.alignment import *
from multicom3.complex_alignment_generation.pipeline import *


class Combine_alignment:

    def __init__(self, hhfilter=""):
        self.hhfilter = hhfilter

    def combine_alignment(self, chain1, chain2, alignments, outdir):
        combine_seqs = {'chain1': [], 'chain2': [], 'combine': [], 'index_1': [], 'index_2': []}
        aln1_start_index = 0
        aln2_start_index = 0
        for alignment in alignments:
            chain1_a3m = alignment['chain1']
            with open(chain1_a3m) as f:
                if chain1_a3m.find('sto') > 0:
                    aln1 = Alignment.from_file(f, format="stockholm")
                else:
                    aln1 = Alignment.from_file(f, format="a3m")
                combine_seqs['chain1'] += list(aln1.seqs)

            chain2_a3m = alignment['chain2']
            with open(chain2_a3m) as f:
                if chain2_a3m.find('sto') > 0:
                    aln2 = Alignment.from_file(f, format="stockholm")
                else:
                    aln2 = Alignment.from_file(f, format="a3m")
                combine_seqs['chain2'] += list(aln2.seqs)

            concatenate_a3m = alignment["concatenate_a3m"]
            with open(concatenate_a3m) as f:
                comb_aln = Alignment.from_file(f, format="a3m")
                combine_seqs['combine'] += list(comb_aln.seqs)

            interaction = pd.read_csv(alignment['interaction_csv'])
            combine_seqs['index_1'] += [str(i + aln1_start_index - 1) for i in interaction.index_1]
            combine_seqs['index_2'] += [str(i + aln2_start_index - 1) for i in interaction.index_2]
            aln1_start_index += len(aln1.seqs)
            aln2_start_index += len(aln2.seqs)

        with open(f"{outdir}/{chain1}.a3m", 'w') as f:
            f.write(f'>{chain1}\n{aln1.main_seq}\n')
            for i in range(len(combine_seqs['chain1'])):
                seq = combine_seqs['chain1'][i]
                f.write(f'>{i}\n{seq}\n')

        with open(f"{outdir}/{chain2}.a3m", 'w') as f:
            f.write(f'>{chain2}\n{aln2.main_seq}\n')
            for i in range(len(combine_seqs['chain2'])):
                seq = combine_seqs['chain2'][i]
                f.write(f'>{i}\n{seq}\n')

        with open(f"{outdir}/{chain1}.a3m") as f:
            aln_1 = Alignment.from_file(f, format="a3m")

        with open(f"{outdir}/{chain2}.a3m") as f:
            aln_2 = Alignment.from_file(f, format="a3m")

        pair_ids = pd.DataFrame({"id_1": combine_seqs['index_1'], "id_2": combine_seqs['index_2']})

        comb_dict = write_dimer_a3ms(pair_ids, aln_1, aln_2, '', outdir)

        with open(comb_dict['aln_file']) as f:
            comb_aln = Alignment.from_file(f, format="a3m")
            for seq in comb_aln.seqs:
                if seq not in combine_seqs['combine']:
                    print(f"Cannot find {seq} in combine alignment file: {comb_dict['aln_file']}")

        os.system(f"{self.hhfilter} -diff 50000 -i {comb_dict['aln_file']} -o {outdir}/{chain1}_{chain2}.90.a3m -id 90")

        filt_interact_dict = {'id_1': [], 'id_2': [], 'index_1': [], 'index_2': []}
        with open(f"{outdir}/{chain1}_{chain2}.90.a3m") as f:
            comb_filt_aln = Alignment.from_file(f, format="a3m")
            for id1, id2, index1, index2 in zip(comb_dict['pair_ids'].id_1, comb_dict['pair_ids'].id_2, comb_dict['pair_ids'].index_1, comb_dict['pair_ids'].index_2):
                comb_seq = aln_1[id1] + aln_2[id2]
                if comb_seq in comb_filt_aln.seqs:
                    filt_interact_dict['id_1'] += [id1]
                    filt_interact_dict['id_2'] += [id2]
                    filt_interact_dict['index_1'] += [index1]
                    filt_interact_dict['index_2'] += [index2]

        df = pd.DataFrame(filt_interact_dict)
        df.to_csv(f"{outdir}/{chain1}_{chain2}_interact.90.csv", index=False)