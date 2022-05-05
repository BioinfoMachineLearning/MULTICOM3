import re, os
from collections import OrderedDict, defaultdict
from copy import deepcopy
import numpy as np
import pandas as pd
from bml_casp15.monomer_alignment_generation.alignment import *
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file


def read_colabfold_a3m(infile):
    aln_indices = []
    headers = []
    starts = []
    ends = []
    alns = []
    index = 0
    monomer_id = ""
    for line in open(infile):
        if line.startswith('>'):
            contents = line.rstrip('\n').split()
            if len(contents) > 1:
                header = contents[0]
                tstart = contents[7]
                tend = contents[8]

                headers += [header]
                starts += [tstart]
                ends += [tend]
                aln_indices += [index]
                index += 1
            else:
                monomer_id = contents[0][1:]
                headers += [contents[0]]
                starts += [-1]
                ends += [-1]
                aln_indices += [index]
                index += 1
        else:
            alns += [line.rstrip('\n')]

    return pd.DataFrame(
        {'aln_index': aln_indices, 'header': headers, 'tstart': starts, 'tend': ends, 'aln': alns}), monomer_id


class colabfold_interact_v2:

    def get_interactions(species_interact_a3m, species_interact_alignment_files,
                         colabfold_alignment_files, outdir, method):
        outdir = outdir + '/' + method
        makedir_if_not_exists(outdir)
        print(outdir)
        monomer_ids = []
        prev_df = None
        for i, alignment_file in enumerate(colabfold_alignment_files):
            curr_df, monomer_id = read_colabfold_a3m(alignment_file)
            curr_df.to_csv(outdir + '/' + monomer_id + '_colab.csv')
            monomer_ids += [monomer_id]
            if prev_df is None:
                prev_df = curr_df
            else:
                prev_df = prev_df.merge(curr_df, on=f"header", suffixes=(f"_{i}", f"_{i + 1}"))

        print(prev_df)
        prev_df.to_csv(outdir + '/colabfold_interact.csv')
        keep_indices = []
        for i in range(len(prev_df)):
            overlap = False
            for j in range(len(colabfold_alignment_files)):
                tstart_1 = prev_df.loc[i, f'tstart_{j + 1}']
                tend_1 = prev_df.loc[i, f'tend_{j + 1}']
                for k in range(j+1, len(colabfold_alignment_files)):
                    tstart_2 = prev_df.loc[i, f'tstart_{k + 1}']
                    tend_2 = prev_df.loc[i, f'tend_{k + 1}']
                    if max(tstart_1, tstart_2) < min(tend_1, tend_2):
                        overlap = True
                        break
                if overlap:
                    break
            if not overlap:
                keep_indices += [i]

        monomer_alignments = {}
        for i in range(len(species_interact_alignment_files)):
            os.system(f"cp {species_interact_alignment_files[i]} {outdir}/{monomer_ids[i]}.spec.a3m")
            monomer_alignments[monomer_ids[i]] = {'headers': [], 'seqs': []}
            for line in open(f"{outdir}/{monomer_ids[i]}.spec.a3m"):
                line = line.rstrip('\n')
                if line.startswith('>'):
                    monomer_alignments[monomer_ids[i]]['headers'] += [line]
                else:
                    monomer_alignments[monomer_ids[i]]['seqs'] += [line]

        seen_complex_seqs = [line.rstrip('\n') for line in open(species_interact_a3m) if not line.startswith('>')]
        for i in keep_indices:
            complex_seq = ''.join([prev_df.loc[i, f'aln_{j + 1}'] for j in range(len(colabfold_alignment_files))])
            if complex_seq not in seen_complex_seqs:
                for j in range(len(colabfold_alignment_files)):
                    header = prev_df.loc[i, f'header']
                    tstart = prev_df.loc[i, f'tstart_{j + 1}']
                    tend = prev_df.loc[i, f'tend_{j + 1}']
                    monomer_alignments[monomer_ids[j]]['headers'] += [f"{header}_{tstart}_{tend}"]
                    monomer_alignments[monomer_ids[j]]['seqs'] += [prev_df.loc[i, f'aln_{j + 1}']]

        pair_ids = {}
        for monomer_idx, monomer_id in enumerate(monomer_ids):
            with open(outdir + '/' + monomer_id + '_con.a3m', 'w') as fw:
                for i in range(len(monomer_alignments[monomer_id]['headers'])):
                    fw.write(f"{monomer_alignments[monomer_id]['headers'][i]}\n"
                             f"{monomer_alignments[monomer_id]['seqs'][i]}\n")
            pair_ids[f'index_{monomer_idx+1}'] = [i for i in range(len(monomer_alignments[monomer_id]['headers']))]

        interact_df = pd.DataFrame(pair_ids)
        interact_df.to_csv(f"{outdir}/{method}_interact.csv", index=False)
        return interact_df
