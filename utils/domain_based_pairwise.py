import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
import pickle
from scipy.stats import pearsonr
from bml_casp15.common.util import check_file, check_dir, check_dirs, makedir_if_not_exists, check_contents, \
    read_option_file, is_file, is_dir
from bml_casp15.common.protein import read_qa_txt_as_df


def get_sequence(contents):
    """Enclosing logic in a function to simplify code"""

    res_codes = [
        # 20 canonical amino acids
        ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
        ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
        ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
        ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
        ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
        # Non-canonical amino acids
        # ('MSE', 'M'), ('SOC', 'C'),
        # Canonical xNA
        ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
        ('  T', 'T'),
    ]

    three_to_one = dict(res_codes)
    # _records = set(['ATOM  ', 'HETATM'])
    _records = set(['ATOM  '])

    sequence = []
    read = set()
    for line in contents:
        line = line.strip()
        if line[0:6] in _records:
            resn = line[17:20]
            resi = line[22:26]
            icode = line[26]
            r_uid = (resn, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
            else:
                continue
            aa_resn = three_to_one.get(resn, 'X')
            sequence.append(aa_resn)

    return ''.join(sequence)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=str, required=True)
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--domain', type=str, required=True)
    parser.add_argument('--af_file', type=str, required=True)

    args = parser.parse_args()

    params = read_option_file(args.option_file)

    domain_indices = args.domain.split(',')

    domain_start = 1

    extract_script = "/home/bml_casp15/BML_CASP15/utils/extract_domain.pl"
    prev_df = None
    for i in range(len(domain_indices)):
        domain_end = int(domain_indices[i])
        domain_dir = args.outdir + '/' + str(domain_end)
        makedir_if_not_exists(domain_dir)

        domain_pdb_dir = domain_dir + '/pdb'
        makedir_if_not_exists(domain_pdb_dir)
        domain_fasta = ""
        for pdb in os.listdir(args.indir):
            os.system(f"perl {extract_script} {args.indir}/{pdb} {domain_pdb_dir}/{pdb} {domain_start} {domain_end}")
            if len(domain_fasta) == 0:
                domain_seq = get_sequence(open(f"{domain_pdb_dir}/{pdb}").readlines())
                fw = open(f"{domain_dir}/domain{i+1}.fasta", 'w')
                fw.write(f">domain{i+1}\n{domain_seq}")
                fw.close()
                domain_fasta = f"{domain_dir}/domain{i+1}.fasta"

        if not os.path.exists(domain_dir + '/pairwise_ranking.tm'):
            os.chdir(domain_dir)
            with open("model.list", 'w') as fw:
                for pdb in os.listdir(domain_pdb_dir):
                    fw.write(f"{domain_pdb_dir}/{pdb}\n")
            print(f"{params['qscore_program']} model.list {domain_fasta} {params['tmscore_program']} . pairwise_ranking")
            os.system(
                f"{params['qscore_program']} model.list {domain_fasta} {params['tmscore_program']} . pairwise_ranking")

        pairwise_ranking_df = read_qa_txt_as_df(domain_dir + '/pairwise_ranking.tm')
        pairwise_ranking_df = pairwise_ranking_df.add_suffix(str(i))
        pairwise_ranking_df['model'] = pairwise_ranking_df[f"model{i}"]
        pairwise_ranking_df = pairwise_ranking_df.drop([f"model{i}"], axis=1)

        if prev_df is None:
            prev_df = pairwise_ranking_df
        else:
            prev_df = prev_df.merge(pairwise_ranking_df, on='model')

        domain_start = domain_end + 1

    avg_scores = []
    for i in range(len(prev_df)):
        avg_score = 0
        for j in range(len(domain_indices)):
            avg_score += prev_df.loc[i, f'score{j}']
        avg_score = avg_score / (len(domain_indices))
        avg_scores += [avg_score]

    prev_df['pairwise_avg_score'] = avg_scores

    af_df = pd.read_csv(args.af_file)
    prev_df = prev_df.merge(af_df, on='model')

    avg_scores = []
    for i in range(len(prev_df)):
        avg_score = (prev_df.loc[i, 'pairwise_avg_score'] + prev_df.loc[i, 'plddt_avg']/100)/2
        avg_scores += [avg_score]

    prev_df['avg_score'] = avg_scores

    print(prev_df)
    prev_df.to_csv(args.outdir + "/summary.csv")
