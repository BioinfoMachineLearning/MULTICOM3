import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool
import dataclasses


@dataclasses.dataclass(frozen=True)
class tmalign_result:
    template: str
    tpdbcode: str
    tmscore: float
    seqid: float
    aligned_length: int
    aln_temp: str
    tstart: int
    tend: int
    aln_query: str
    qstart: int
    qend: int


def tmscore_file_reader(gen_dir):
    if os.path.exists(gen_dir):
        contents = open(gen_dir).readlines()
        for line in contents:
            line = line.rstrip()
            if "Aligned length=" in line:
                line_contents = line.split(",")
                aligned_length = line_contents[0].split("=")[1].strip()
                rmsd = line_contents[1].strip().split("=")[1]
                id = line_contents[2].strip().split("=")[2].strip()
            if "(if normalized by length of Chain_1)" in line:
                line = line[0:line.index("(if normalized by length of Chain_1)")].strip()
                tmscore = line.split("=")[1].strip()
        aln_temp = contents[-2].rstrip()
        aln_query = contents[-4].rstrip()
        return tmalign_result(template=os.path.basename(gen_dir).split(".")[0].strip(),
                              tpdbcode=os.path.basename(gen_dir).split(".")[0].strip()[0:4],
                              tmscore=float(tmscore),
                              seqid=float(id),
                              aligned_length=int(aligned_length),
                              aln_temp=aln_temp,
                              tstart=1,
                              tend=len(aln_temp.replace("-", "")),
                              aln_query=aln_query,
                              qstart=1,
                              qend=len(aln_query.replace("-", "")))
    return None


def tmalign(inparams):
    model1, model2, tmalign_program, logdir = inparams
    modelname = os.path.basename(model2)
    cmd = f"{tmalign_program} {model1} {model2} -c > {logdir}/{modelname}.txt"
    # print(cmd)
    os.system(cmd)
    return f"{logdir}/{modelname}.txt"


def create_df(protein_results):
    row_list = []
    for name in protein_results:
        row_dict = dict(template=protein_results[name].template,
                        tpdbcode=protein_results[name].tpdbcode,
                        tmscore=protein_results[name].tmscore,
                        seqid=protein_results[name].seqid,
                        aligned_length=protein_results[name].aligned_length,
                        aln_temp=protein_results[name].aln_temp,
                        tstart=protein_results[name].tstart,
                        tend=protein_results[name].tend,
                        aln_query=protein_results[name].aln_query,
                        qstart=protein_results[name].qstart,
                        qend=protein_results[name].qend)
        row_list += [row_dict]
    return pd.DataFrame(row_list)


class Monomer_iteractive_template_search_pipeline:

    def __init__(self, params):

        self.params = params

    def search(self, monomer_pdb, otudir):

        outdir = os.path.abspath(outdir) + "/"

        makedir_if_not_exists(outdir)

        all_scores = {}

        monomers_paths = []
        for i in range(len(monomers)):
            monomer = monomers[i]
            os.system("cp " + monomer + " " + outdir + os.path.basename(monomer))
            monomers_paths += [outdir + os.path.basename(monomer)]

        seen_pdbcodes = set()

        os.chdir(self.params['template_atom_dir'])

        all_pdb_files = os.listdir(self.params['template_atom_dir'])

        for i in range(len(monomers_paths)):

            monomer = monomers_paths[i]

            print(f"Searching templates for {monomer}")

            tm_score_dir = f'{outdir}/tmscore/{os.path.basename(monomer)}/'

            makedir_if_not_exists(tm_score_dir)

            if len(seen_pdbcodes) == 0:
                pdbs_to_compare = all_pdb_files
            else:
                pdbs_to_compare = [pdb for pdb in all_pdb_files if pdb[0:4] in seen_pdbcodes]

            process_list = []
            for pdb in pdbs_to_compare:
                if not os.path.exists(f"{tm_score_dir}/{pdb}.txt") or len(
                        open(f"{tm_score_dir}/{pdb}.txt").readlines()) == 0:
                    process_list.append([monomer, pdb, self.params['tmalign_program'], tm_score_dir])

            pool = Pool(processes=30)
            results = pool.map(tmalign, process_list)
            pool.close()
            pool.join()

            temp_score = {}
            seen_pdbcodes = set()
            for protein in os.listdir(tm_score_dir):
                protein_result = tmscore_file_reader(tm_score_dir + '/' + protein)
                if len(protein_result.template) > 0:
                    temp_score[protein_result.template] = protein_result
                    seen_pdbcodes.add(protein_result.tpdbcode)
            # os.system("rm -rf " + tm_score_dir)
            all_scores[i] = temp_score

            print(temp_score)

        print(f"Extracting complex templates")

        prev_df = create_df(all_scores[0])
        print(prev_df)
        for i in range(1, len(all_scores)):
            curr_df = create_df(all_scores[i])
            prev_df = prev_df.merge(curr_df, how="inner", on='tpdbcode', suffixes=(str(i), str(i + 1)))

        avg_tmscores = []
        for i in range(len(prev_df)):
            avg_tmscore = 0
            for j in range(len(monomers)):
                avg_tmscore += float(prev_df.loc[i, f"tmscore{j + 1}"])
            avg_tmscore = avg_tmscore / len(monomers)
            avg_tmscores += [avg_tmscore]

        prev_df['avg_tmscore'] = avg_tmscores

        prev_df_sorted = prev_df.sort_values(by='avg_tmscore', ascending=False)

        prev_df_sorted.head(int(self.params['template_count'])).to_csv(outdir + "/structure_templates.csv", index=False)

        print("The structure based template searching for dimers has finished!")

        os.system(f"rm -rf {outdir}/tmscore")
