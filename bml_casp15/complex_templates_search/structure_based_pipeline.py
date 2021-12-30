import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool


def read_pair_file_into_array(_input):
    output_array = []
    file = open(_input, "r")
    if file.mode == 'r':
        output_array = file.read().splitlines()
        file.close()
    return output_array


def read_pair_file_into_string(_input):
    output_array = []
    file = open(_input, "r")
    if file.mode == 'r':
        output_array = file.read()
        file.close()
    return output_array


def specific_filters(_input):
    os.chdir(_input)
    fileNames = []
    # fileNames = glob.glob(fasta_dir + "/*fasta")
    for root, directories, files in os.walk(_input):
        for file in files:
            if ".txt" in file:
                fileNames.append(os.path.abspath(file))
    return fileNames


def write2file(file, contents):
    with open(file, "w") as f:
        f.writelines(contents)


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
        return {'template': os.path.basename(gen_dir).split(".")[0].strip(),
                'tmscore': tmscore,
                'seqid': id,
                'aligned_length': aligned_length,
                'aln_temp': aln_temp,
                'tstart': 1,
                'tend': len(aln_temp.replace("-", "")),
                'aln_query': aln_query,
                'qstart': 1,
                'qend': len(aln_query.replace("-", ""))}
    return None


def only_protein_names(_input):
    output = []
    for name in _input:
        temp_array = name.split()
        output.append(temp_array[1].strip())
        output.append(temp_array[2].strip())
    output = list(dict.fromkeys(output))
    return output


def tmalign(inparams):
    model1, model2, tmalign_program, logdir = inparams

    modelname = os.path.basename(model2)

    cmd = f"{tmalign_program} {model1} {model2} -c > {logdir}/{modelname}.txt"

    print(cmd)

    os.system(cmd)

    return f"{logdir}/{modelname}.txt"


def aln_output_formatter(_query, _temp, _aligned_len, _query_name, _temp_name):
    str_format = f'C; Aligned length= {_aligned_len}\n'
    str_format = str_format + f'>P1;{_temp_name}\n'
    str_format = str_format + f'Structure: {_temp_name}: 1  :  :{len(_temp.replace("-", ""))} : : : : :\n'
    str_format = str_format + f'{_temp}*\n'
    str_format = str_format + f'>P1;{_query_name}\n'
    str_format = str_format + f'Structure: {_query_name}: 1    : :{len(_query.replace("-", ""))} : : :  :\n'
    str_format = str_format + f'{_query}*'
    # print(str_format)
    return str_format


class Complex_structure_based_template_search_pipeline:

    def __init__(self, params):

        self.params = params

    def search(self, monomers, outdir, is_homodimer=False):

        outdir = os.path.abspath(outdir) + "/"

        dimer_list = os.path.abspath(self.params["heterodimer_list"])
        if is_homodimer:
            dimer_list = os.path.abspath(self.params["homodimer_list"])

        temp_specific_dimers = read_pair_file_into_array(dimer_list)

        pdbs_to_compare = only_protein_names(temp_specific_dimers)

        makedir_if_not_exists(outdir)

        all_scores = {}

        for i in range(len(monomers)):

            os.chdir(self.params['complex_atom_dir'])

            monomer = monomers[i]

            os.system("cp " + monomer + " " + outdir + os.path.basename(monomer))

            if is_homodimer and i > 0:

                all_scores[i] = copy.deepcopy(temp_score)

            else:
                tm_score_dir = f'{outdir}/tmscore/{os.path.basename(monomer)}/'
                makedir_if_not_exists(tm_score_dir)

                aln_score_dir = f'{outdir}/aln/{os.path.basename(monomer)}/'
                makedir_if_not_exists(aln_score_dir)

                process_list = []
                for pdb in pdbs_to_compare:
                    if not os.path.exists(f"{tm_score_dir}/{pdb}.atom.txt") or len(
                            open(f"{tm_score_dir}/{pdb}.atom.txt").readlines()) == 0:
                        process_list.append([monomer, pdb + '.atom', self.params['tmalign_program'], tm_score_dir])

                pool = Pool(processes=10)
                results = pool.map(tmalign, process_list)
                pool.close()
                pool.join()

                txt_files = specific_filters(tm_score_dir)

                temp_score = {}
                for proteins in txt_files:
                    val = tmscore_file_reader(proteins)
                    if len(val) > 0:# and float(val[1]) > float(self.params["tmalign_threshold"]):
                        # aln_out = aln_output_formatter(_query=val[len_val - 1], _query_name=current_pdb_name,
                        #                                _aligned_len=val[len_val - 3], _temp_name=template_name,
                        #                                _temp=val[len_val - 2])
                        temp_score[val['template']] = val
                        # write2file(current_aln_dir + template_name + ".pir", aln_out)

                all_scores[i] = temp_score

        dimer_to_be_compare = []
        with open(dimer_list) as f:
            for line in f:
                line = line.rstrip()
                pdbcode, chain1, chain2 = line.split()
                dimer_to_be_compare += [f"{chain1}_{chain2}"]
                dimer_to_be_compare += [f"{chain2}_{chain1}"]

        temp_scores = []
        for dimer in dimer_to_be_compare:
            chain1, chain2 = dimer.split('_')

            chain1_tmscore = all_scores[0][chain1]['tmscore']
            chain2_tmscore = all_scores[1][chain2]['tmscore']
            avg_score = (float(chain1_tmscore) + float(chain2_tmscore)) / 2

            aln_temp1 = all_scores[0][chain1]['aln_temp']
            tstart1 = all_scores[0][chain1]['tstart']
            tend1 = all_scores[0][chain1]['tend']

            aln_query1 = all_scores[0][chain1]['aln_query']
            qstart1 = all_scores[0][chain1]['qstart']
            qend1 = all_scores[0][chain1]['qend']

            aln_temp2 = all_scores[1][chain2]['aln_temp']
            tstart2 = all_scores[1][chain2]['tstart']
            tend2 = all_scores[1][chain2]['tend']

            aln_query2 = all_scores[1][chain2]['aln_query']
            qstart2 = all_scores[1][chain2]['qstart']
            qend2 = all_scores[1][chain2]['qend']

            temp_scores.append([chain1, chain2, avg_score,
                               chain1_tmscore, aln_temp1, tstart1, tend1, aln_query1, qstart1, qend1,
                               chain2_tmscore, aln_temp2, tstart2, tend2, aln_query2, qstart2, qend2])

        temp_scores_sorted = sorted(temp_scores, key=lambda energy: float(energy[2]), reverse=True)

        df = pd.DataFrame(temp_scores_sorted[0:int(self.params['template_count'])],
                          columns=['chain1', 'chain2', 'avg_tmscore',
                                   'tmscore1', 'aln_temp1', 'tstart1', 'tend1', 'aln_query1', 'qstart1', 'qend1',
                                   'tmscore2', 'aln_temp2', 'tstart2', 'tend2', 'aln_query2', 'qstart2', 'qend2'])

        df.to_csv(outdir + "/structure_based_templates.csv", index=False)
