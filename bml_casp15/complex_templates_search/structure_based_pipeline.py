import copy
import os
import sys
import time
from bml_casp15.common.util import makedir_if_not_exists, check_dirs
import pandas as pd
from multiprocessing import Pool

def loadFastaDictionary(dict_file):
    fasta_dict = {}
    with open(dict_file, "r") as f:
        for line in f:
            fasta_dict[line.strip().split(":")[0].strip()] = line.strip().split(":")[1].strip()
    return fasta_dict


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
                contents = item.split(",")
                aligned_lenght = contents[0].strip().split("=")[1]
                rmsd = contents[1].strip().split("=")[1]
                id = contents[2].strip().split("=")[2]
            if "(if normalized by length of Chain_1)" in line:
                line = line[0:line.index("(if normalized by length of Chain_1)")].strip()
                tmscore = line.split("=")[1].strip()
        aln_temp = contents[-3].rstrip()
        aln_query = contents[-1].rstrip()
        return [os.path.basename(gen_dir).split(".")[0].strip(), tmscore, id, aligned_lenght, aln_temp, aln_query]
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

    def search(self, monomers, outdir, is_homodimer = False):

        outdir = os.path.abspath(outdir) + "/"

        dimer_list = os.path.abspath(self.params["heterodimer_list"])
        if is_homodimer:
            dimer_list = os.path.abspath(self.params["homodimer_list"])

        temp_specific_dimers = read_pair_file_into_array(dimer_list)

        pdbs_to_compare = only_protein_names(temp_specific_dimers)

        makedir_if_not_exists(outdir)

        os.chdir(self.params['complex_atom_dir'])

        print(self.params['complex_atom_dir'])

        all_scores = {}

        for i in range(len(monomers)):

            monomer = monomers[i]

            os.system("cp " + monomer + " " + outdir + os.path.basename(monomer))

            if is_homodimer and i > 0:

                all_scores[i] = copy.deepcopy(temp_score)

            else:
                tm_score_dir = f'{outdir}/tmscore/{monomer}/'
                makedir_if_not_exists(tm_score_dir)

                aln_score_dir = f'{outdir}/aln/{monomer}/'
                makedir_if_not_exists(aln_score_dir)

                process_list = []
                for pdb in pdbs_to_compare:
                    if not os.path.exists(f"{tm_score_dir}/{pdb}.atom.txt"):
                        process_list.append([monomer, pdb + '.atom', self.params['tmalign_program'], tm_score_dir])

                pool = Pool(processes=10)
                results = pool.map(tmalign, process_list)
                pool.close()
                pool.join()

                txt_files = specific_filters(tm_score_dir)

                temp_score = []
                for proteins in txt_files:
                    val = tmscore_file_reader(proteins)
                    print(val)
                    len_val = len(val)
                    if len(val) > 0 and float(val[1]) > float(self.params["tmalign_threshold"]):
                        template_name = os.path.basename(proteins).split(".")[0]
                        # aln_out = aln_output_formatter(_query=val[len_val - 1], _query_name=current_pdb_name,
                        #                                _aligned_len=val[len_val - 3], _temp_name=template_name,
                        #                                _temp=val[len_val - 2])
                        temp_score.append(val)
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

            chain1_tmscore = all_scores[0]['chain1']
            chain2_tmscore = all_scores[1]['chain2']

            avg_score = (float(chain1_tmscore) + float(chain2_tmscore))/2

            temp_scores += [chain1, chain2, chain1_tmscore, chain2_tmscore, avg_score, min(chain1_tmscore, chain2_tmscore),
                            max(chain1_tmscore, chain2_tmscore)]

        temp_scores_sorted = sorted(temp_scores, key=lambda energy: float(energy[4]), reverse=True)

        df = pd.DataFrame(temp_scores_sorted, columns=['chain1', 'chain2', 'tmscore1', 'tmscore2', 'avg_tmscore', 'min_tmscore', 'max_tmscore'], dtype=float)

        df.to_csv(outdir + "/result_contact_" + str(self.params["tmalign_threshold"]) + ".csv")

