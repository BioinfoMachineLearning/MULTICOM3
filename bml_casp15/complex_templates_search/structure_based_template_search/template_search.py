import copy
import concurrent.futures
import os
import sys
import time

def loadFastaDictionary(dict_file):
    fasta_dict = {}
    with open(dict_file, "r") as f:
        for line in f:
            fasta_dict[line.strip().split(":")[0].strip()] = line.strip().split(":")[1].strip()
    return fasta_dict


def aln_output_formatter(_query, _temp, _aligned_len, _query_name, _temp_name):
    str_format = f'C; Aligned length= {_aligned_len}\n'
    str_format = str_format + f'>P1;{_temp_name}\n'
    str_format = str_format + f'Structure: {_temp_name}: 1  :  :{len(_temp.replace("-",""))} : : : : :\n'
    str_format = str_format + f'{_temp}*\n'
    str_format = str_format + f'>P1;{_query_name}\n'
    str_format = str_format + f'Structure: {_query_name}: 1    : :{len(_query.replace("-",""))} : : :  :\n'
    str_format = str_format + f'{_query}*'
    # print(str_format)
    return str_format

template_settings = loadFastaDictionary("template_search_settings.txt")

# INPUTS
ligand = os.path.abspath(sys.argv[1])
receptor = os.path.abspath(sys.argv[2])
is_homodimer = int(sys.argv[3])
db = sys.argv[4]
# interacting_files = int(sys.argv[4])
output_dir =os.path.abspath( sys.argv[5])+"/"

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

        f = open(gen_dir, "r")
        if f.mode == 'r':
            contents = f.read()
            f.close()

        for item in contents.split("\n"):
            temp_array = []
            if "Aligned length=" in item:
                aligned_lenght = item.strip().split(",")[0].strip().split("=")[1]
                rmsd = item.strip().split(",")[1].strip().split("=")[1]
                tmscore = item.strip().split(",")[2].strip().split("=")[1]
                id = item.strip().split(",")[3].strip().split("=")[1]
                aln_temp = contents.strip().split("\n")[-3]
                aln_query = contents.strip().split("\n")[-1]
        return [os.path.basename(gen_dir).split(".")[0].strip(), tmscore, id, aligned_lenght,aln_temp, aln_query]

def only_protein_names(_input):
    output = []
    for name in _input:
        temp_array = name.split(",")
        output.append(temp_array[0].strip())
        output.append(temp_array[1].strip())
    output = list(dict.fromkeys(output))
    return output


def tmalign(_name):
    path_b = "./" + _name + '.atom'
    name_a = _name + ".atom"
    name_b = os.path.basename(model)
    if os.path.isfile(path_b):
        # print(name_a + "," + name_b)
        cmd = template_settings["TM_ALIGN_DIR"]  + " " + path_b + " " + model+" -c"
        print(cmd + " > " + LOG_DIR + '/' + _name + ".txt" + "\n")
        counter_val = os.system(cmd + " > " + LOG_DIR + '/' + _name + ".txt" + "\n")

# temp_specific_dimers = read_pair_file_into_array(interacting_files)
if is_homodimer == 0:
    temp_specific_dimers = read_pair_file_into_array(os.path.abspath(template_settings["HETERO_FILE"]))
else:
    # print(os.path.abspath(template_settings["HOMO_FILE"]))
    temp_specific_dimers = read_pair_file_into_array(os.path.abspath(template_settings["HOMO_FILE"]))

filter_all = only_protein_names(temp_specific_dimers)
input_pdbs = []
LOG_DIR = ""

if not os.path.exists(output_dir):
    os.system("mkdir -p " + output_dir)
tm_score_a = output_dir + 'tmscore/a/'

if not os.path.exists(tm_score_a):
    os.system("mkdir -p " + tm_score_a)

aln_score_a = output_dir + 'aln/a/'
if not os.path.exists(aln_score_a):
    os.system("mkdir -p " + aln_score_a)


if is_homodimer == 0:
    tm_score_b = output_dir + 'tmscore/b/'
    if not os.path.exists(tm_score_b):
        os.system("mkdir -p " + tm_score_b)

    aln_score_b = output_dir + 'aln/b/'
    if not os.path.exists(aln_score_b):
        os.system("mkdir -p " + aln_score_b)

all_score_a = []
all_score_b = []
if is_homodimer == 0:
    contact_list_string = read_pair_file_into_string(os.path.abspath(template_settings["HETERO_FILE"]))
else:
    contact_list_string = read_pair_file_into_string(os.path.abspath(template_settings["HOMO_FILE"]))



# copyin input files
os.system("cp " + ligand + " " + output_dir + os.path.basename(ligand))
os.system("cp " + receptor + " " + output_dir + os.path.basename(receptor))
input_pdbs.append(output_dir + os.path.basename(ligand))
input_pdbs.append(output_dir + os.path.basename(receptor))

start = time.perf_counter()

LOG_DIR = tm_score_a
model = ligand
os.chdir(db)
# os.chdir(template_settings["ATOM_DIR"])
with concurrent.futures.ProcessPoolExecutor() as executor:
    results = executor.map(tmalign, filter_all)

# Just To make sure no overlapping of the model variable
time.sleep(int(template_settings["SLEEP_TIME"]))
if is_homodimer == 0:
    model = receptor
    LOG_DIR = tm_score_b
    with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
        results = executor.map(tmalign, filter_all)

for pdb in input_pdbs:
    txt_files = []
    current_aln_dir = ""
    current_pdb_name = os.path.basename(pdb).split(".")[0]
    if input_pdbs.index(pdb) == 0:
        txt_files = specific_filters(tm_score_a)
        current_aln_dir = aln_score_a
    else:
        if is_homodimer == 0:
            txt_files = specific_filters(tm_score_b)
            current_aln_dir = aln_score_b

    temp_score = []
    counter = 0
    for proteins in txt_files:
        if os.path.exists(proteins):
            val = tmscore_file_reader(proteins)
            len_val = len(val)
            if len(val) > 0 and float(val[1]) > float(template_settings["THRESHOLD"]):
                template_name = os.path.basename(proteins).split(".")[0]
                aln_out = aln_output_formatter(_query=val[len_val - 1], _query_name=current_pdb_name,
                                               _aligned_len=val[len_val - 3], _temp_name=template_name,
                                               _temp=val[len_val - 2])
                temp_score.append(val)
                write2file(current_aln_dir + template_name + ".pir", aln_out)
                # tmscore
                counter = counter + 1

    if input_pdbs.index(pdb) == 0:
        all_score_a = copy.deepcopy(temp_score)
    else:
        all_score_b = copy.deepcopy(temp_score)

    temp_score = []
    txt_files = []

# sorted
all_score_a = sorted(all_score_a, key=lambda energy: float(energy[1]), reverse=True)

if is_homodimer == 0:
    all_score_b = sorted(all_score_b, key=lambda energy: float(energy[1]), reverse=True)
else:
    all_score_b = copy.deepcopy(all_score_a)


output_array = []
for a in all_score_a:
    for b in all_score_b:
        if a[0][0:4] == b[0][0:4]:
            avg_v = (float(a[1]) + float(b[1])) / 2
            min_v = min(float(a[1]), float(b[1]))
            max_v = max(float(a[1]), float(b[1]))
            temp = a[0][0:4] + "," + a[0] + "," + b[0] + "," + str(a[1]) + "," + str(b[1]) + "," + str(
                round(avg_v, 3)) + "," + str(round(min_v, 3)) + "," + str(round(max_v, 3)) + "\n"
            output_array.append(temp)


output_string = list(dict.fromkeys(output_array))

final_array=[]
for val in output_string:
    temp = val.split(",")
    if temp[1] + "," + temp[2]  in contact_list_string or\
            temp[2] + "," + temp[1]  in contact_list_string:
        final_array.append(val)
# write2file(output_dir + "/result_contact.txt", value)

final_array = sorted(final_array, key=lambda mean_value: float(mean_value.split(",")[5]), reverse=True)[0:int(template_settings["INTERMEDIATE_TEMPLATE_COUNT"])]
final_output = "COMPLEX_ID,CHAIN_ID_1,CHAIN_ID_2,TMSCORE_1,TMSCORE_2,TMSCORE_AVERAGE,TMSCORE_MIN,TMSCORE_MAX"+"\n"
value = ""
for values in final_array:
    final_output=final_output+values
write2file(output_dir + "/result_contact_" + str(template_settings["THRESHOLD"]) + ".csv", final_output)
finish = time.perf_counter()

print(f'Finished in {round(finish - start, 2)} second(s)')