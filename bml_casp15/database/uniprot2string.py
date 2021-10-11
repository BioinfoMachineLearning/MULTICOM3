import os, sys, argparse, time
from bml_casp15.common.util import is_dir, is_file, die, makedir_if_not_exists, clean_dir
from multiprocessing import Pool
from tqdm import tqdm

class uniprot2string:

    def __init__(self,
                 string_links_file: str,
                 string_aliases_file: str,
                 string2uniprot_map: str,
                 uniprot2pdb_mapping_file: str):
        self.string_links_file = string_links_file,
        self.string_aliases_file = string_aliases_file
        self.string2uniprot_map = string2uniprot_map


    def update(self):

        t1 = time.time()

        string_aliases_uniprot = self.string_aliases_file + "_uniprot"

        print("greping...")

        os.system(f"grep BLAST_UniProt_AC {self.string_aliases_file} > {string_aliases_uniprot}")

        print("reading aliases file")
        string_uniprot_map = {}

        process_list = []

        # read through the data table line-by-line to improve speed
        num_lines = sum(1 for line in open(string_aliases_uniprot))
        with open(string_aliases_uniprot) as f:
            for line in tqdm(f, total=num_lines):
                string_id, uniprot_id, _ = line.rstrip('\n').split()
                if string_id in string_uniprot_map:
                    string_uniprot_map[string_id] += f",UniRef100_{uniprot_id}"
                else:
                    string_uniprot_map[string_id] = f"UniRef100_{uniprot_id}"

        print("converting all string ids")

        with open(self.string2uniprot_map, "w") as fw:
            with open(self.string_links_file) as f:
                for line in f:
                    # ena_data is formatted as 'genome1:cds1,genome2:cds2'
                    string_id1, string_id2, score = line.rstrip('\n').split()
                    if score == "combined_score":
                        continue
                    if int(score) < 500:
                        continue
                    if string_id1 in string_uniprot_map and string_id2 in string_uniprot_map:

                        uniprot_ids1 = string_uniprot_map[string_id1].split(',')
                        uniprot_ids2 = string_uniprot_map[string_id2].split(',')

                        for uniprotid1 in uniprot_ids1:

                            for uniprotid2 in uniprot_ids2:

                                fw .write(f"{uniprotid1}\t{uniprotid2}\t{score}\n")


        t2 = time.time()

        # print(f"Total pdb count: {total_count}\n Total match count: {match_count}")
        print("Uniprot ids to string ids mapping file generation is done. Total time:" + (t2 - t1).__str__())
