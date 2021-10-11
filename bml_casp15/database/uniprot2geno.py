import os, sys, argparse, time
from bml_casp15.common.util import is_dir, is_file, die, makedir_if_not_exists, clean_dir
from multiprocessing import Pool
from tqdm import tqdm


class uniprot2geno:

    def __init__(self,
                 con_pro_cds_dir: str,
                 uniprot_official_mapping_file: str,
                 uniprot_to_embl_table: str,
                 ena_genome_location_table: str):
        self.con_pro_cds_dir = con_pro_cds_dir,
        self.uniprot_official_mapping_file = uniprot_official_mapping_file
        self.uniprot_to_embl_table = uniprot_to_embl_table
        self.ena_genome_location_table = ena_genome_location_table

    def update(self):

        t1 = time.time()

        if os.path.exists(self.uniprot_to_embl_table):
            os.system("rm " + self.uniprot_to_embl_table)

        print("1. Processing uniprot_to_embl_table")
        num_lines = sum(1 for line in open(self.uniprot_official_mapping_file))
        with open(self.uniprot_official_mapping_file) as f:
            with open(self.uniprot_to_embl_table, "w") as fw:
                for line in tqdm(f, total=num_lines):
                    contents = line.rstrip('\n').split('\t')
                    if len(contents[16].lstrip()) != 0 and len(contents[17].lstrip()) != 0:
                        uniprot_ac = contents[0]
                        uniprot_id = contents[1]

                        EMBL = contents[16]
                        EMBL_CDS = contents[17]

                        embls = EMBL.split(';')
                        embl_cdss = EMBL_CDS.split(';')

                        cds = ""
                        for i in range(len(embls)):
                            embl = embls[i].lstrip().rstrip()
                            embl_cds = embl_cdss[i].lstrip().rstrip()
                            cds += f"{embl}:{embl_cds},"
                        cds = cds.rstrip(",")
                        fw.write(f"{uniprot_ac} {uniprot_id} {cds}\n")
                        # print(f"{uniprot_ac} {uniprot_id} {cds}")

        print("2. Processing ena_genome_location_table")
        con_pro_files = os.listdir(self.con_pro_cds_dir)
        with open(self.ena_genome_location_table, 'w') as fw:
            for con_pro_file in con_pro_files:
                print("Processing " + con_pro_file)
                num_lines = sum(1 for line in open(self.con_pro_cds_dir + '/' + con_pro_file))
                with open(self.con_pro_cds_dir + '/' + con_pro_file) as f:
                    id = ""
                    subid = ""
                    cds = ""
                    start = -1
                    end = -1
                    bjoin = False
                    for line in tqdm(f, total=num_lines):
                        contents = line.rstrip('\n').split()
                        if contents[0] == "ID":
                            id = contents[1][0:-1]
                            bjoin = False
                        if len(contents) > 1:
                            if contents[0] == "FT":
                                if contents[1] == "CDS" and len(contents) == 3:
                                    cds = contents[2]
                                    if cds.find('join') >= 0:
                                        bjoin = True
                                        # print(cds)
                                    else:
                                        if cds.find('complement') >= 0:
                                            cds = cds[cds.index('(') + 1:cds.index(')')]
                                        # print(id)
                                        subcontents = cds.split(':')
                                        subid = subcontents[0]
                                        start, end = subcontents[1].split('..')
                                        start = start.lstrip("<>").rstrip("<>")
                                        end = end.lstrip("<>").rstrip("<>")
                                        start = int(start) - 1
                                elif bjoin:
                                    if contents[1].find('codon_start') < 0:
                                        cds += contents[1]
                                    else:
                                        # print(cds)
                                        cds = cds[cds.rindex('(') + 1:cds.index(')')]

                                        cds_contents = cds.split(',')
                                        start_min = -1
                                        end_max = 0
                                        subid = ""
                                        for cds_content in cds_contents:
                                            subcontents = cds_content.split(':')
                                            if len(subid) > 0 and subid != subcontents[0]:
                                                print(f"The ids in join are different: {id}")
                                            subid = subcontents[0]
                                            start, end = subcontents[1].split('..')
                                            start = int(start.lstrip("<>").rstrip("<>"))
                                            end = int(end.lstrip("<>").rstrip("<>"))
                                            start -= 1
                                            if start < start_min or start_min == -1:
                                                start_min = start
                                            if end > end_max:
                                                end_max = end
                                        bjoin = False
                                        start = start_min
                                        end = end_max
                            if contents[1].find('db_xref="UniProtKB/TrEMBL:') >= 0:
                                embl_id = contents[1].split(':')[1].rstrip('"')
                                fw.write(f"{id}.1\t{subid}\t{embl_id}\t{str(start)}\t{str(end)}\n")

        t2 = time.time()

        print("Geno map generation is done. Total time:" + (t2 - t1).__str__())

        return self.uniprot_to_embl_table, self.ena_genome_location_table
