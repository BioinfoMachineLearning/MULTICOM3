import os, sys, argparse, time
from bml_casp15.common.util import is_dir, is_file, die, makedir_if_not_exists, clean_dir
from multiprocessing import Pool
from tqdm import tqdm


def compare_seq(inparams):
    seq_src_file, seq_dest_path, seq_dest_files, uniprot_id = inparams
    # print(f"Processing {uniprot_id}")
    contents = open(seq_src_file, "r").readlines()
    src_pdbcode = contents[0].rstrip('\n')
    seq_src = contents[1].rstrip('\n')

    result = []
    for seq_dest_file in seq_dest_files:
        contents = open(seq_dest_path + '/' + seq_dest_file, "r").readlines()
        dest_pdbcode = contents[0].rstrip('\n')
        seq_dest = contents[1].rstrip('\n')
        result.append([seq_src == seq_dest, src_pdbcode == dest_pdbcode, dest_pdbcode, uniprot_id])
    return result


class uniprot2pdb:

    def __init__(self,
                 pdb_fasta_dir,
                 complex_fasta_dir,
                 uniprot_official_mapping_file,
                 uniprot2pdb_mapping_file):
        self.pdb_fasta_dir = pdb_fasta_dir
        self.complex_fasta_dir = complex_fasta_dir
        self.uniprot_official_mapping_file = uniprot_official_mapping_file
        self.uniprot2pdb_mapping_file = uniprot2pdb_mapping_file
        print(self.pdb_fasta_dir)

    def update(self):

        t1 = time.time()


        print(f"1. Filtering {self.uniprot_official_mapping_file} ")

        # with open(f"{self.uniprot_official_mapping_file}_filtered", "w") as fw:
        #
        #     for line in open(self.uniprot_official_mapping_file):
        #
        #         tmp = line.split('\t')
        #
        #         if len(tmp[5].lstrip()) != 0:
        #             uniprot_id = tmp[7]
        #
        #             fw.write(f"{uniprot_id}\t{tmp[5]}\n")

        if os.path.exists(self.uniprot2pdb_mapping_file):
            os.system("rm " + self.uniprot2pdb_mapping_file)

        match_count = 0

        total_count = 0

        process_list = []

        print(f"2. Generating {self.uniprot2pdb_mapping_file}")

        processing_list_file = "process.list"
        fw = open(processing_list_file, 'w')

        for content in open(f"{self.uniprot_official_mapping_file}_filtered"):

            content = content.rstrip('\n')

            uniprot_id, pdbcodes = content.split('\t')

            pdbcodes = pdbcodes.split(';')

            for pdbcode in pdbcodes:

                pdbcode = pdbcode.lstrip().rstrip('\n')

                if pdbcode.find(':') < 0:
                    continue

                pdbcode, chainid = pdbcode.split(':')

                chain_seq_file = f"{self.pdb_fasta_dir}/{pdbcode}{chainid}.fasta"
                print(chain_seq_file)
                if not os.path.exists(chain_seq_file):
                    continue

                pdbchains = [pc for pc in os.listdir(self.complex_fasta_dir) if pc.find(pdbcode) == 0]

                process_list.append([chain_seq_file, self.complex_fasta_dir, pdbchains, uniprot_id])

                print(f"{uniprot_id}\t{pdbcode}")
                fw.write(f"{chain_seq_file}\t{self.complex_fasta_dir}\t{'_'.join(pdbchains)}\t{uniprot_id}\n")

        pool = Pool(processes=50)

        results = []

        results = pool.map(compare_seq, process_list)

        pool.close()

        pool.join()

        with open(self.uniprot2pdb_mapping_file, "w") as fw:

            for result in results:

                for subres in result:

                    seq_same, chain_same, pdbcode, uniprotid = subres

                    if seq_same:
                        fw.write(f"{uniprotid}\t{pdbcode}\n")
                        total_count += 1
                        match_count += chain_same

        t2 = time.time()

        print(f"Total pdb count: {total_count}\n Total match count: {match_count}")
        print("Uniprot ids to pdb codes mapping file generation is done. Total time:" + (t2 - t1).__str__())
