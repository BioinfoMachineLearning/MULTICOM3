import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm


def direct_download(tools_dir, address, tool):  ####Tools don't need to be configured after downloading and configuring
    if not os.path.exists(tools_dir):
        os.system("mkdir " + tools_dir)
        os.chdir(tools_dir)
        os.system("wget "+address)

        print("Decompressing "+tool)
        os.system(f"tar -zxf {tool} && rm {tool}")
        os.system("chmod -R 755 "+tools_dir)
        print("Downloading "+tool+"....Done")
    else:
        print(tool+" has been installed "+tool+"....Skip....")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--option_file', type=is_file, required=True)

    args = parser.parse_args()

    #### Download UniRef30_2020_01
    uniref_db = params['uniref_db']
    print(f"Download {uniref_db}\n")

    db_date = uniref_db.split('_', 2)[1]

    direct_download(params['uniref_db_dir'],
                    f"http://gwdu111.gwdg.de/~compbiol/uniclust/{db_date}/{uniref_db}_hhsuite.tar.gz",
                    f"{uniref_db}_hhsuite.tar.gz")

    os.chdir(params['uniref_db_dir'])
    if os.path.exists(f"{uniref_db}_a3m_db"):
        os.system(f"unlink {uniref_db}_a3m_db")
    if os.path.exists(f"{uniref_db}_hhm_db"):
        os.system(f"unlink {uniref_db}_hhm_db")
    os.system(f"ln -s {uniref_db}_a3m.ffdata {uniref_db}_a3m_db")
    os.system(f"ln -s {uniref_db}_hhm.ffdata {uniref_db}_hhm_db")

    #### Download UniRef30_2020_01.fasta

    direct_download(params['uniref_db_dir'],
                    f"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-{db_date}/"
                    f"uniref/uniref{db_date}.tar.gz",
                    f"uniref/uniref{db_date}.tar.gz")

    os.system("rm uniref100.fasta uniref50.fasta")
    os.system(f"mv uniref90.fasta {uniref_fasta}")

    #### Download Uniclust30_2018_08
    uniclust_db = params['uniclust_db']
    print(f"Download {uniclust_db}\n")

    db_date = uniclust_db.split('_', 2)[1]

    direct_download(params['uniclust_db_dir'],
                    f"http://gwdu111.gwdg.de/~compbiol/uniclust/{db_date}/{uniclust_db}_hhsuite.tar.gz",
                    f"{uniclust_db}_hhsuite.tar.gz")

    os.chdir(params['uniclust_db_dir'])
    if os.path.exists(f"{uniclust_db}_a3m_db"):
        os.system(f"unlink {uniclust_db}_a3m_db")
    if os.path.exists(f"{uniclust_db}_hhm_db"):
        os.system(f"unlink {uniclust_db}_hhm_db")
    os.system(f"ln -s {uniclust_db}_a3m.ffdata {uniclust_db}_a3m_db")
    os.system(f"ln -s {uniclust_db}_hhm.ffdata {uniclust_db}_hhm_db")


    #### Generate mapping files
    os.system(f"python SOFTWARE_PATH/bml_casp15/database/generate_mapping_geno_dist.py "
              f"--option_file {args.option_file}")

    os.system(f"python SOFTWARE_PATH/bml_casp15/database/generate_mapping_uniprot2pdb.py "
              f"--option_file {args.option_file}")

    os.system(f"python SOFTWARE_PATH/bml_casp15/database/generate_mapping_uniprot2string.py "
              f"--option_file {args.option_file}")

