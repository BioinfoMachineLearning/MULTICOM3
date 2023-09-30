#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import glob
import subprocess, argparse

def makedir_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    directory = os.path.abspath(directory)
    return directory

def rm_if_exists(directory):
    if os.path.exists(directory):
        directory = os.path.abspath(directory)
        if os.path.isdir(directory):
            os.system("rm -r "+directory)
        else:
            os.system("rm "+directory)

def direct_download(tool, address, tools_dir):  ####Tools don't need to be configured after downloading and configuring
    os.chdir(tools_dir)
    tool_dir = os.path.join(tools_dir, tool)
    if not os.path.exists(tool_dir):
        rm_if_exists(tool_dir)
        os.system("wget "+address)
        print("Decompressing "+tool_dir)
        os.system("tar -zxf "+tool+".tar.gz && rm "+tool+".tar.gz")
        os.system("chmod -R 755 "+tool_dir)
        print("Downloading "+tool_dir+"....Done")
    else:
        print(tool+" has been installed "+tool_dir+"....Skip....")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--envdir', type=str, required=True)
    parser.add_argument('--af_dir', type=str, required=True)
    parser.add_argument('--afdb_dir', type=str, required=True)
    args = parser.parse_args()

    # Set directory of multicom databases and tools
    install_dir = os.path.dirname(os.path.realpath(__file__))
    database_dir = os.path.join(install_dir, "databases")
    tools_dir = os.path.join(install_dir, "tools")
    bin_dir = os.path.join(install_dir, "bin")
    log_dir = os.path.join(install_dir, "installation/log")
    makedir_if_not_exists(database_dir)
    makedir_if_not_exists(tools_dir)
    makedir_if_not_exists(bin_dir)
    makedir_if_not_exists(log_dir)

    args.af_dir = os.path.abspath(args.af_dir)
    args.afdb_dir = os.path.abspath(args.afdb_dir)
    
    print("MULTICOM3 database path : "+ database_dir)
    print("MULTICOM3 tool path : "+ tools_dir)

    ### (1) Download basic tools
    os.chdir(tools_dir)
    tools_lst = ["DockQ", "foldseek", "mmalign", "pairwiseQA", "tmalign", "tmscore", "deepmsa", "colabfold", "hhsuite-3.2.0"]
    for tool in tools_lst:
        logfile = os.path.join(log_dir, tool + '.done')
        if os.path.exists(logfile):
            print(os.path.join(log_dir, tool) +" installed....skip")
        else:
            runfile = os.path.join(log_dir, tool + '.running')
            os.system(f"touch {runfile}")
            address = "http://sysbio.rnet.missouri.edu/multicom_cluster/multicom3_db_tools/tools/"+tool+".tar.gz"
            direct_download(tool, address, tools_dir)
            os.system(f"mv {runfile} {logfile}")
            print(os.path.join(log_dir, tool) + " installed")

    ### (2) Download databases
    os.chdir(database_dir)

    img_fasta_url = "http://sysbio.rnet.missouri.edu/multicom_cluster/multicom3_db_tools/databases/img/img_new.fasta"
    img_dir = "img"
    makedir_if_not_exists(img_dir)
    os.chdir(img_dir)
    if not os.path.exists("img_new.fasta"):
        os.system("wget " + img_fasta_url)
    print("Downloading img ....Done")

    #### Download db_lst
    db_lst = ["af_pdbs","Metaclust_2018_06","myg_uniref100_04_2020","pdb_complex","pdb_sort90","string","uniprot2pdb"]
    for db in db_lst:
        print("Download "+db)
        if os.path.exists(os.path.join(database_dir, db)):
            continue
        direct_download(db,"http://sysbio.rnet.missouri.edu/multicom_cluster/multicom3_db_tools/databases/"+db+".tar.gz",database_dir)

    #### Download uniclust30_2018_08_hhsuite
    print("Download uniclust30_2018_08_hhsuite\n")
    if not os.path.exists(os.path.join(database_dir, 'uniclust30_2018_08')):
        direct_download("uniclust30_2018_08_hhsuite","https://gwdu111.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz",database_dir)

    #### Download ColabFold database
    print("Download ColabFold database\n")
    colabfold_db_dir = os.path.join(database_dir, 'colabfold')
    if not os.path.exists(colabfold_db_dir):
        os.makedirs(colabfold_db_dir)
        os.system(f"sh {os.path.join(tools_dir, "/colabfold/setup_databases.sh")} {colabfold_db_dir}")

    ### (3) copy the alphafold-addon scripts
    alphafold_addon_dir = os.path.join(install_dir, 'alphafold_addon')
    if not os.path.exists(alphafold_addon_dir):
        raise Exception(f"Cannot find alphafold_addon_dir: {alphafold_addon_dir}")
    
    src_to_trg_dict = {'customized/data_custom': 'alphafold/data_custom',
                        'customized/run_alphafold_custom.py': 'run_alphafold_custom.py',
                        'customized/run_alphafold_multimer_custom.py': 'run_alphafold_multimer_custom.py',
                        'default/run_alphafold_pre.py': 'run_alphafold_pre.py',
                        }

    for srcfile in src_to_trg_dict:
        trgfile = os.path.join(args.af_dir, src_to_trg_dict[srcfile])
        if os.path.exists(trgfile):
            os.system(f'rm -rf {trgfile}')
        os.system(f"cp -r {os.path.join(alphafold_addon_dir, srcfile)} {trgfile}")

    # (4) configure db_option file
    db_option_file_template = os.path.join(install_dir, 'bin/.db_option.default')
    newlines = []
    keywords_dict = {'YOUR_ENV': os.path.abspath(args.envdir.rstrip('/')),
                    'INSTALLDIR_TOOLS': os.path.join(install_dir, 'tools'),
                    'INSTALLDIR_DATABASES': os.path.join(install_dir, 'databases'),
                    'AFDB_DIR': args.afdb_dir.rstrip('/'),
                    'AF_DIR': args.af_dir.rstrip('/')}

    for line in open(db_option_file_template):
        newline = line
        for keyword in keywords_dict:
            newline = newline.replace(keyword, keywords_dict[keyword])
        newlines += [newline]
    
    with open(os.path.join(install_dir, '/bin/db_option'), 'w') as fw:
        fw.writelines(''.join(newlines))

    print("\nConfiguration....Done")
