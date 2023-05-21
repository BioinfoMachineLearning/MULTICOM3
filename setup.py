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
    if not os.path.exists(tools_dir+"/"+tool):
        rm_if_exists(tools_dir+"/"+tool)
        os.system("wget "+address)
        print("Decompressing "+tools_dir+"/"+tool)
        os.system("tar -zxf "+tool+".tar.gz && rm "+tool+".tar.gz")
        os.system("chmod -R 755 "+tools_dir+"/"+tool)
        print("Downloading "+tools_dir+"/"+tool+"....Done")
    else:
        print(tool+" has been installed "+tools_dir+"/"+tool+"....Skip....")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

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

    print("MULTICOM3 database path : "+ database_dir)
    print("MULTICOM3 tool path : "+ tools_dir)

    ### (1) Download databases
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
        direct_download(db,"http://sysbio.rnet.missouri.edu/multicom_cluster/multicom3_db_tools/databases/"+db+".tar.gz",database_dir)

    #### Download uniclust30_2018_08_hhsuite
    print("Download uniclust30_2018_08_hhsuite\n")
    direct_download("uniclust30_2018_08","https://gwdu111.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz",database_dir)
   
    

    ### (2) Download basic tools
    os.chdir(tools_dir)
    tools_lst = ["DockQ", "foldseek", "mmalign", "mmseqs", "pairwiseQA", "tmalign", "tmscore", "deepmsa"]
    for tool in tools_lst:
        if os.path.exists(log_dir+"/"+tool+".done"):
            print(log_dir+"/"+tool+" installed....skip")
        else:
            os.system("touch "+log_dir+"/"+tool+".running")
            address = "http://sysbio.rnet.missouri.edu/multicom_cluster/multicom3_db_tools/tools/"+tool+".tar.gz"
            direct_download(tool, address, tools_dir)
            os.system("mv "+log_dir+"/"+tool+".running "+log_dir+"/"+tool+".done")
            print(log_dir+"/"+tool+" installed")


    print("\nConfiguration....Done")