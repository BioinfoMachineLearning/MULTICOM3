#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import glob
import subprocess, argparse
import shutil

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
    parser.add_argument('--install_dir', type=str, required=True)
    parser.add_argument('--multicom3_db_dir', type=str, required=True)
    parser.add_argument('--afdb_dir', type=str, required=True)
    parser.add_argument('--outfile', type=str, required=True)
    args = parser.parse_args()

    # configure db_option file
    db_option_file_template = args.install_dir + '/docker/.db_option.default'
    newlines = []
    keywords_dict = {'INSTALLDIR_DATABASES': args.multicom3_db_dir.rstrip('/'),
                    'AFDB_DIR': args.afdb_dir.rstrip('/')}

    for line in open(db_option_file_template):
        newline = line
        for keyword in keywords_dict:
            newline = newline.replace(keyword, keywords_dict[keyword])
        newlines += [newline]
    
    with open(args.outfile, 'w') as fw:
        fw.writelines(''.join(newlines))

    print("\nConfiguration....Done")
