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
    parser.add_argument('--multicom3db_dir', type=str, required=True)
    parser.add_argument('--afdb_dir', type=str, required=True)
    args = parser.parse_args()

    # Set directory of multicom databases and tools
    install_dir = os.path.dirname(os.path.realpath(__file__))
    args.multicom3db_dir = os.path.abspath(args.multicom3db_dir)
    args.afdb_dir = os.path.abspath(args.afdb_dir)
    
    # configure db_option file
    db_option_file_template = os.path.join(install_dir, 'bin/.db_option.default')
    newlines = []
    keywords_dict = {'YOUR_ENV': os.path.abspath(args.envdir.rstrip('/')),
                    'INSTALLDIR_TOOLS': os.path.join(install_dir, 'tools'),
                    'INSTALLDIR_DATABASES': args.multicom3db_dir.rstrip('/'),
                    'AFDB_DIR': args.afdb_dir.rstrip('/'),
                    'AF_DIR': os.path.join(install_dir, 'tools/alphafold-v2.3.2')}

    for line in open(db_option_file_template):
        newline = line
        for keyword in keywords_dict:
            newline = newline.replace(keyword, keywords_dict[keyword])
        newlines += [newline]
    
    with open(os.path.join(install_dir, 'bin/db_option'), 'w') as fw:
        fw.writelines(''.join(newlines))

    print("\nConfiguration....Done")
