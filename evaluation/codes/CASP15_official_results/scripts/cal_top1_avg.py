import argparse, os
import sys
import logging


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-i", "--indir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()

    group_infos = {}
    targets = []
    for result_file in os.listdir(args.indir):
        if result_file[0:2] == "TS":
            continue
        targetname = result_file[0:result_file.find('.results')]
        targets += [targetname]
        print(result_file)
        for line in open(args.indir + '/' + result_file):
            line = line.rstrip('\n')
            contents = line.split()
            groupname = contents[0][2:]
            tmscores = []
            for i in range(1, len(contents)):
                tmscores += [contents[i]]
            while len(tmscores) < 5:
                tmscores += ['0.00']

            if groupname not in group_infos:
                group_infos[groupname] = {}

            group_infos[groupname][targetname] = tmscores

    groups = sorted(group_infos)
    for targetname in targets:
        for groupname in groups:
            if targetname not in group_infos[groupname]:
                group_infos[groupname][targetname] = ['0.00'] * 5

    print(','.join(groups))
    print(','.join(targets))
    for targetname in targets:
        for i in range(5):
            scores = []
            for groupname in groups:
                scores += [group_infos[groupname][targetname][i]]
            print(','.join(scores))



