import argparse, os
import sys
import logging

MMalign_program = '/home/bml_casp15/BML_CASP15/tools/MMalign'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-t", "--nativedir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-i", "--tarballdir", help="pdb name in lower case", type=str, required=True)
    parser.add_argument("-o", "--outdir", help="output folder for the pdb files", type=str, required=True)
    
    args = parser.parse_args()


    for nativepdb in os.listdir(args.nativedir):

        targetname = nativepdb.rstrip('_filtered.pdb')

        print(f"Processing {targetname}:")

        for predpdb in os.listdir(args.tarballdir + '/' + targetname):

            nativepdb_full = f"{args.nativedir}/{nativepdb}"
            predpdb_full = args.tarballdir + '/' + targetname + '/' + predpdb
            outdir = args.outdir + '/' + targetname
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            cmd = f"{MMalign_program} {args.nativedir}/{nativepdb} {predpdb_full} > {outdir}/{predpdb}_out"
            os.system(cmd)
