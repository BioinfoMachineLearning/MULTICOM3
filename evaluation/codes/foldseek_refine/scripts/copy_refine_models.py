import argparse, os
import sys
import logging

valid_server_dir = '/bmlfast/bml_casp15/TS_run/valid_new_mgy_ref/'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.description = "Download pdb structure from pdb bank for monomer or dimer by list"
    parser.add_argument("-o", "--outdir", help="output folder for the pdb files", type=str, required=True)
    
    args = parser.parse_args()

    print("server")
    for targetname in os.listdir(valid_server_dir):
        refinement_dir = f"{valid_server_dir}/{targetname}/N10_multimer_structure_refinement_final"
        if not os.path.exists(refinement_dir):
            continue
        print(f"cp -r {refinement_dir} {args.outdir}/{targetname}")
        os.system(f"cp -r {refinement_dir} {args.outdir}/{targetname}")
        os.system(f"cp {valid_server_dir}/{targetname}/N9_multimer_structure_evaluation/pairwise_af_avg.ranking {args.outdir}/{targetname}")
