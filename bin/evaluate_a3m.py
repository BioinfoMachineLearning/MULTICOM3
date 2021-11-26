import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm

def combine_pdb(pdb1, pdb2, combine_pdb):
    with open(combine_pdb, 'w') as out:
        for line in open(pdb1):
            if not line.startswith('ATOM'): continue
            out.write(line[:21] + 'A' + line[22:])
        for line in open(pdb2):
            if not line.startswith('ATOM'): continue
            out.write(line[:21] + 'B' + line[22:])


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dimerlist', type=is_file, required=True)
    parser.add_argument('--atom_dir', type=is_dir, required=True)
    parser.add_argument('--output_dir', type=str, required=True)

    args = parser.parse_args()

    for dimer in args.output_dir:

        print(f"Processing {dimer}")
        chain1, chain2 = dimer.split('_')

        chain1_atom = f"{args.atom_dir}/{chain1}.atom"

        chain2_atom = f"{args.atom_dir}/{chain2}.atom"

        combine_atom = f"{args.output_dir}/{dimer}/{chain1}_{chain2}.atom"

        combine_pdb(chain1_atom, chain2_atom, combine_atom)

        for method in ['uniclust_oxmatch_a3m',
                                'pdb_interact_uniref_a3m',
                                'species_interact_uniref_a3m',
                                'uniprot_distance_uniref_a3m',
                                # 'string_interact_uniref_a3m',
                                # 'geno_dist_uniref_a3m',
                                'pdb_interact_uniref_sto',
                                'species_interact_uniref_sto',
                                'uniprot_distance_uniref_sto',
                                # 'string_interact_uniref_sto',
                                # 'geno_dist_uniref_sto',
                                'pdb_interact_uniprot_sto',
                                'species_interact_uniprot_sto',
                                'uniprot_distance_uniprot_sto']:
            dimer_chainA = chain1[4]
            dimer_chainB = chain2[4]
            dockQ_scores = []
            for i in range(1,5):

                model = f'{args.ouput_dir}/{dimer}/{method}/unrelaxed_model_{i}_multimer.pdb'

                score = os.popen(f"python /data/bml_casp15/BML_CASP15/tools/DockQ/DockQ.py "
                                 f"-native_chain1 {dimer_chainA} -model_chain1 A "
                                 f"-native_chain2 {dimer_chainB} -model_chain2 B")

                dockQ_scores += [score]

            mean_score = np.mean(np.array(dockQ_scores))

            print(f"{method}\t{mean_score}")