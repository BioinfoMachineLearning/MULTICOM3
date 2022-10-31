import os
import sys, argparse

option_file='/storage/htc/bdm/jianliu/BML_CASP15/bin/db_option_newest_mgy_test_lewis'
fasta_dir='/storage/htc/bdm/jianliu/BML_CASP15/test_datasets/PDB/work_200/fasta'
output_dir='/storage/htc/bdm/jianliu/BML_CASP15/test_datasets/PDB/work_200/pdb'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--targetlist', type=str, required=True)
    parser.add_argument('--sbatchdir', type=str, required=True)
    args = parser.parse_args()

    for target in open(args.targetlist):
        target = target.rstrip('\n')
        fasta_file = fasta_dir + '/' + target + '.fasta'
        outdir = output_dir + '/' + target
        with open(args.sbatchdir + '/' + target + '.sh', 'w') as fw:
            fw.write('#!/bin/sh\n')
            fw.write('#SBATCH --partition gpu4\n')
            fw.write('#SBATCH --cpus-per-task=1\n')
            fw.write('#SBATCH --mem-per-cpu=100G\n')
            fw.write('#SBATCH --time 2-00:00 \n')
            fw.write('#SBATCH --account=engineering-gpu\n')
            fw.write('#SBATCH --gres gpu:1\n')
            fw.write(f'#SBATCH --job-name={target}\n')
            fw.write(f'#SBATCH --output={target}-%j.out\n')

            fw.write('__conda_setup=\"$(\'conda\' \'shell.bash\' \'hook\' 2> /dev/null)\"\n')
            fw.write('eval \"$__conda_setup\"\n')
            fw.write('unset __conda_setup\n')
            fw.write('conda activate /storage/htc/bdm/jianliu/anaconda3/envs/bml_casp15/\n')

            fw.write('module load cuda/11.1.0\n') 
            fw.write('module load cudnn/8.0.4.30-11.0-linux-x64\n') 

            fw.write(f'option_file={option_file}\n')
            fw.write(f'fasta_file={fasta_file}\n')
            fw.write(f'output_dir={outdir}\n')

            fw.write("export TF_FORCE_UNIFIED_MEMORY='1'\n")
            fw.write("export XLA_PYTHON_CLIENT_MEM_FRACTION='2.0'\n")
            fw.write("export PYTHONPATH=/storage/htc/bdm/jianliu/BML_CASP15\n")

            fw.write(f"python /storage/htc/bdm/jianliu/BML_CASP15/pipelines/test/multimer_structure_notemplates.py --option_file $option_file --fasta_file $fasta_file --output_dir $output_dir\n")
