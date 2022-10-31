#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition gpu4
#SBATCH --cpus-per-task=1  # cores per task
#SBATCH --mem-per-cpu=100G  # memory per core (default is 1GB/core)
#SBATCH --time 2-00:00     # days-hours:minutes time 
#SBATCH --account=engineering-gpu  # investors in gpu4 will replace this (e.g. engineering-gpu)
#SBATCH --gres gpu:1
## labels and outputs
#SBATCH --job-name=CASP_GPU
#SBATCH --output=CASP_GPU-%j.out  # %j is the unique jobID

__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
conda activate /storage/htc/bdm/jianliu/anaconda3/envs/bml_casp15/
echo 'Activate alphafold environment!'
module load cuda/11.1.0
module load cudnn/8.0.4.30-11.0-linux-x64

option_file=/storage/htc/bdm/jianliu/BML_CASP15/bin/db_option_newest_mgy_test_lewis
fasta_file=/storage/htc/bdm/jianliu/BML_CASP15/test_datasets/PDB/fasta/7N1J.fasta
output_dir=/storage/htc/bdm/jianliu/BML_CASP15/test_datasets/PDB/outdir_200/7N1J/
# TensorFlow control
export TF_FORCE_UNIFIED_MEMORY='1'

# JAX control
export XLA_PYTHON_CLIENT_MEM_FRACTION='2.0'

export PYTHONPATH=/storage/htc/bdm/jianliu/BML_CASP15

# sh $alphafold_dir/run_alphafold.sh -d $data_dir -o $output_dir -f $fasta_file -t 2022-10-10 -g true -m monomer -c 8 -p true
python /storage/htc/bdm/jianliu/BML_CASP15/pipelines/test/multimer_structure.py --option_file $option_file --fasta_file $fasta_file --output_dir $output_dir
