import os
import sys

tar_list='/storage/htc/bdm/zhiye/CASP15/alphafold/example/test.lst'
fasta_dir = '/storage/htc/bdm/zhiye/CASP15/alphafold/example/'
a3m_dir = '/storage/htc/bdm/zhiye/CASP15/alphafold/example/T1064/msas/'
output_dir = '/storage/htc/bdm/zhiye/CASP15/Test/recycle_test/'
sbatch_dir = '/storage/htc/bdm/zhiye/CASP15/Test/recycle_test/sbatch/'

tar_list='/storage/htc/bdm/zhiye/CASP15/alphafold/example/test.lst'
fasta_dir = '/storage/htc/bdm/zhiye/CASP15/alphafold/example/'
a3m_dir = '/storage/htc/bdm/zhiye/CASP15/alphafold/example/7V5MF_7V5ME/msas/'
output_dir = '/storage/htc/bdm/zhiye/CASP15/Test/recycle_test2/'
sbatch_dir = '/storage/htc/bdm/zhiye/CASP15/Test/recycle_test2/sbatch/'

tar_list='/storage/htc/bdm/zhiye/DNCON4_db_tools/features/CASP14/test_fm.lst'
fasta_dir = '/storage/htc/bdm/zhiye/DNCON4_db_tools/features/CASP14/fasta/'
output_dir = '/storage/htc/bdm/zhiye/CASP15/Test/new_pdb70_test/'
sbatch_dir = '/storage/htc/bdm/zhiye/CASP15/Test/new_pdb70_test/sbatch/'

tar_list='/storage/htc/bdm/zhiye/CASP15/Test/new_pdb70_test/casp_capri.lst'
fasta_dir = '/storage/htc/bdm/zhiye/CASP15/Test/new_pdb70_test/casp_capri/fasta/'
output_dir = '/storage/htc/bdm/zhiye/CASP15/Test/new_pdb70_test/casp_capri_error2/'
sbatch_dir = '/storage/htc/bdm/zhiye/CASP15/Test/new_pdb70_test/casp_capri_error2/sbatch/'

if_gen_sbatch = True
def chkdirs(fn):
    '''create folder if not exists'''
    dn = os.path.dirname(fn)
    if not os.path.exists(dn): os.makedirs(dn)

f = open(tar_list, 'r')
chkdirs(output_dir)
chkdirs(sbatch_dir)

repeat_time = 1
use_gpu=True
use_precomputed_msas = False
model_preset = 'multimer'
for line in f.readlines():
    line_list = line.strip('\n').split(' ')
    fasta_name = line_list[0]
    print("process %s"%fasta_name)

    fasta_file = f'{fasta_dir}/{fasta_name}.fasta'
    for max_recycle in range(8, 9):
        for run_count in range(repeat_time):
            # sbatch_file = f'{sbatch_dir}/{fasta_name}_recycle{max_recycle}_{run_count}.sh'
            sbatch_file = f'{sbatch_dir}/{fasta_name}.sh'
            print(sbatch_file)
            # sub_outdir = f'{output_dir}/recycle{max_recycle}_{run_count}/{fasta_name}/'
            sub_outdir = f'{output_dir}/{fasta_name}/'
            chkdirs(sub_outdir)
            if use_precomputed_msas == True:
                af_msadir = f'{sub_outdir}/{fasta_name}/msas/'
                chkdirs(af_msadir)
                os.system(f'cp -r {a3m_dir}* {af_msadir}')
            
            if use_gpu == False:
                with open(sbatch_file, 'w') as myfile:
                    myfile.write('#!/bin/bash -l\n')
                    myfile.write('#SBATCH -J  %s-%s\n'%(fasta_name, max_recycle))
                    myfile.write('#SBATCH -o %s-%%j.out\n'%fasta_name)
                    myfile.write('#SBATCH -p Lewis,hpc4,hpc5\n')
                    myfile.write('#SBATCH -N 1\n')
                    myfile.write('#SBATCH -n 8\n')
                    myfile.write('#SBATCH -t 2-00:00\n')
                    myfile.write('#SBATCH --mem 80G\n')   
                    myfile.write('module load cuda/11.1.0\n') 
                    myfile.write('module load cudnn/8.0.4.30-11.0-linux-x64\n') 
                    myfile.write('module load miniconda3\n')
                    myfile.write('__conda_setup=\"$(\'conda\' \'shell.bash\' \'hook\' 2> /dev/null)\"\n')
                    myfile.write('eval \"$__conda_setup\"\n')
                    myfile.write('unset __conda_setup\n')
                    myfile.write('conda activate /storage/htc/bdm/zhiye/CASP15/alphafold/env/\n')

                    myfile.write('fasta_name=%s\n'%fasta_name)
                    myfile.write('fasta_file=%s\n'%fasta_file)
                    myfile.write('use_precomputed_msas=%s\n'%use_precomputed_msas)
                    myfile.write('model_preset=%s\n'%model_preset)
                    myfile.write('max_recycle=%s\n'%max_recycle)
                    myfile.write('outdir=%s\n'%sub_outdir)
                    myfile.write('scrip_dir=/storage/htc/bdm/zhiye/CASP15/alphafold/\n')
                    myfile.write('cd $scrip_dir\n')
              
                    # myfile.write('python $scrip_dir/run_alphafold_custom.py -f $fasta_file -a $aln_file  -o $outdir\n')
                    myfile.write(f'sh $scrip_dir/run_alphafold.sh -d /storage/htc/bdm/tools/alphafold/database/ -o $outdir -f $fasta_file -m $model_preset -t 2020-05-14 -p $use_precomputed_msas -c $max_recycle -g false') # -u use_precomputed_msas
            else:
                with open(sbatch_file, 'w') as myfile:
                    myfile.write('#!/bin/bash -l\n')
                    myfile.write('#SBATCH -J  %s-%s\n'%(fasta_name, max_recycle))
                    myfile.write('#SBATCH -o %s-%%j.out\n'%fasta_name)
                    myfile.write('#SBATCH --partition gpu4\n')
                    myfile.write('#SBATCH --cpus-per-task=4\n')
                    myfile.write('#SBATCH --mem-per-cpu=60G\n')
                    myfile.write('#SBATCH --time 2-00:00\n')
                    myfile.write('#SBATCH --account=engineering-gpu\n')   
                    myfile.write('#SBATCH --gres gpu:1\n')   
                    myfile.write('#SBATCH --qos gpu4\n')   
                    myfile.write('module load cuda/11.1.0\n') 
                    myfile.write('module load cudnn/8.0.4.30-11.0-linux-x64\n') 
                    myfile.write('module load miniconda3\n')
                    myfile.write('__conda_setup=\"$(\'conda\' \'shell.bash\' \'hook\' 2> /dev/null)\"\n')
                    myfile.write('eval \"$__conda_setup\"\n')
                    myfile.write('unset __conda_setup\n')
                    myfile.write('conda activate /storage/htc/bdm/zhiye/CASP15/alphafold/env/\n')

                    myfile.write('fasta_name=%s\n'%fasta_name)
                    myfile.write('fasta_file=%s\n'%fasta_file)
                    myfile.write('use_precomputed_msas=%s\n'%use_precomputed_msas)
                    myfile.write('model_preset=%s\n'%model_preset)
                    myfile.write('max_recycle=%s\n'%max_recycle)
                    myfile.write('outdir=%s\n'%sub_outdir)
                    myfile.write('scrip_dir=/storage/htc/bdm/zhiye/CASP15/alphafold/\n')
                    myfile.write('cd $scrip_dir\n')
              
                    # myfile.write('python $scrip_dir/run_alphafold_custom.py -f $fasta_file -a $aln_file  -o $outdir\n')
                    myfile.write(f'sh $scrip_dir/run_alphafold.sh -d /storage/htc/bdm/tools/alphafold/database/ -o $outdir -f $fasta_file -m $model_preset -t 2020-05-14 -p $use_precomputed_msas -c $max_recycle -g true') # -u use_precomputed_msas
                    # myfile.write(f'sh $scrip_dir/run_alphafold.sh -d /storage/htc/bdm/tools/alphafold/database/ -o $outdir -f $fasta_file -m $model_preset -t 2020-05-14 -p $use_precomputed_msas -c $max_recycle -g false') # -u use_precomputed_msas

