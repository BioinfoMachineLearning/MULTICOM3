#!/bin/bash
name1=$1
fasta1=$(readlink -f $2)
name2=$3
fasta2=$(readlink -f $4)
outdir_=$5
#outdir=$(readlink -f $22)


if [ ! -d $outdir_ ]; then
    mkdir -p $outdir_
    echo "$?"
fi

outdir=$(readlink -f $5)

source /home/multicom4s_tool/anaconda3/bin/activate
conda activate alphafold

FASTAFILE=$outdir/"${name1}-${name2}.fasta"  #Path to file with concatenated fasta sequences.
python3 /data/bml_casp15/BML_CASP15/tools/alphafold/merge_fasta.py --fasta1 $FASTAFILE1 --fasta2 $FASTAFILE2 --mergefasta $FASTAFILE

echo "FASTA $FASTAFILE"
echo "OUTDIR $outdir"

#exit

python3 /data/bml_casp15/BML_CASP15/tools/alphafold/docker/run_docker.py \
        --fasta_paths=$FASTAFILE \
        --is_prokaryote_list=false \
        --max_template_date=2020-05-14 \
        --model_preset=multimer \
        --data_dir /data/bml_casp15/BML_CASP15/databases/alphafold_databases \
        --output_dir $outdir
