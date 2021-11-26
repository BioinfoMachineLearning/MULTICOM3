#!/bin/sh

. "/data/bml_casp15/BML_CASP15/tools/colabfold/conda/etc/profile.d/conda.sh"
conda activate /data/bml_casp15/BML_CASP15/tools/colabfold/colabfold-conda
export NVIDIA_VISIBLE_DEVICES="all"
export TF_FORCE_UNIFIED_MEMORY="1"
export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"
export COLABFOLD_PATH="/data/bml_casp15/BML_CASP15/tools/colabfold"

if [ $# != 7 ]
then
        echo "USAGE: chainA name, chainA fasta file, chainB name, chainB fasta file, fused msa, paired msa, output folder"
        exit 1
fi

CHAINID1=$1
FASTAFILE1=$2 #Path to fasta file of chain1
CHAINID2=$3
FASTAFILE2=$4 #Path to fasta file of chain2
WORKDIR=$7 #output folder

FUSEDMSA=$WORKDIR/"${CHAINID1}-${CHAINID2}_fused.a3m"
if [ ! -s $FUSEDMSA ]
then
        cp $5 $FUSEDMSA
fi

PAIREDMSA=$WORKDIR/"${CHAINID1}-${CHAINID2}_paired.a3m"
if [ ! -s $PAIREDMSA ]
then
        cp $6 $PAIREDMSA
fi

MSAS="$PAIREDMSA,$FUSEDMSA" #Comma separated list of msa paths

FASTAFILE=$WORKDIR/"${CHAINID1}-${CHAINID2}.fasta"  #Path to file with concatenated fasta sequences.
python3 /data/bml_casp15/BML_CASP15/tools/colabfold/merge_fasta.py --fasta1 $FASTAFILE1 --fasta2 $FASTAFILE2 --mergefasta $FASTAFILE

UNIREF='/data/bml_casp15/BML_CASP15/databases/uniref90/uniref90_2021_06.fasta'
MGNIFY='/data/bml_casp15/BML_CASP15/databases/mgnify/mgy_clusters.fa'
SMALLBFD='/data/bml_casp15/BML_CASP15/databases/smallbfd/bfd-first_non_consensus_sequences.fasta'

cd /data/bml_casp15/BML_CASP15/tools/colabfold
python3.7 /data/bml_casp15/BML_CASP15/tools/colabfold/runner_af2advanced.py \
          --pair_mode unpaired+paired \
          --add_custom_msa 1 \
          --custom_msas $MSAS \
          --only_use_custom_msas 1 \
          --input $FASTAFILE \
          --output_dir $WORKDIR \
          --homooligomer 1:1 \
          --use_ptm --use_turbo \
          --max_recycle 3 \
          --num_relax Top5 \
          --msa_method jackhmmer \
          --uniref90_fasta $UNIREF \
          --smallbfd $SMALLBFD \
          --mgnify $MGNIFY