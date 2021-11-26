#!/bin/bash -x

if [ $# != 7 ]
then
        echo "USAGE: chainA name, chainA fasta file, chainB name, chainB fasta file, fused msa, paired msa, output folder"
        exit 1
fi


BASEDIR='/data'

WORKDIR=$7 #output folder
mkdir $WORKDIR

CHAINID1=$1
FASTAFILE1=$2 #Path to fasta file of chain1
CHAINID2=$3
FASTAFILE2=$4 #Path to fasta file of chain2

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

## Predicting
#**Chain Break and Fasta**
#E.g. seq1='AAA', seq2='BBB', catseq=AAABBB (the sequence that should be in the fasta file) and CB=3
FASTAFILE=$WORKDIR/"${CHAINID1}-${CHAINID2}.fasta"  #Path to file with concatenated fasta sequences.
echo "python3 $BASEDIR/bml_casp15/BML_CASP15/tools/FoldDock/data/merge_fasta_custom.py --fasta1 $FASTAFILE1 --fasta2 $FASTAFILE2 --mergefasta $FASTAFILE"
python3 $BASEDIR/bml_casp15/BML_CASP15/tools/FoldDock/data/merge_fasta_custom.py --fasta1 $FASTAFILE1 --fasta2 $FASTAFILE2 --mergefasta $FASTAFILE

L1=`${BASEDIR}/bml_casp15/BML_CASP15/tools/FoldDock/src/seqlen.sh $FASTAFILE1` # Length of chain1 (where to introduce chain break)
L2=`${BASEDIR}/bml_casp15/BML_CASP15/tools/FoldDock/src/seqlen.sh $FASTAFILE2` # Length of chain2
CB=$L1

MSAS="$PAIREDMSA,$FUSEDMSA" #Comma separated list of msa paths
echo "Chain break is: "$CB
echo $MSAS
#**AF2 CONFIGURATION**
AFHOME=${BASEDIR}/bml_casp15/BML_CASP15/tools/FoldDock/src/alphafold/ # Path of alphafold directory in FoldDock
SINGULARITY=${BASEDIR}/bml_casp15/BML_CASP15/tools/FoldDock/AF_data/AF_environment.sif
PARAM=${BASEDIR}/bml_casp15/BML_CASP15/tools/FoldDock/AF_data/  # Path to AF2 params
OUTFOLDER=$WORKDIR/   #Path where AF2 generates its output folder structure
PRESET="full_dbs" # Configuration1 - no ensembling (full_dbs) and (reduced_dbs) or 8 model ensemblings (casp14).
MODEL_NAME="model_1" # Configuration2
MAX_RECYCLES=3 # recycles number (default=3)

echo "Run Alphafold2"

UNICL30=$BASEDIR'/bml_casp15/BML_CASP15/databases/uniclust30_2018_08/uniclust30_2018_08'
UNIREF=$BASEDIR'/bml_casp15/BML_CASP15/databases/uniref90/uniref90_2020_03.fasta'
MGNIFY=$BASEDIR'/multicom4s_tool/alphafold_databases/data/mgnify/mgy_clusters.fa'
PDB70=$BASEDIR'/multicom4s_tool/alphafold_databases/data/pdb70/pdb70'
MMCIF=$BASEDIR'/multicom4s_tool/alphafold_databases/data/pdb_mmcif/mmcif_files/'
OBS=$BASEDIR'/multicom4s_tool/alphafold_databases/data/pdb_mmcif/obsolete.dat'
BFD=$BASEDIR'/multicom4s_tool/alphafold_databases/data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt'

 singularity exec --nv --bind $BASEDIR:$BASEDIR $SINGULARITY \
         python3 $AFHOME/run_alphafold.py \
                 --fasta_paths=$FASTAFILE \
                 --msas=$MSAS \
                 --chain_break_list=$CB \
                 --output_dir=$OUTFOLDER \
                 --model_names=$MODEL_NAME \
                 --data_dir=$PARAM \
                 --fold_only \
                 --preset=$PRESET \
                 --max_recycles=$MAX_RECYCLES

### split modelled chains in 2, renumber and rename them A and B
python3 ${BASEDIR}/bml_casp15/BML_CASP15/tools/FoldDock/src/split_chains.py --structure \
${OUTFOLDER}${CHAINID1}-${CHAINID2}/unrelaxed_model_1.pdb --outname ${OUTFOLDER}${CHAINID1}-${CHAINID2}_unranked