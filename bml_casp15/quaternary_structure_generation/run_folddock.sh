#!/bin/bash -x

if [ $# != 5 ]
then
        echo "USAGE: chainA name, chainA fasta file, chainB name, chainB fasta file, output folder"
        exit 1
fi


BASEDIR='/data'

HHBLITS=$BASEDIR/bml_casp15/BML_CASP15/tools/hhsuite-3.1.0-SSE2-Linux/bin/hhblits #Path to hhblits version 3.1.0 \
UNICLUST30=$BASEDIR/bml_casp15/BML_CASP15/databases/uniclust30_2018_08/uniclust30_2018_08

#For each of the two chains, run HHblits against Uniclust30
echo "Generating a3m for chainA"
WORKDIR=$5 #output folder
mkdir $WORKDIR
CHAINID1=$1
FASTAFILE1=$2 #Path to fasta file of chain1
OUTNAME1=$WORKDIR/$CHAINID1".a3m"
if [ ! -s $OUTNAME1 ]
then
        $HHBLITS -i $FASTAFILE1 -d $UNICLUST30 -E 0.001 -all -oa3m $OUTNAME1
fi

echo "Generating a3m for chainB"
CHAINID2=$3
FASTAFILE2=$4 #Path to fasta file of chain2
OUTNAME2=$WORKDIR/$CHAINID2".a3m"
if [ ! -s $OUTNAME2 ]
then
        $HHBLITS -i $FASTAFILE2 -d $UNICLUST30 -E 0.001 -all -oa3m $OUTNAME2
fi

# #Create two input MSAs (paired and fused) from the HHblits results for each chain
echo "Paired a3m for chainA and chainB"
A3M1=$OUTNAME1 #Path to a3m from chain 1 (from step 1)
A3M2=$OUTNAME2 #Path to a3m from chain 2 (from step 1)
MGF=0.9        #The max gap fraction allowed in the sequences
PAIREDMSA=$WORKDIR/"${CHAINID1}-${CHAINID2}_paired.a3m"
echo $PAIREDMSA

if [ ! -s $PAIREDMSA ]
then
        python3 ${BASEDIR}/bml_casp15/BML_CASP15/tools/FoldDock/src/oxmatch.py --a3m1 $A3M1 --a3m2 $A3M2 --max_gap_fraction $MGF --outname $PAIREDMSA
fi

echo "Fused a3m for chainA and chainB"
A3M1=$OUTNAME1 #Path to a3m from chain 1 (from step 1)
A3M2=$OUTNAME2 #Path to a3m from chain 2 (from step 1)
MGF=0.9        #The max gap fraction allowed in the sequences
FUSEDMSA=$WORKDIR/"${CHAINID1}-${CHAINID2}_fused.a3m"
if [ ! -s $FUSEDMSA ]
then
        python3 ${BASEDIR}/bml_casp15/BML_CASP15/tools/FoldDock/src/fuse_msas.py --a3m1 $A3M1 --a3m2 $A3M2 --max_gap_fraction $MGF --outname $FUSEDMSA
fi


## Predicting
#**Chain Break and Fasta**
#E.g. seq1='AAA', seq2='BBB', catseq=AAABBB (the sequence that should be in the fasta file) and CB=3
FASTAFILE=$WORKDIR/"${CHAINID1}-${CHAINID2}.fasta"  #Path to file with concatenated fasta sequences.
python3 $BASEDIR/bml_casp15/BML_CASP15/tools/FoldDock/data/merge_fasta_custom.py --fasta1 $FASTAFILE1 --fasta2 $FASTAFILE2 --mergefasta $FASTAFILE
CB=$?
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

echo "Run Alpahfold2"

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