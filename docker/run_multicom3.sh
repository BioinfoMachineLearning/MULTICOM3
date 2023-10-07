#!/bin/bash

ldconfig

usage() {
        echo ""
        echo "Please make sure all required parameters are given"
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-m <mode>                          monomer, heteromer or homomer"
        #echo "-op <option_file>   db_option file"
        echo "-d <multicom3_db_dir>  multicom3 database dir"
        echo "-a <af_db_dir>                alphafold database dir"
        echo "-f <fasta_path>                       Path to fasta file" 
        echo "-o <output_dir>                       output directory"
        echo "-i <run_img>                        Option, whether to use IMG alignment to generate models"
        echo ""
        exit 1
}

while getopts "m:d:a:f:o:i:" arg; do
        case $arg in
            m) mode=$OPTARG;;
            d) multicom3_db_dir=$OPTARG;;
            a) af_db_dir=$OPTARG;;
            f) fasta_path=$OPTARG;;
            o) output_dir=$OPTARG;;
            i) run_img=$OPTARG;;
        esac
done

mode=$(echo "$mode" | xargs)
multicom3_db_dir=$(echo "$multicom3_db_dir" | xargs)
af_db_dir=$(echo "$af_db_dir" | xargs)
fasta_path=$(echo "$fasta_path" | xargs)
output_dir=$(echo "$output_dir" | xargs)
run_img=$(echo "$run_img" | xargs)

# Parse input and set defaults
if [[ "$mode" == "" || "$multicom3_db_dir" == "" || "$af_db_dir" == "" || "$fasta_path" == "" || "$output_dir" == "" ]] ; then
    usage
fi

configure_script="/app/MULTICOM3/docker/configure.py"
python $configure_script --install_dir /app/MULTICOM3 --multicom3_db_dir $multicom3_db_dir --afdb_dir $af_db_dir --outfile $output_dir/db_option_docker

option_file="$output_dir/db_option_docker"
monomer_script="/app/MULTICOM3/bin/monomer.py"
heteromer_script="/app/MULTICOM3/bin/heteromer.py"
homomer_script="/app/MULTICOM3/bin/homomer.py"

export PYTHONPATH=/app/MULTICOM3

if [[ "$mode" =~ "monomer" ]] ; then
    echo "Predicting structure for monomer"
    python $monomer_script --option_file=$option_file --fasta_path=$fasta_path --output_dir=$output_dir --run_img=$run_img
fi

if [[ "$mode" =~ "heteromer" ]] ; then
    echo "Predicting structure for heteromer"
    python $heteromer_script --option_file=$option_file --fasta_path=$fasta_path --output_dir=$output_dir --run_img=$run_img
fi

if [[ "$mode" =~ "homomer" ]] ; then
    echo "Predicting structure for homomer"
    python $homomer_script --option_file=$option_file --fasta_path=$fasta_path --output_dir=$output_dir --run_img=$run_img
fi


