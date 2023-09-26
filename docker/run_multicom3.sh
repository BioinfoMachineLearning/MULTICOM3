#!/bin/bash

ldconfig

usage() {
        echo ""
        echo "Please make sure all required parameters are given"
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-mode <mode>                          monomer, heteromer or homomer"
        #echo "-op <option_file>   db_option file"
        echo "-multicom3_db_dir <multicom3_db_dir>  multicom3 database dir"
        echo "-af_db_dir <af_db_dir>                alphafold database dir"
        echo "-fasta_path <fasta_path>                       Path to fasta file" 
        echo "-output_dir <output_dir>                       output directory"
        echo "-run_img <run_img>                        Option, whether to use IMG alignment to generate models"
        echo ""
        exit 1
}

while getopts ":mode:multicom3_db_dir:af_db_dir:fasta_path:output_dir:run_img" i; do
        case "${i}" in
        mode)
                mode=$OPTARG
        ;;
        multicom3_db_dir)
                multicom3_db_dir=$OPTARG
        ;;
        af_db_dir)
                af_db_dir=$OPTARG
        ;;
        fasta_path)
                fasta_path=$OPTARG
        ;;
        output_dir)
                output_dir=$OPTARG
        ;;
        run_img)
                run_img=$OPTARG
        ;;
        esac
done

# Parse input and set defaults
if [[ "$mode" == "" || "$multicom3_db_dir" == "" || "$af_db_dir" == "" || "$fasta_path" == "" || "$output_dir" == "" ]] ; then
    usage
fi

if [[ "$run_img" == "" ]] ; then
    run_img="False"
fi

configure_script = "/app/MULTICOM3/configure.py"
python $configure_script --install_dir=/app/MULTICOM3 --multicom3_db_dir=$multicom3_db_dir --afdb_dir=$af_db_dir

option_file="/app/MULTICOM3/bin/option_file"
monomer_script="/app/MULTICOM3/bin/monomer.py"
heteromer_script="/app/MULTICOM3/bin/heteromer.py"
homomer_script="/app/MULTICOM3/bin/homomer.py"

if [ "$mode" == "monomer" ]
then  
    python $monomer_script --option_file=$option_file --fasta_path=$fasta_path --run_img=$run_img --output_dir=$output_dir
elif [ "$mode" == "heteromer" ]
then
    python $heteromer_script --option_file=$option_file --fasta_path=$fasta_path --run_img=$run_img --output_dir=$output_dir
else
    python $homomer_script --option_file=$option_file --fasta_path=$fasta_path --run_img=$run_img --output_dir=$output_dir
fi
 

