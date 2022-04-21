source /home/bml_casp15/anaconda3/bin/activate base

conda activate enqa

cd /home/bml_casp15/tools/EnQA/

if [[ $# -eq 5 ]] && [[ $5 -eq 'false' ]]
then
  python predict_batch.py --input $1 --output $2 --method $3 --alphafold_prediction $4 --cpu
else
  python predict_batch.py --input $1 --output $2 --method $3 --alphafold_prediction $4
fi


