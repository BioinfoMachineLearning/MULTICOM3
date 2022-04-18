source /home/bml_casp15/anaconda3/bin/activate

conda activate enqa

cd /home/bml_casp15/tools/EnQA/

python predict_batch.py --input $1 --output $2 --method $3 --alphafold_prediction $4
