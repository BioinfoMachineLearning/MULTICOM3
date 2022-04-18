source /home/bml_casp15/anaconda3/bin/activate base

conda activate enqa

python /home/bml_casp15/tools/EnQA/predict_batch.py --input $1 --output $2 --model $3 --alphafold_prediction $4
