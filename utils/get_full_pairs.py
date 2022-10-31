import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--innpy', type=str, required=True)
    args = parser.parse_args()

    with open(args.innpy, 'rb') as f:
        paired_rows = np.load(f)

    full_pairs_count = 0
    for row_index in range(paired_rows.shape[0]):
        print(paired_rows[row_index, :])
        full_pair = True
        for col_index in range(paired_rows.shape[1]):
            if paired_rows[row_index, col_index] == -1:
                full_pair = False
                break
        if full_pair:
            full_pairs_count += 1

    print(full_pairs_count)
