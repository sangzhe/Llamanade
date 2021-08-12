#!/usr/bin/python3
from Bio import SeqIO
from tqdm import tqdm
from multiprocessing import Pool
import numpy as np
import pandas as pd
import os
import sys
import re
import simplejson as json

infile = sys.argv[1]
excluded_cols = ['Id','domain_no','hmm_species','chain_type','e-value','score','seqstart_index','seqend_index','identity_species','v_gene','v_identity','j_gene','j_identity']
def count_gap(x):
    n = 0
    for i in x:
        if i =='-':
            n+=1
    return n
data = pd.read_csv(infile)

data.drop(excluded_cols,axis=1,inplace=True)
basename = os.path.basename(infile)
dirname = os.path.dirname(infile)
outname = basename.split(".")[0]+".json"
full_out_path = os.path.join(dirname,outname)
result = dict()
with open(full_out_path, 'w') as f:
    for col in tqdm(data.columns):
        gap_ratio =  count_gap(data[col])/data.shape[0]
        #if gap_ratio>0.9999:
        #    print(f'Ignore positioon {col},gap frequency {gap_ratio}')
        #    continue
        count_result = data[col].value_counts(normalize=True)
        values = count_result.keys().tolist()
        counts = count_result.tolist()
        value_dict = dict(zip(values, counts))
        result[col] = value_dict
    json.dump(result, f)
