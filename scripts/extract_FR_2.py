#!/usr/bin/python3
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import os
import sys
import argparse
from Scheme import SCHEMES,DEFINITIONS,NUMBERING


def main():
    parser = argparse.ArgumentParser(description='extract FR')
    parser.add_argument('--i', help=f"input file.")
    parser.add_argument('--s', help=f"Need to specify input numbering scheme. {NUMBERING}")
    parser.add_argument('--d', help=f"Need to specify input CDR definition.{DEFINITIONS}")
    parser.add_argument('--r', help=f"regions to extract.e.g FR1+FR2+FR3+FR4")
    args = parser.parse_args()
    print(args)

    fin = args.i

    scheme = args.s

    if not scheme  in NUMBERING:
        print(f"scheme must be from {NUMBERING}")
        exit(-1)
    definition =args.d
    if not definition in DEFINITIONS:
        print(f"CDR definition must be from {DEFINITIONS}")
        exit(-1)
    regions = ["FR1","CDR1","FR2","CDR2","FR3","CDR3","FR4"]
    regions_requested=dict(zip(regions,[0,0,0,0,0,0,0]))
    for region in args.r.split("+"):
        if not region in regions:
            print(f"region:{region} not recognized")
            exit(-1)
        regions_requested[region] = 1

    definition_method = SCHEMES[definition][scheme]
    dirname = os.path.dirname(fin)
    basename = os.path.basename(fin)
    
    fout = os.path.join(dirname,f'{definition}_{args.r}_'+basename.split(".")[0]+".fasta")

    data = pd.read_csv(fin)

    fout_handle = open(fout,"w+")
    seqs_pool = set()
    invalid_seq_num = 0
    redundant_seq_num = 0

    for i,row in tqdm(data.iterrows()):
        valid_seq = True
        if row['chain_type'] != "H":
            continue
        _id = row['id']
        _aas = row['aas'].split("#")
        region_seqs = [""]*7
        for ind,r in enumerate(regions):
            region_seq=""
            if regions_requested[r]==0:
                continue
            region_seq ="".join(_aas[definition_method[ind]:definition_method[ind+1]])   
            if len(region_seq) <4 or region_seq.count("X")>4:
                invalid_seq_num+=1
                valid_seq=False
                break
            region_seqs[ind] = region_seq
        if not valid_seq:
            continue
        concated_seq = "".join(region_seqs)
        if concated_seq in seqs_pool:
            redundant_seq_num+=1
            continue
        else:
            seqs_pool.add(concated_seq)
        record = SeqRecord(id=str(_id),description="",seq = Seq(concated_seq))
        SeqIO.write(record,fout_handle,"fasta")

    print(f"in total : {data.shape[0]}")
    print(f"invalid seqs : {invalid_seq_num}")
        









if __name__ == "__main__":
    main()

