#!/usr/bin/python3
from Bio import SeqIO
from tqdm import tqdm
from multiprocessing import Pool
import numpy as np
import os
import sys
import re
import logging
import argparse
def annotate(seq_obj):
    seq = seq_obj.seq
    text = os.popen('ANARCI -i {} --r heavy --ncpu 4 --scheme martin --assign_germline'.format(seq)).read()
    if len(text.splitlines())<7:
        logger.warning(f"cannot align to heavy chain:{seq_obj.id}")
        return None
    #_1:Seq Id
    #_2:(text)ANARCI numbered
    #_3:(text)Domain
    #_4:(text)Most significant hit
    #_5:(text)|species|chain_type|e-value|score|seqstart_index|seqend_index|
    #_6:||$|$|$|$|$
    #_7:(text)Most sequence-identical germlines
    #_8:(text)|species|v_gene|v_identity|j_gene|j_identy
    #_9:||$|$|$|$|$
    #_10:Scheme =$
    #*_11: numbering
    # H 1   Q
    # H 2   V
    _1,_2,_3,_4,_5,_6,_7,_8,_9,_10,*_11 = text.splitlines()

    blank, _hmm_species, _chain_type, _e_value, _score, _seqstart_index, _seqend_index,b = _6.split("|")
    blank,_germline_species,_v_gene,_v_identity,_j_gene,_j_identity,_t = _9.split("|")
    # if SPECIES!="" and _germline_species!=SPECIES:
    #     logger.warning(f"wrong (germline)species:{_germline_species},expected:{SPECIES}:{seq_obj.id}")
    #     return None
    if float(_score)<70:
        logger.warning(f"score < 70:{seq_obj.id}")
        return None
    prefix = seq[:int(_seqstart_index)]
    suffix = seq[int(_seqend_index)+1:]
    sequence = [""] * 114
    for resi in _11:
        if resi == '//':
            break
        _aa = resi[-1]
#        if _aa =='-':
#            continue
        try:
            _resn = int(re.sub(r"\D","",resi[1:-1]))
            _is_insertion = resi[-3].isalpha()
            if _aa == "-":
                _aa = "X" if not _is_insertion else ""
            sequence[_resn] = sequence[_resn]+_aa
        except:
            logger.warning(f'{_11}:{seq_obj.id}')
            return None
    return (_chain_type,_germline_species,_v_gene,_v_identity,_j_gene,_j_identity,prefix,sequence,suffix,_score,_e_value)

"""
    FR1 =len("".join(sequence[:25]))
    CDR1=len("".join(sequence[25:41]))
    FR2 =len("".join(sequence[41:58]))
    CDR2=len("".join(sequence[58:78]))
    FR3 =len("".join(sequence[78:109]))
    CDR3=len("".join(sequence[109:138]))
    FR4 =len("".join(sequence[138:]))
"""
    #return (_score,np.cumsum([int(_seqstart_index)+FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4+len(seq)-int(_seqend_index)-1]))
    #return (_score,np.cumsum([int(_seqstart_index),FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4,int(_seqend_index)]))

class Sequence:
    def __init__(self,_id,_seq):
        self.id = _id
        self.seq = _seq
    def add_result(self,result):
        self.chain=result[0]
        self.germline_species = result[1]
        self.v_gene = result[2]
        self.v_identity = result[3]
        self.j_gene = result[4]
        self.j_identity = result[5]
        self.lead_seq = result[6]
        self.aas=result[7]
        self.tail_seq = result[8]
        self.score=result[9]
        self.evalue = result[10]

    def __str__(self):
        row = f'{self.id},{self.chain},{self.germline_species},{self.score},{self.evalue},{self.v_gene},{self.v_identity},{self.j_gene},{self.j_identity},{self.lead_seq},{"#".join(self.aas)},{self.tail_seq}\n'
        return row
    def __repr__(self):
        return self.__str__()

def read_fasta(fin):
    for record in SeqIO.parse(fin,'fasta'):
        sequence = Sequence(str(record.id),str(record.seq))
        yield sequence
def write_fasta(seqs,fout):
    for seq in seqs:
        if seq == None:
            continue
        fout.write(str(seq))
    fout.flush()

def process(seq):
    result = annotate(seq)
    if result == None:
        return 

    seq.add_result(result)
    return seq

logging.basicConfig(\
            stream=sys.stdout,\
            level=logging.INFO,\
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',\
            datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('CDR Annotation')

def main():
    fin = sys.argv[1]
    seqs = list(read_fasta(fin))
    logger.info("finish reading {} sequences\n".format(len(seqs)))
    dirname = os.path.dirname(fin)
    basename = os.path.basename(fin).split(".")[0]+".csv"

    fout_path = os.path.join(dirname,f"Martin_annotated_{basename}")

    fout=open(fout_path,"w+")
    fout.write("id,chain_type,germline_species,score,e_value,v_gene,v_idneity,j_gene,j_identity,lead_seq,aas,tail_seq\n")
    SEQ_NUM=10000
    p = Pool(12)
    bins = list(range(0,len(seqs),SEQ_NUM))+[len(seqs)]
    for l in range(len(bins)-1):
        _left = bins[l]
        _right = bins[l+1]
        batch = seqs[_left:_right]
        logger.info("process {} sequences".format(_right))
        annotated_seqs = p.map(process,batch)
        write_fasta(annotated_seqs,fout)
    p.close()

main()
