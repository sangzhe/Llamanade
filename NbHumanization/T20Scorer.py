from NbHumanization.Sequence import NumberedSequence
import tempfile
import os
from NbHumanization import params
import shutil
from NbHumanization.Annotation import annotate
from NbHumanization.IOUtils import get_SeqRecord_from_fasta,create_output_folder
import subprocess
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import sys
logging.basicConfig(\
            filename="NbHumanization.log",\
            #stream=sys.stdout,\
            level=logging.INFO,\
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',\
            datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger('T20Scorer')
logger.setLevel(logging.INFO)


def run_blastp_T20(seq:str,dest_dir:str,blastdb:str):
    input_tmp = os.path.join(dest_dir,"tmp_input_T20.fasta")
    output_blastp = os.path.join(dest_dir,"Blastp_T20.result")
    with open(input_tmp,"w+") as f:
        f.write(f">seq\n{seq}")
    #input_tmp.close()
    args = [params.BLAST_EXEC,"-db",blastdb,"-query",input_tmp,"-max_target_seqs","20","-outfmt","7 pident bitscore ppos score","-num_threads", "10","-out",output_blastp]
    subprocess.check_call(args)
    os.remove(input_tmp)
    return output_blastp

def parse_blastp(blastp_result:str):
    parsed_T20_result = os.path.join(os.path.dirname(blastp_result),"blastp_processed_result.csv")
    args = [params.LLAMANADE+"/NbHumanization/blastp_parser.py","--input",blastp_result,"--out",parsed_T20_result]
    subprocess.check_call(args)
    return parsed_T20_result

def score_fasta(file:str,by="FR"):
    try:
        seqrecord = get_SeqRecord_from_fasta(file)
        dest_dir = create_output_folder()
        numbered_seq = annotate(seqrecord,params.ANNOTATION_SCHEME,params.ANNOTATION_DEFINITION,dest_dir)
        s = score(numbered_seq,dest_dir=dest_dir,by=by)
        print(s)
    finally:
        #pass
        shutil.rmtree(dest_dir)
def score_plain_seq(seq:str,dest_dir:str,by="FR"):
    seqrecord = SeqRecord(id="tmp",seq=Seq(seq))
    numbered_seq = annotate(seqrecord,params.ANNOTATION_SCHEME,params.ANNOTATION_DEFINITION,dest_dir)
    return score(numbered_seq,dest_dir,by)

def score(numbered_seq:NumberedSequence,dest_dir:str,by="FR"):
    """This functioni score a nanobody sequence using T20 score. The score can be based on FR or full VH

    Args:
        file (str): The fasta file containiing nanobody sequences
        dest_dir (str): The output directory
        by (str, optional): Score based on "FR" or "VH". Defaults to "FR".

    Returns:
        score (float): the T20 score
    """
    queried_seq = None
    blastdb = None
    if by =="FR":
        queried_seq = numbered_seq.get_fr()
        blastdb = params.T20_HUM_FR_BLAST_DB
    elif by =="VH":
        queried_seq = numbered_seq.get_vh()
        blastdb = params.T20_HUM_VH_BLAST_DB
    else:
        raise ValueError(f'by={by} was not in ("FR","VH")')

    blastp_output = run_blastp_T20(queried_seq,dest_dir,blastdb)
    processed_blastp_output = parse_blastp(blastp_output)
    os.remove(blastp_output)
    df = pd.read_csv(processed_blastp_output)
    os.remove(processed_blastp_output)

    return df['% Identity'].values[0]

if __name__ =="__main__":
    test_file = "./test/Nb21.fasta"
    score_fasta(test_file)


    





