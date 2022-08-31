from Bio import SeqIO
import tempfile
def get_SeqRecord_from_fasta(file:str):
    seqrecord = None
    with open(file) as f:
        seqrecord = next(SeqIO.parse(file,"fasta"))
    return seqrecord

def create_output_folder():
    dirpath = tempfile.mkdtemp(prefix="result_")
    return dirpath
