import subprocess
import argparse
import os
import sys
import logging
from Bio.SeqRecord import SeqRecord
from Sequence import NumberedSequence
from Scheme import NUMBERING,DEFINITIONS
import params

PDB_PATH="/opt/NbModeling/resources/modeller_data/PDBs"
BLASTP_DB_PATH = "/opt/NbModeling/resources/modeller_data/Modeller_VH"

logging.basicConfig(\
            filename="NbHumanization.log",\
            #stream=sys.stdout,\
            level=logging.INFO,\
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',\
            datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger('Annotation')
logger.setLevel(logging.INFO)

def run_anarci_numbering(sequence:SeqRecord,scheme:str,task_dir:str):
    """ This function calls ANARCI to number the sequence

    Args:
        sequence (SeqRecord): Antibody/Nanobody heavy chain sequence
        scheme (str): IMGT/AHo/Kabat/Chothia
        task_dir (str):  path to task folder

    Returns:
        outfile (str): path of output file
    """


    outfile = os.path.join(task_dir, f"{sequence.id}.numbering")
    subprocess.check_call(['ANARCI','-i',str(sequence.seq),'--r','heavy','--ncpu','4','--scheme',scheme,'--outfile',outfile,'--assign_germline']) 
    logger.info(f"Finished runnning ANARCI using scheme:{scheme},output in {outfile}")
    return outfile



def run_blastp(query_file:str,blast_db:str,max_target_seqs:int,blastp_outfmt:str,task_dir:str):
    """ This function will call blastp to conduct search

    Args:
        query_file (str): the file path to antibody/Nanobody heavy chain sequence
        blast_db (str): path to blast database
        max_target_seqs (int): maximum target sequences to retrieve
        blastp_outfmt (str): format of blastp output
        task_dir (str): path to task folder

    Returns:
        str: the output file path of blastp result
    """
    basename = os.path.basename(query_file).split(".")[0]
    outfile = os.path.join(task_dir,f"{basename}.blastp")
    subprocess.check_call([params.BLAST_EXEC,'-db',blast_db,'-query',query_file,'-max_target_seqs',str(max_target_seqs),'-outfmt',blastp_outfmt,'-out',outfile]) 
    return outfile

def query_template_by_blastp(query_file:str,task_dir:str,exclude_identical:False):
    """ This function process blastp result to retrieve template id from PDB

    Args:
        query_file (str): the file path to antibody/Nanobody heavy chain sequence
        blast_db (str): path to blast database
        task_dir (str): path to task folder
        exclude_identical (False): whether to exclude identical seuquence found in the databse

    Returns:
        (tuple): (template_id,chain_id)
    """
    max_target_seqs = 1
    if exclude_identical:
        max_target_seqs = 2
    max_target_seqs = 5

    #TODO: check the output format
    blastp_outfmt = "6 sseqid pident sstart send"
    
    blast_result_file = run_blastp(query_file,BLASTP_DB_PATH,max_target_seqs,blastp_outfmt,task_dir)
    logger.info("Finished running blastp")

    template_id = None
    chain_id=None
    with open(blast_result_file) as f:
        ''' Example output
        1QDO:A  100.0   1   119
        5LHR:B  96.0    1   121
        5IMK:B  95.181  1   123
        '''
        hits = f.read().splitlines()
        hits = [i.split("\t") for i in hits]
        hits.sort(key=lambda x:x[1],reverse=True)
        print(hits)
        
        if len(hits) ==0:
            raise RuntimeError(f"Not hit was found by blastp")
        if exclude_identical:
                template_id,chain_id = hits[1][0].split(":")
                template_start = int(hits[1][2])
                template_end = int(hits[1][3])
        else:
            template_id,chain_id = hits[0][0].split(":")
            template_start = int(hits[0][2])
            template_end = int(hits[0][3])
        


    return (template_id,chain_id,template_start,template_end)

def annotate(sequence:SeqRecord,scheme:str,definition_method:str,task_dir:str):
    """This function number the input VH sequence according desinated scheme

    Args:
        sequence (SeqRecord): [description]
        scheme (str): imgt/aho/kabat/chothia
        definition_method (str): IMGT/AHo/AbM/Kabat/Chothia
        task_dir (str):  path to task folder

    Returns:
        [type]: [description]
    """
    if not (scheme in NUMBERING):
        scheme = 'imgt'
        logger.warn(f"{scheme} is not valid, use IMGT in default")
    if not(definition_method in DEFINITIONS):
        definition_method = 'IMGT'
        logger.warn(f"{definition_method} is not valid, use IMGT in default")

    result_file = run_anarci_numbering(sequence,scheme,task_dir)
    numbered_seq = NumberedSequence(seq = str(sequence.seq),id=sequence.id)
    numbered_seq.set_definition(definition_method)
    try:
        with open(result_file) as f:
            content = f.read().rstrip()
            numbered_seq.parse_anarci_numbering(content)
    except:
        logger.error(f"error occurred while processing {sequence.seq}")
        return None
    
    return numbered_seq


def search_template(numbered_seq:NumberedSequence,task_dir:str,exclude_identical:False):
    """This function search the most similar template structure for input sequence

    Args:
        sequence (NumberedSequence): Antibody/Nanobody heavy chain sequence
        task_id (str): path to task folder
        exclude_identical (False): whether to exclude identical seuquence found in the databse
    Returns:
        (str): template id in PDB code
    """
    query_seq = numbered_seq.get_vh()
    '''
    if len(numbered_seq.CDR3) > 30:
        query_seq = numbered_seq.get_vh()
    else:
        query_seq = numbered_seq.get_pseudovh()
    '''
    print(f"Input Seq:{numbered_seq.seq}")
    print(f"Query with :{query_seq}")
    query_fname = os.path.join(task_dir,numbered_seq.id+"_query.fasta")
    with open(query_fname,'w+') as out:
        out.write(f">{numbered_seq.id}\n{query_seq}\n")
    template_id,chain_id,template_start,template_end = query_template_by_blastp(query_fname,task_dir,exclude_identical)
    template_file_path = os.path.join(PDB_PATH,f"{template_id}_{chain_id}.pdb")
    return (template_id,chain_id,template_start,template_end,template_file_path)

    






