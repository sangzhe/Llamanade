from prody import parsePDB,writePDB
from Bio.PDB.Polypeptide import three_to_one 
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Annotation import annotate,search_template
from Modeling import ComparativeModellingParameters,Modeling
from IOUtils import get_SeqRecord_from_fasta
import params
import os
import logging

logging.basicConfig(\
            filename="NbHumanization.log",\
            #stream=sys.stdout,\
            level=logging.INFO,\
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',\
            datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger('StructUtils')
logger.setLevel(logging.INFO)

def extract_heavy_chain_from_structure(pdb_file:str,chain_id:str,dest_dir:str):
    """This function extract the nanobody heavy chain from a given pdb file and chain id.
    And save the extracted structure to the designated directory
    The lead and tail sequence does not belong to VH will be trimmed.

    Args:
        pdb_file (str): Structure file provided by user
        chain_id (str): chain id specifying the nanobody
        dest_dir (str): output directory
    Returns:
        numbered_seq (NumberedSequence): the numbered sequence object
        output_path (str):  The absolute path to the extracted structure
    """
    pdb_name = os.path.basename(pdb_file).split(".")[0]
    struct = parsePDB(pdb_file,chain=chain_id)
    resname_3codes = struct.ca.getResnames()
    chain_sequence = "".join([three_to_one(aa) for aa in resname_3codes])

    seqrecord = SeqRecord(id=pdb_name,seq=Seq(chain_sequence))
    numbered_seq = annotate(seqrecord,params.ANNOTATION_SCHEME,params.ANNOTATION_DEFINITION,dest_dir)
    vh_start  = numbered_seq.vh_start
    vh_end = numbered_seq.vh_end

    vh_struct = struct.select(f"resnum {vh_start} to {vh_end}")
    output_path = os.path.join(dest_dir,f"VH_{pdb_name}.pdb")
    writePDB(output_path,vh_struct,renumber=True)
    logger.info(f"extracted VH structure to {output_path}")
    
    return numbered_seq,output_path

def model(file:str,dest_dir:str):
    """This function generate homology model for the nanobody

    Args:
        file (str): The fasta file containing nanobody sequence
        dest_dir (str): The output directory
    Returns:
        numbered_seq (NumberedSequence): the numbered sequence object
        output_path (str):  The absolute path to the modeled structure
    """
    seqrecord = get_SeqRecord_from_fasta(file)
    numbered_seq = annotate(seqrecord,params.ANNOTATION_SCHEME,params.ANNOTATION_DEFINITION,dest_dir)

    pir_file_path = numbered_seq.to_pir(dest_dir)
    template_id,chain_id,template_start,template_end,template_file_path = search_template(numbered_seq,dest_dir,True)

    paramters = ComparativeModellingParameters(dest_dir)
    paramters.add_template(template_id,chain_id,template_start,template_end,template_file_path)
    paramters.add_target(numbered_seq,pir_file_path)

    m = Modeling(paramters)
    m.model()
   # _m.loop_refinement()
    return numbered_seq,os.path.join(dest_dir,m.best_model)
    









