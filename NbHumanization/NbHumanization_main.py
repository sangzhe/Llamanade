import argparse
import shutil
from humanizer import humanize
from StructUtils import model_by_Modeller,model_by_NanoNet,extract_heavy_chain_from_structure
import os
import logging
import traceback
from IOUtils import create_output_folder
from T20Scorer import score,score_plain_seq
logging.basicConfig(\
            filename="NbHumanization.log",\
            #stream=sys.stdout,\
            level=logging.INFO,\
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',\
            datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger('Main')
logger.setLevel(logging.INFO)

def main():
    parser = argparse.ArgumentParser(prog="Nb Humanization",description='Program for comprehensive nanobody humanization')
    parser.add_argument("--fa","-f",help="sequence file of a nanobody in fasta format")
    parser.add_argument("--pdb","-p",help="structure file of a nanobody in pdb format")
    parser.add_argument("--chain","-c",help="chain of nanobody in structure file")
    parser.add_argument("--modeling","-m",default="NanoNet",choices=['NanoNet', 'Modeller'],help="structural modeling tools")

    args = parser.parse_args()
    do_modeling=False

    fa_file = args.fa
    pdb_file = args.pdb
    chain_id = args.chain
    model_tool = args.modeling
    

    dest_dir = create_output_folder()
    try:
        numbered_seq = None
        nb_struct_file=None

        if fa_file!=None:
            if model_tool=="NanoNet":
                numbered_seq,nb_struct_file = model_by_NanoNet(fa_file,dest_dir)
            else:
                numbered_seq,nb_struct_file = model_by_Modeller(fa_file,dest_dir)
            
        elif pdb_file!= None and chain_id!=None:
            numbered_seq,nb_struct_file  = extract_heavy_chain_from_structure(pdb_file,chain_id,dest_dir)
        
        else:
            parser.print_help()
            raise(OSError("arguments are not complete. Provide either a [fasta] or [pdb and chain] information"))
        logger.info(numbered_seq)
        logger.info(nb_struct_file)
        #basename = os.path.basename(nb_struct_file).split(".")[0]
        #anarci_file = f"{basename}.numbering"
        original_seq,df,humanized_seq = humanize(numbered_seq,nb_struct_file)

        original_score = score(numbered_seq,dest_dir)
        humanized_score = score_plain_seq(humanized_seq,dest_dir)

        logger.info(f"{original_seq}:{original_score}")
        logger.info(f"{humanized_seq}:{humanized_score}")

        print(f"{original_seq}:{original_score}")
        print(f"{humanized_seq}:{humanized_score}")

    



        
    except Exception as e:
        traceback.print_exc()
    finally:
        print(dest_dir)
        #shutil.rmtree(dest_dir)


if __name__ == "__main__":
    main()