import numpy as np
import subprocess
from NbHumanization.Scheme import SCHEMES
from NbHumanization import params
from NbHumanization.Sequence import NumberedSequence
import os
import json
import pandas as pd
import logging
logging.basicConfig(\
            filename="NbHumanization.log",\
            #stream=sys.stdout,\
            level=logging.INFO,\
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',\
            datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger('humanizer')
logger.setLevel(logging.INFO)

def load_human_ab_profile(file:str):
    """ This function read pre-computed human antibody heavy chain sequence profile

    Args:
        file (str): the pre-computed profile
    Returns:
        hum_ab_profile (dict): the dict of frequencies
            {'1': {'E': 0.5955518945634267,
                    'Q': 0.24233937397034597,
                    'X': 0.10247116968698518,
                    ...
                    }
            '2': {'V': 0.9609555189456342,
                    ...
                    }
            ...
            }
    """
    hum_ab_profile = json.load(open(file))
    return hum_ab_profile
def get_pos_aa(numbered_seq:NumberedSequence):
    Pos = []
    AA = []
    for p , a in zip(numbered_seq.numbering,numbered_seq.gapped_seq):
        if a =="-" or p =="113":
            continue
        Pos.append(p)
        AA.append(a)
    return Pos,AA

def parse_anarci_file(file:str):
    """This read numbered sequence from anarci output

    Args:
        file (str): The fasta file containing nanobody sequence
    Returns:
        Pos (np.array): the array of amino acid position idices
        AA (np.array): the array of amino acids
    """
    Pos =[]
    AA = []
    with open(file) as f:
        for line in f:
            if line[0]!="H":
                continue
            items = line.rstrip().split(" ")
            if items[-1] =='-':
                continue
            if items[1]=='113':
                #drop last residue
                break
            Pos.append("".join(items[1:-1]))
            AA.append(items[-1])
    #Pos = np.array(Pos)
    #AA = np.array(AA)
    return (Pos,AA)

def getProfileForNb(Pos:np.array,AA:np.array,hum_profile:dict,freq_threshold=0.1):
    """This function evaluate the residue prevalence in human antibodies for the input nanobody

    Args:
        Pos (np.array): The array of amino acids indices 
        AA (np.array): The array of amino acids
        hum_profile (dict): The human antibody heavy chain profile
        freq_threshold (float, optional): The threshold to determine if a residue need to be humanized. Defaults to 0.1.

    Returns:
        [type]: [description]
    """
    scheme = SCHEMES[params.ANNOTATION_DEFINITION][params.ANNOTATION_SCHEME]

    cdr_pos = Pos[Pos.index(str(scheme[0])):Pos.index(str(scheme[1]))+1]\
         + Pos[Pos.index(str(scheme[2])):Pos.index(str(scheme[3]))+1]\
              + Pos[Pos.index(str(scheme[4])):Pos.index(str(scheme[5]))+1]
    freqs = []
    labels = []
    substitutions = []
    for idx,(pos,aa) in enumerate(zip(Pos,AA)):
        if pos in hum_profile.keys() and aa in hum_profile[pos].keys():
            freqs.append(float(hum_profile[pos][aa]))
        else:
            freqs.append(np.min(list(hum_profile[pos].values())))
        if pos in cdr_pos:
            labels.append(-1)
            substitutions.append("")
            
        else:
            if freqs[idx]==None or freqs[idx] >=freq_threshold:
                labels.append(0)
                substitutions.append("")
            else:
                labels.append(1)
                substitution = next(iter(hum_profile[pos].keys()))
                substitutions.append(substitution)
                        
    return (freqs,labels,substitutions) 
def run_protinter(pdb_file:str):
    """ This functioni generate intra-interaction predictions

    Args:
        pdb_file (str): The structure file for the nanobody
    """
    dest_dir = os.path.dirname(pdb_file)
    #basename = os.path.basename(pdb_file).split(".")[0]
    os.chdir(dest_dir)
    ionic_output = os.path.join(dest_dir,f"result_ionic.csv")
    aroaro_output = os.path.join(dest_dir,f"result_aroaro.csv")
    catpi_output = os.path.join(dest_dir,f"result_cationpi.csv")
    subprocess.check_call([params.PROTINTER_EXEC,pdb_file,'-csv','--ionic'], stdout=subprocess.DEVNULL) 
    subprocess.check_call([params.PROTINTER_EXEC,pdb_file,'-csv','--catpi'], stdout=subprocess.DEVNULL) 
    subprocess.check_call([params.PROTINTER_EXEC,pdb_file,'-csv','--aroaro'], stdout=subprocess.DEVNULL)
    logger.info(os.listdir(dest_dir))
    return ionic_output,aroaro_output,catpi_output


def parseIntraInteractions(*files):
    """This function parse interactions in nanobody by Protinter.
    Args:
    *file: a list file for ionic, cation-pi, pi-pi interactions

    Returns:
        interacting_res (set): the resnums of interacting residues
    """
    interacting_res = set()
    for file in files:
        interaction = pd.read_csv(file)
        for _pos in interaction[' idRES1 ']:
            interacting_res.add(int(_pos))
        for _pos in interaction[' idRES2 ']:
            interacting_res.add(int(_pos))
    return interacting_res

def humanize(numbered_seq:NumberedSequence,pdb_file:str,freq_threshold=0.1):
    """This function generate humanized nb sequence

    Args:
        numbered_file (str): The output from ANARCI
        pdb_file (str): The pdb file for the nanobody
        freq_threshold (float, optional): The threshold to determine if a residue need to be humanized . Defaults to 0.1.

    Returns:
        original_seq (str): the original sequence
        df(pd.DataFrame): The dataframe of the humanization data
        humanized_seq (str): the humanized sequence
    """
    hum_ab_profile = load_human_ab_profile(params.HUM_AB_PROFILE)
    pos,aa = get_pos_aa(numbered_seq)
    #pos,aa = parse_anarci_file(numbered_file)
    freqs,labels,substitutions = getProfileForNb(pos,aa,hum_ab_profile,freq_threshold)
    interaction_files = run_protinter(pdb_file)
    interacting_res = parseIntraInteractions(*interaction_files)
    df = pd.DataFrame({"pos":pos,"aa":aa,"freq":freqs,"label":labels,"sub":substitutions})
    df['idx'] = np.arange(1,df.shape[0]+1)
    df['intra'] = df['idx'].apply(lambda x: int(x in interacting_res))
    df['humanized'] = df.apply(lambda x: 1 if x['label']==1 and x['intra']==0 else 0,1)
    df['humanized_aa'] = df.apply(lambda x: x['sub'] if x['humanized']==1 else x['aa'],1)
    original_seq = "".join(df['aa'].values)
    humanized_seq = "".join(df['humanized_aa'].values)
    df.to_csv(os.path.join(os.path.dirname(pdb_file),"humanized_data.csv"))
    return original_seq,df,humanized_seq
