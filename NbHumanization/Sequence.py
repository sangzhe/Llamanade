from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import defaultdict
import numpy as np
import re
import os
from NbHumanization.Scheme import SCHEMES,NUMBERING

class NumberedSequence(object):
    def __init__(self, seq:str,id:str):
        self.id = id
        self.seq = seq
    

    def set_VH_range(self,start,end):
        self.vh_start=start
        self.vh_end = end
    def set_numbering_scheme(self,scheme):
        self.numbering_scheme = scheme 
    def set_definition(self,definition):
        self.definition = definition
    def set_species(self ,species:str):
        """ Set species according to ANARCI result

        Args:
            species (str): human/mouse/alpaca/......
        """
        self.species = species

    def set_IGHV_germline(self,germline):
        self.ighv_germline = germline
    def set_IGHJ_germline(self,germline):
        self.ighj_germline = germline
    def set_IGHV_germline_identity(self,germline_identity):
        self.ighv_germline_identity = germline_identity
    def set_IGHJ_germline_identity(self,germline_identity):
        self.ighj_germline_identity = germline_identity

    def set_numbered_seqeunce(self, numbering:list,gapped_sequence:str):
        self.numbering = numbering
        self.gapped_seq = gapped_sequence

    # def __getitem__(self, index):

    #     if isinstance(index,int):
    #         real_start = self.numbering_index_dict[str(index)]
    #         real_stop = self.numbering_index_dict[str(index+1)]
    #         return self.gapped_seq[real_start:real_stop]
    #     elif isinstance(index,slice):
    #         if index.start == None or index.start < 1:
    #             start = 1
    #         else:
    #             start = index.start

    #         if index.stop == None or index.stop > int(self.numbering[-1]):
    #             stop = int(self.numbering[-1])
    #         else:
    #             stop = index.stop

    #         real_start = self.numbering_index_dict[str(start)]
    #         real_stop = self.numbering_index_dict[str(stop)]
    #         return self.gapped_seq[real_start:real_stop]
    #     else:
    #         raise ValueError(f"Can not slice by {index}")

    def parse_anarci_numbering(self,result:str):
        """ This function parses result from ANARCI
        Example:
        # Input sequence
        # ANARCI numbered
        # Domain 1 of 1
        # Most significant HMM hit
        #|species|chain_type|e-value|score|seqstart_index|seqend_index|
        #|alpaca|H|8.6e-64|204.6|2|127|
        # Most sequence-identical germlines
        #|species|v_gene|v_identity|j_gene|j_identity|
        #|alpaca|IGHV3-3*01|0.93|IGHJ4*01|1.00|
        # Scheme = aho
        H 1       D
        H 2       V
        H 3       Q
        H 4       L
        H 5       V
        .
        .
        .
        H 146     T
        H 147     V
        H 148     S
        H 149     S

        Args:
            result (str): content of result
        """
        _1,_2,_3,_4,_5,_6,_7,_8,_9,_10,*_11,_end = result.splitlines()
            
        blank, _hmm_species, _chain_type, _e_value, _score, _seqstart_index, _seqend_index,b = _6.split("|")
        self.set_VH_range(int(_seqstart_index),int(_seqend_index))

        blank,_germline_species,_v_gene,_v_identity,_j_gene,_j_identity,_t = _9.split("|")

        self.set_species(_germline_species)
        self.set_IGHV_germline(_v_gene)
        self.set_IGHV_germline_identity(_v_identity)
        self.set_IGHJ_germline(_j_gene)
        self.set_IGHJ_germline_identity(_j_identity)

        _scheme = _10.split("=")[1].strip()
        self.set_numbering_scheme(_scheme)

        #prefix = seq[:int(_seqstart_index)]
        #suffix = seq[int(_seqend_index)+1:]
        numbering = []
        resi = []
        _seq = [""] * 150

        if len(_11) <100:
            raise ValueError(f"Poor sequence quality")
        for line in _11:           
            aa = line[-1]
            _num = line[1:-1].replace(" ","")
            numbering.append(_num)
            resi.append(aa)

            if aa!='-':
                _resn = int(re.sub(r"\D","",_num))
                _seq[_resn] = _seq[_resn]+aa

        gapped_seq = "".join(resi)
        self.determine_FR_CDR(_seq)

        self.set_numbered_seqeunce(numbering,gapped_seq)
    def determine_FR_CDR(self,concated_seq):
        delimiters = SCHEMES[self.definition][self.numbering_scheme]
        self.FR1 ="".join(concated_seq[:delimiters[0]])
        self.CDR1="".join(concated_seq[delimiters[0]:delimiters[1]])
        self.FR2 ="".join(concated_seq[delimiters[1]:delimiters[2]])
        self.CDR2="".join(concated_seq[delimiters[2]:delimiters[3]])
        self.FR3 ="".join(concated_seq[delimiters[3]:delimiters[4]])
        self.CDR3="".join(concated_seq[delimiters[4]:delimiters[5]])
        self.FR4 ="".join(concated_seq[delimiters[5]:])
    def get_fr(self):
        return f"{self.FR1}{self.FR2}{self.FR3}{self.FR4}"

    def get_vh(self):
        return f"{self.FR1}{self.CDR1}{self.FR2}{self.CDR2}{self.FR3}{self.CDR3}{self.FR4}"
    def get_pseudovh(self):
        return f'{self.FR1}{"".join(["X"]*len(self.CDR1))}{self.FR2}{"".join(["X"]*len(self.CDR2))}{self.FR3}{"".join(["X"]*len(self.CDR3))}{self.FR4}'

    # def get_fr(self):
    #     if self.scheme == "aho":
    #         fr1 = self.__getitem__(slice(0,25))
    #         fr2 = self.__getitem__(slice(41,58))
    #         fr3 = self.__getitem__(slice(78,109))
    #         fr4 = self.__getitem__(slice(138,len(self.numbering_index_dict)))

    #     elif self.scheme == "imgt":
    #         fr1 = self.__getitem__(slice(0,27))
    #         fr2 = self.__getitem__(slice(39,56))
    #         fr3 = self.__getitem__(slice(66,105))
    #         fr4 = self.__getitem__(slice(118,len(self.numbering_index_dict)))
    #     else:
    #         raise ValueError(f"{self.scheme} can not be used temporarily")
    #     concat_seq =  f"{fr1}{fr2}{fr3}{fr4}"
    #     concat_seq = "".join([i for i in concat_seq if i.isalpha()])
    #     return concat_seq

    # def get_vh(self):
    #     if self.scheme == "aho":
    #         fr1 = self.__getitem__(slice(1,25))
    #         cdr1 = self.__getitem__(slice(25,41))
    #         fr2 = self.__getitem__(slice(41,58))
    #         cdr2 = self.__getitem__(slice(58,78))
    #         fr3 = self.__getitem__(slice(78,109))
    #         cdr3 = self.__getitem__(slice(109,138))
    #         fr4 = self.__getitem__(slice(138,len(self.numbering_index_dict)))
    #     elif self.scheme == "imgt":
    #         fr1 = self.__getitem__(slice(1,27))
    #         cdr1 = self.__getitem__(slice(27,39))
    #         fr2 = self.__getitem__(slice(39,56))
    #         cdr2 = self.__getitem__(slice(56,66))
    #         fr3 = self.__getitem__(slice(66,105))
    #         cdr3 = self.__getitem__(slice(105,118))
    #         fr4 = self.__getitem__(slice(118,len(self.numbering_index_dict)))
    #     else:
    #         raise ValueError(f"{scheme} can not be used temporarily")
    #     concat_seq =  f"{fr1}{cdr1}{fr2}{cdr2}{fr3}{cdr3}{fr4}"
    #     concat_seq = "".join([i for i in concat_seq if i.isalpha()])
    #     return concat_seq

    def to_pir(self,target_folder):
        fname = os.path.join(target_folder,str(self.id)+".pir")
        out = open(fname,'w+')
        out.write(f">P1;{self.id}\n")
        out.write(f"sequence:{self.id}:::::::0.00: 0.00\n")
        degapped_seq = self.get_vh()
        out.write(f"{degapped_seq}*")
        out.close()
        return fname
    def __str__(self):
        r_str = f'raw seq:{self.seq}\nid:{self.id}\n'
        r_str += f'vh_start:{self.vh_start}\nvh_end:{self.vh_end}\n'
        r_str += f'Lead:{self.seq[:self.vh_start]}\n'
        r_str += f'FR1:{self.FR1}\nFR2:{self.FR2}\nFR3:{self.FR3}\nFR4:{self.FR4}\n'
        r_str += f'CDR1:{self.CDR1}\nCDR2:{self.CDR2}\nCDR3:{self.CDR3}\n'
        r_str += f'Tail:{self.seq[self.vh_end+1:]}'
        return r_str

    def __repr__(self):
        return self.__str__()

    def serialize(self,target_folder):
        fname = os.path.join(target_folder,str(self.id)+".txt")
        out = open(fname,'w+')
        out.write(str(self))
        out.close()
