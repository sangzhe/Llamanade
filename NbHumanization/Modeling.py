from modeller import *
from Sequence import NumberedSequence
from modeller.automodel import *
from modeller.parallel import *
from modeller import soap_protein_od
import subprocess
import logging
#log.verbose()
import os
import sys
import params


logging.basicConfig(\
            filename="NbHumanization.log",\
            #stream=sys.stdout,\
            level=logging.INFO,\
            format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',\
            datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger('Modelling')
logger.setLevel(logging.INFO)
log.none()
class AbLoop(loopmodel):
    def add_cdr3_position(self,cdr3_start:int,cdr3_end:int):
        self.cdr3_start = str(cdr3_start)+":A"
        self.cdr3_end = str(cdr3_end)+":A"
        logger.info(f'CDR3 start:{cdr3_start},CDR3 end:{cdr3_end}')

    def select_loop_atoms(self):
        return selection(self.residue_range(self.cdr3_start,self.cdr3_end))


class ComparativeModellingParameters(object):
    def __init__(self,task_dir:str):
        self.task_dir = task_dir
    def add_template(self,pdb_id:str,chain_id:str,template_start:int,template_end:int,pdb_file_path:str):
        self.template_pdb_id = pdb_id
        self.template_chain_id = chain_id
        self.template_start=template_start
        self.template_end = template_end
        self.template_pdb_file_path = pdb_file_path
    def add_target(self,seq:NumberedSequence,seq_file_path:str):
        self.target_seq = seq
        self.target_seq_file_path = seq_file_path
    def get_template_file(self):
        return self.template_pdb_file_path
    def get_template_model_segment(self):
        return (f"FIRST:{self.template_chain_id}",f"LAST:{self.template_chain_id}")
        #return (f"{self.template_start}:H",f"{self.template_end}:H")
    def get_template_align_codes(self):
        return f"{self.template_pdb_id.lower()}{self.template_chain_id}"
        #return f"{self.template_pdb_id.lower()}H"
    def get_target_file(self):
        return self.target_seq_file_path
    def get_target_align_codes(self):
        return self.target_seq.id

class NanoNetModeling(object):   
    def __init__(self,task_dir):
        self.models = None
        self.best_model = None
        full_dest_dir = os.path.join(task_dir,"NanoNet")
        os.mkdir(full_dest_dir)
        self.task_dir = full_dest_dir
        
    def model(self,filename):
        cmd = ["python3",f"{params.NANONET_EXEC}/NanoNet.py","-n",f"{params.NANONET_EXEC}/NanoNet","-s",filename,"-o",self.task_dir,"-p",f'{params.NANONET_EXEC}/pulchra']
        subprocess.check_call(cmd)
        
        self.best_model = os.path.join(self.task_dir,os.listdir(self.task_dir)[0])


class Modeling(object):
    def __init__(self,parameters:ComparativeModellingParameters):
        self.params = parameters
        self.models = None
        self.best_model = None
        self.environ = environ()
        #self.jobs = job(host='localhost')
        #for i in range(10):
        #    self.jobs.append(local_slave())
        os.chdir(self.params.task_dir)

    def model(self):
        env = self.environ
        aln = alignment(env)
        mdl = model(env, file=self.params.get_template_file(), model_segment=self.params.get_template_model_segment())
        #extracted_template_name = os.path.join(parameters.task_dir,f'{parameters.get_target_align_codes()}-{parameters.get_template_align_codes()}-template.pdb')
        #s = selection()
        #s.add(mdl.chains[0].residues[:parameters.template_end+5])
        #s.write(extracted_template_name)
        aln.append_model(mdl, align_codes=self.params.get_template_align_codes(), atom_files=self.params.get_template_file())
        #aln.append_model(mdl, align_codes=parameters.get_template_align_codes(), atom_files=extracted_template_name)
        aln.append(file=self.params.get_target_file(), align_codes=self.params.get_target_align_codes())
        aln.align()

        ali_output = os.path.join(self.params.task_dir,f"{self.params.get_target_align_codes()}-{self.params.get_template_align_codes()}.ali")
        #aln_output = f"{parameters.get_target_align_codes()}-{parameters.get_template_align_codes()}.pap"

        aln.write(file=ali_output, alignment_format='PIR')
        #aln.write(file=aln_output, alignment_format='PAP')
        
        a = automodel(env, alnfile=ali_output,
                    knowns=self.params.get_template_align_codes(), sequence=self.params.get_target_align_codes(),
                    assess_methods=(assess.DOPE,
                                    #soap_protein_od.Scorer(),
                                    assess.GA341))
        
        a.blank_single_chain=False
        a.starting_model = 1
        a.ending_model = 10
        #a.use_parallel_job(self.jobs)
        a.make()
        key = 'DOPE score'
        # Get a list of all successfully built models from a.outputs
        ok_models = [x for x in a.outputs if x['failure'] is None]
        ok_models.sort(key=lambda a: a[key])


        # Get top model
        self.models = [_mdl['name'] for _mdl in ok_models]
        self.best_model = ok_models[0]['name']
        logger.info(f'Best Model:{self.best_model}')
    def loop_refinement(self):
        #
        # First try FREAD
        # if failed
        # ab-inito loop modelling
        if not self.loop_refinement_by_pyfread():
            self.loop_refinement_by_ab_initio()
        

    def loop_refinement_by_pyfread(self):
        logger.info("FREAD loop prediction")
        cdr3 = self.params.target_seq.CDR3
        cdr3_start = self.params.target_seq.seq.find(cdr3)
        subprocess.check_call(['pyfread','/home/shilab/FREAD/db_h3',self.best_model,str(cdr3_start),cdr3]) 
        pyfread_summary_size = os.path.getsize("summary.table")
        if pyfread_summary_size == 0:
            return False
        else:
            loop_refined_result = open("summary.table").read()
            best_loop_refined_model = loop_refined_result.splitlines()[0].split("\t")[0]+".model.pdb"
            logger.info(f'best loop refined model:{best_loop_refined_model}')
            self.loop_refined_model = best_loop_refined_model
            return True

    def loop_refinement_by_ab_initio(self):
        logger.info("Ab initio modeling")
        cdr3 = self.params.target_seq.CDR3
        cdr3_start = self.params.target_seq.seq.find(cdr3)
        cdr3_end = cdr3_start + len(cdr3)
        env = self.environ
        m = AbLoop(env,
           inimodel=self.best_model, # initial model of the target
           sequence=self.params.get_target_align_codes(),
           loop_assess_methods=(assess.DOPE
                            #soap_protein_od.Scorer()
                            ))# code of the target
        m.add_cdr3_position(cdr3_start=cdr3_start,cdr3_end=cdr3_end)          

        m.loop.starting_model= 1           # index of the first loop model 
        m.loop.ending_model  = 2          # index of the last loop model
        m.loop.md_level = refine.fast # loop refinement method; this yields
                                   # models quickly but of low quality;
                                   # use refine.slow for better models

        m.make()

        key = 'DOPE score'
        # Get a list of all successfully built models from a.outputs
        ok_models = [x for x in m.loop.outputs if x['failure'] is None]
        ok_models.sort(key=lambda a: a[key])


        # Get top model
        self.models = [_mdl['name'] for _mdl in ok_models]
        self.best_model = ok_models[0]['name']
        logger.info(f'Best Model:{self.best_model}')

