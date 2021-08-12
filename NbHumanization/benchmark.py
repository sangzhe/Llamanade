from Annotation import *
from Modeling import ComparativeModellingParameters,Modeling
from Sequence import NumberedSequence
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import shutil
import hashlib

TASKS_DIR="/home/shilab/NanobodyHumanization"

def create_task(parent_dir,plain_sequence):
    seuqence_hash = str(hashlib.md5(plain_sequence.encode('UTF-8')).hexdigest())
    task_dir_path = os.path.join(parent_dir,seuqence_hash)
    if os.path.exists(task_dir_path):
        shutil.rmtree(task_dir_path)
    os.mkdir(task_dir_path)
    return task_dir_path

def test_modeling(parent_dir:str,seq:SeqRecord):

    task_dir = create_task(parent_dir,str(seq.seq))
    numbered_seq = annotate(seq,"martin","AbM",task_dir)
    if numbered_seq == None:
        shutil.rmtree(task_dir)
        return
    pir_file_path = numbered_seq.to_pir(task_dir)
    numbered_seq.serialize(task_dir)
    template_id,chain_id,template_start,template_end,template_file_path = search_template(numbered_seq,task_dir,True)

    paramters = ComparativeModellingParameters(task_dir)
    paramters.add_template(template_id,chain_id,template_start,template_end,template_file_path)
    paramters.add_target(numbered_seq,pir_file_path)
    for k,v in paramters.__dict__.items():
            print(f"{k}\t{v}")

    _m = Modeling(paramters)
    _m.model()
    best_model_name = _m.best_model
    shutil.move(best_model_name,task_dir)
    for f in os.listdir():
        if f.startswith(numbered_seq.id):
            os.remove(f)

def main():
    nb_dir = "/home/shilab/NanobodyHumanization/Nb_sampling_models"
    #human_dir = "/home/shilab/NanobodyHumanization/Human_benchmark"
    Nb_test_files = "/home/shilab/NanobodyHumanization/Nb_sampling_modelling.fasta"
    #Human_test_files = "/home/shilab/NanobodyHumanization/Benchmark_Human.fasta"
    
    for record in SeqIO.parse(Nb_test_files,"fasta"):
       test_modeling(nb_dir,record) 
       print("====================")
    print("All Nbs modelling finished")
    #for record in SeqIO.parse(Human_test_files,"fasta"):
    #   test_modeling(human_dir,record) 
    #   print("====================")

if __name__ == "__main__":
    main()


