# Server is now avaialble at [llamanade.app](www.llamande.app). 

Send a request to zhe[dot]sang[at]icahn[dot]mssm[dot]edu if you encounter any difficulty.



# Llamanade: an open-source computational pipeline for robust nanobody humanization 

Zhe Sang,Yufei Xiang, Ivet Bahar and Yi Shi

University of Pittsburgh
## Abstract 

Nanobodies (Nbs) have recently emerged as a promising class of antibody fragments for biomedical and therapeutic applications. Despite having marked physicochemical properties, Nbs are derived from camelids and may require “humanization” to improve translational potentials for clinical trials. Here we have systematically analyzed the sequence and structural properties of Nbs based on NGS (next-generation sequencing) databases and high-resolution structures. Our analysis reveals substantial framework diversities and underscores the key differences between Nbs and human Immunoglobulin G (IgG) antibodies. We identified conserved residues that may contribute to enhanced solubility, structural stability, and antigen-binding, providing insights into Nb humanization. Based on big data analysis, we developed “Llamanade'', a user-friendly, open-source to facilitate rational humanization of Nbs. Using Nb sequence as input, Llamanade provides information on the sequence features, model structures, and optimizes solutions to humanize Nbs. The full analysis for a given Nb takes less than a minute on a local computer. To demonstrate the robustness of this tool, we applied it to successfully humanize a cohort of structurally diverse and highly potent SARS-CoV-2 neutralizing Nbs. Llamanade is freely available and will be easily accessible on a web server to support the development of a rapidly expanding repertoire of therapeutic Nbs into safe and effective trials.

for citation, please cite our paper: https://doi.org/10.1016/j.str.2021.11.006
<p align="left"><img src="https://ars.els-cdn.com/content/image/1-s2.0-S0969212621004184-fx1.jpg" height="500"/></p>

How to run Llamanade locally:

    1. Clone the git repository : git clone "https://github.com/sangzhe/Llamanade.git"
    2. Make sure you have the following libraries installed in your environment:
            - HMMER(http://hmmer.org)(required by ANARCI)
            - ANARCI(https://github.com/oxpig/ANARCI)
            - NCBI Blastp(https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST) or (sudo apt-get install ncbi-blast+)
            - Bio (pip install Bio)
            - ProDy(http://prody.csb.pitt.edu/)(pip install prody)
            - Protinter(https://github.com/maxibor/protinter)
            - Modeller (requires license - https://salilab.org/modeller/)
    3. Configure NanoNet
            - NanoNet is included here, please install required library accordingly.
            - NanoNet only predict CA atoms, sidechains will be constructed using pulchr. The source code is also included under NanoNet folder, please compile first and move the executable under the NanoNet folder
            - Test if NanoNet works
    4. Unzip recourses.zip
    5. Edit params.py
            - Varaibles in this file refer to the executables or resources, please edit accordingly especially the path for LLAMANADE. 


Run Llamanade
```bash
usage: Nb Humanization [-h] [--fa FA] [--pdb PDB] [--chain CHAIN] [--modeling {NanoNet,Modeller}]

Program for comprehensive nanobody humanization

optional arguments:
  -h, --help            show this help message and exit
  --fa FA, -f FA        sequence file of a nanobody in fasta format
  --pdb PDB, -p PDB     structure file of a nanobody in pdb format
  --chain CHAIN, -c CHAIN
                        chain of nanobody in structure file
  --modeling {NanoNet,Modeller}, -m {NanoNet,Modeller}
                        structural modeling tools

```

Run with a sequence in FASTA
```bash
(base) zhesang@zhe-Alienware:~/Llamanade$ python NbHumanization_main.py --fa test/Nb21.fa
Initial coordinates will be preserved.
reconstruct
/home/zhesang/Llamanade
NanoNet ended successfully, models are located in directory:'/tmp/result_zvdp6c7w/NanoNet', total time : 1.934037835104391.
1 processed
1 processed
QVQLVESGGGLVQAGGSLRLSCAVSGLGAHRVGWFRRAPGKEREFVAAIGANGGNTNYLDSVKGRFTISRDNAKNTIYLQMNSLKPQDTAVYYCAARDIETAEYTYWGQGTQVTVS:82.0931
QVQLVESGGGLVQPGGSLRLSCAASGLGAHRVGWVRRAPGKGLEFVAAIGANGGNTNYADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCARRDIETAEYTYWGQGTLVTVS:94.7688
result saved to /home/zhesang/Llamanade/test/result_zvdp6c7w
```

Run with a structure in PDB
```bash
(base) zhesang@zhe-Alienware:~/Llamanade$ python NbHumanization_main.py --pdb test/Nb20.pdb --chain C
@> 852 atoms and 1 coordinate set(s) were parsed in 0.02s.
1 processed
1 processed
QVQLVESGGGLVQAGGSLRLSCAVSGAGAHRVGWFRRAPGKEREFVAAIGASGGMTNYLDSVKGRFTISRDNAKNTIYLQMNSLKPQDTAVYYCAARDIETAEYIYWGQGTQVTVS:82.1499
QVQLVESGGGLVQPGGSLRLSCAASGAGAHRVGWFRQAPGKEREFVAAIGASGGMTNYADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCARRDIETAEYIYWGQGTLVTVS:92.4397
result saved to /home/zhesang/Llamanade/test/result_g5f9lf2f
```
