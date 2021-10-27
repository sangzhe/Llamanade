wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-src.tar.gz
tar -xzvf ncbi-blast-2.11.0+-src.tar.gz
rm ncbi-blast-2.11.0+-src.tar.gz
cd ncbi-blast-2.11.0+-src
cd c++
./configure
cd ReleaseMT/build
make all_r
