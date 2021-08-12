wget --no-check-certificate https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar -xzvf muscle3.8.31_i86linux64.tar.gz
cp muscle3.8.31_i86linux64 /usr/bin/muscle

wget http://eddylab.org/software/hmmer/hmmer.tar.gz 
tar zxf hmmer.tar.gz
cd hmmer-3.3.2
./configure
make
make check
make install
cd ..
rm hmmer.tar.gz

