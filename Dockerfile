FROM ubuntu:18.04

# install necessary dependencies
RUN apt-get update && \
    apt-get install -y wget python3-dev
RUN apt-get install -y build-essential
RUN apt-get install -y git
WORKDIR /opt

#Install Modeller
RUN mkdir modeller
WORKDIR /opt/modeller
RUN wget --no-check-certificate https://salilab.org/modeller/9.25/modeller_9.25-1_amd64.deb
RUN env KEY_MODELLER=MODELIRANJE dpkg -i /opt/modeller/modeller_9.25-1_amd64.deb

#Install muscle
WORKDIR /opt
COPY install_ANARCI.sh /opt
RUN sh install_ANARCI.sh
COPY install_ANARCI2.sh /opt
RUN sh install_ANARCI2.sh

COPY install_blast.sh /opt
RUN sh install_blast.sh

RUN mkdir /NbModeling
COPY NbHumanization /opt/NbHumanization

RUN mkdir /data
RUN mkdir /data/BlastDB
COPY BlastDB /data/BlastDB

RUN apt-get install -y vim
RUN pip3 install pandas

RUN mkdir /opt/protinter
COPY protinter /opt/protinter

RUN pip3 install prody

#WORKDIR /data
#CMD ["ANARCI","-h"]
