FROM ubuntu
MAINTAINER Jeff Cole <coleti16@students.ecu.edu>

# Install all the software needed to run the pipeline
RUN apt-get -qq update
RUN apt-get install -y wget git build-essential cmake unzip curl
RUN apt-get install -qqy python3-setuptools python3-docutils python3-flask
RUN easy_install3 snakemake pip



ARG package
WORKDIR /home/usr


RUN mkdir data
RUN  apt-get -y install libboost-dev
RUN apt-get -y install libboost-program-options-dev
RUN wget ftp://pbil.univ-lyon1.fr/pub/logiciel/silix/silix-1.2.11.tar.gz
RUN tar zxvf silix-1.2.11.tar.gz
WORKDIR /home/usr/silix-1.2.11
RUN ./configure
RUN make
RUN make check
RUN make install




WORKDIR /home/usr

RUN git clone https://github.com/veg/hyphy.git

WORKDIR /home/usr/hyphy
RUN cmake . ;cmake .; make HYPHYMP; make install


WORKDIR /home/usr
RUN git clone https://github.com/tijeco/FUSTr.git

RUN wget -qO- -O tmp.zip https://sourceforge.net/projects/evolveagene/files/EvolvAGene4%20Package.zip && \
      unzip tmp.zip && rm tmp.zip
RUN echo 'for name in *\ *; do mv -v "$name" "${name// /}"; done' > tmp.sh
RUN ["bash", "tmp.sh"]
RUN cp EvolvAGene4Package/EvolveAGene4-linux-x86-64 /usr/bin/EvolveAGene

WORKDIR /home/usr/FUSTr/Simulations/
RUN git pull
RUN ["python3","/home/usr/FUSTr/Simulations/makeSeed.py","8"]

RUN snakemake
# -d /home/usr/FUSTr/Simulations/seeds/
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-4.3.14-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

#RUN conda install numpy biopython scipy -y

RUN pip install numpy biopython scipy
##
ENV PATH /opt/conda/bin:$PATH
# RUN wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_nNklB/gmst_linux_64.tar.gz
RUN ln -sf /bin/bash /bin/sh
ADD $package /home/usr/data
