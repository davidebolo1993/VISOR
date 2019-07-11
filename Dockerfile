#Built with:
#sudo docker build -t davidebolo1993/visor .

FROM ubuntu:18.04

# File author/maintainer info
MAINTAINER Davide Bolognini <davidebolognini7@gmail.com>

# Install dependencies
RUN apt-get update && apt-get install -y curl git
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
RUN bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda create -y -n visorenv python=3.6
RUN echo "source activate visorenv" > ~/.bashrc
ENV PATH /miniconda/envs/visorenv/bin:$PATH
RUN conda install -n visorenv -c bioconda samtools wgsim pbsim bwa minimap2 pybedtools pysam pyfaidx numpy
RUN git clone https://github.com/davidebolo1993/VISOR.git && cd VISOR && python setup.py install

#Pull with:
#sudo docker pull davidebolo1993/visor

#Then run:
#sudo docker run davidebolo1993/visor VISOR --help

#Or load the environment
#sudo docker run -ti davidebolo1993/visor
#$(visorenv) VISOR --help
