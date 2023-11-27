#Build with:

FROM ubuntu:latest
LABEL description="VISOR"
LABEL base_image="ubuntu:latest"
LABEL software="VISOR"
LABEL about.home="https://github.com/davidebolo1993/VISOR"
LABEL about.license="GPLv3"

# Install dependencies
ENV DEBIAN_FRONTEND=noninteractive 

#workdir
WORKDIR /opt

RUN apt-get update

RUN apt-get install -y nano \
    curl \
    git \
    build-essential \
    g++ \
    cmake \
    libz-dev \
    libcurl4-openssl-dev \
    libssl-dev libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    python3-distutils python3-dev python3-pip python3-edlib \
    && apt-get -y clean all \
    && rm -rf /var/cache

RUN ln -s /usr/bin/python3 /usr/bin/python

#get htslib
RUN curl -LO https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
    && tar -vxjf htslib-1.18.tar.bz2 \
    && cd htslib-1.18 \
    && make \
    && make install

#get samtools
RUN curl -LO https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
    && tar -vxjf samtools-1.18.tar.bz2 \
    && cd samtools-1.18 \
    && make \
    && make install

#get minimap2
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -
ENV PATH /opt/minimap2-2.26_x64-linux:$PATH

#get bedtools
RUN curl -LO https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static \
    && mv bedtools.static bedtools \
    && chmod +x bedtools

#get VISOR and the required python dependencies
RUN git clone https://github.com/davidebolo1993/VISOR.git \
    && cd VISOR \
    && pip install -r requirements.txt \
    && pip install --upgrade cython \
    && pip install .
