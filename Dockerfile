# Dockerfile based on GTEx's rnaseq Dockerfile:
# https://github.com/broadinstitute/gtex-pipeline/blob/b53b734a9b096caed237952b34cbce88b38485bd/rnaseq/Dockerfile
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
        build-essential \
        software-properties-common \
        cmake \
        curl \
        git \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        libhdf5-serial-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && apt-get clean \
    && apt-get autoremove -y \
    && rm -rf /var/lib/{apt,dpkg,cache,log}/

# htslib 1.10.2
RUN cd /tmp \
    && wget --no-check-certificate -q https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2 \
    && tar -xf htslib-1.10.2.tar.bz2 && rm htslib-1.10.2.tar.bz2 && cd htslib-1.10.2 \
    && ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs \
    && make -j && make install && ldconfig \
    && cd /tmp && rm -rf htslib-1.10.2

# samtools 1.10
RUN cd /tmp \
    && wget  --no-check-certificate -q https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
    && tar -xf samtools-1.10.tar.bz2 && rm samtools-1.10.tar.bz2 && cd samtools-1.10 \
    && ./configure --with-htslib=system \
    && make -j && make install \
    && cd /tmp && rm -rf samtools-1.10

# STAR 2.6.1d
RUN cd /tmp && \
    wget --no-check-certificate -q https://github.com/alexdobin/STAR/archive/2.6.1d.tar.gz \
    && tar -xf 2.6.1d.tar.gz && rm 2.6.1d.tar.gz \
    && install STAR-2.6.1d/bin/Linux_x86_64/STAR /usr/local/bin/ \
    && install STAR-2.6.1d/bin/Linux_x86_64/STARlong /usr/local/bin/ \
    && rm -rf STAR-2.6.1d

# RSEM 1.3.1
RUN cd /tmp && \
    wget --no-check-certificate -q https://github.com/deweylab/RSEM/archive/v1.3.1.tar.gz \
    && tar -xvf v1.3.1.tar.gz && rm v1.3.1.tar.gz && cd RSEM-1.3.1 \
    && make -j && make install \
    && cd /tmp && rm -rf RSEM-1.3.1

# Picard tools 2.21.4
RUN wget --no-check-certificate -q -P /usr/local/lib/picard-tools/ \
    https://github.com/broadinstitute/picard/releases/download/2.21.4/picard.jar

# RNA-SeQC v2.3.6
RUN cd /tmp \
    && wget --no-check-certificate -q https://github.com/getzlab/rnaseqc/releases/download/v2.3.6/rnaseqc.v2.3.6.linux.gz \
    && gunzip rnaseqc.v2.3.6.linux.gz \
    && install -m755 rnaseqc.v2.3.6.linux /usr/local/bin/rnaseqc \
    && rm rnaseqc.v2.3.6.linux

# Snakemake and other Python packages
RUN python3 -m pip install snakemake==5.20.1 numpy pandas \
    && rm -rf /root/.cache/pip
