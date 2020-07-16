FROM continuumio/miniconda3:4.8.2

# Configure locale and timezone
RUN echo "America/Chicago" > /etc/timezone && \
    rm /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata

RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y python=3.8 \
        snakemake-minimal=5.20.1 \
        pandas=1.0.5 \
        star=2.6.1d \
        samtools=1.10 htslib=1.10 \
    && conda clean -y --all