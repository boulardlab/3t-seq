# syntax=docker/dockerfile:1

FROM curlimages/curl:7.81.0
WORKDIR /home/curl_user
RUN curl -LO https://github.com/hyunhwan-jeong/SalmonTE/archive/refs/tags/v0.4.tar.gz && \
    tar xf v0.4.tar.gz


FROM continuumio/miniconda3:23.3.1-0
RUN conda update -y -n base -c defaults conda; \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda config --add channels anaconda
RUN conda install -n base -y conda-libmamba-solver ; \
    conda config --set solver libmamba
RUN conda install -y  \
    snakemake \
    pandas \
    docopt \
    r-scales \
    r-writexls \
    r-tidyverse \
    bioconductor-deseq2
    # r-base \
    # r-latticeextra \
    # r-hmisc \
    # r-rms \
    
WORKDIR /opt
COPY --from=0 /home/curl_user/SalmonTE-0.4 ./SalmonTE
