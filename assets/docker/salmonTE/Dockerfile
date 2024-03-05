# syntax=docker/dockerfile:1

FROM curlimages/curl:7.81.0
WORKDIR /home/curl_user
RUN curl -LO "https://github.com/hyunhwan-jeong/SalmonTE/archive/refs/heads/main.zip" && \
    unzip main.zip && \
    rm -r SalmonTE-main/example SalmonTE-main/test


FROM continuumio/miniconda3:23.3.1-0-alpine
RUN conda update -y -n base -c defaults conda; \
    conda config --add channels conda-forge; \
    conda config --add channels bioconda; \
    conda config --add channels anaconda
RUN conda install -n base -y conda-libmamba-solver ; \
    conda config --set solver libmamba
RUN conda install -y  \
    snakemake=7.28.3 \
    pandas \
    docopt \
    r-base=4.2.3 \
    r-scales \
    r-writexls \
    r-tidyverse \
    bioconductor-deseq2
    
WORKDIR /opt
COPY --from=0 /home/curl_user/SalmonTE-main ./SalmonTE
