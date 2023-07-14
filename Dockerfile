FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="307216c2e491acb1e8c991883cf9f63f8886a4ffb4236f00589097f6db6d6d94"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: env/R.yml
#   prefix: /conda-envs/40cdd57d2470dfd817a34e1ec7edeaa4
#   channels:
#     - conda-forge
#     - r
#     - bioconda
#   dependencies:
#     - r-base=4.0.3
#     - r-hexbin
#     - r-data.table
#     - r-pheatmap
#     - r-hmisc
#     - r-latticeextra
#     - r-r.utils
#     - r-scales
#     - r-rcurl
#     - r-ggplot2
#     - r-openxlsx
#     - r-knitr
#     - r-rmarkdown
#     - r-tidyverse
#     - bioconductor-deseq2
#     - bioconductor-rtracklayer
#     - bioconductor-topgo
#     - bioconductor-vsn
#     - bioconductor-apeglm
#     - bioconductor-reactomepa
RUN mkdir -p /conda-envs/40cdd57d2470dfd817a34e1ec7edeaa4
COPY env/R.yml /conda-envs/40cdd57d2470dfd817a34e1ec7edeaa4/environment.yaml

# Conda environment:
#   source: env/alignment.yml
#   prefix: /conda-envs/8e96037ab9b9dd95318e6dde69e1b470
#   channels:
#     - bioconda
#   dependencies:
#     - star=2.7.6a
#     - subread=2.0.1
RUN mkdir -p /conda-envs/8e96037ab9b9dd95318e6dde69e1b470
COPY env/alignment.yml /conda-envs/8e96037ab9b9dd95318e6dde69e1b470/environment.yaml

# Conda environment:
#   source: env/bedtools.yml
#   prefix: /conda-envs/7548059a7c044c6fa179ed2c582570cb
#   channels:
#     - bioconda
#   dependencies:
#     - bedtools=2.30.0
RUN mkdir -p /conda-envs/7548059a7c044c6fa179ed2c582570cb
COPY env/bedtools.yml /conda-envs/7548059a7c044c6fa179ed2c582570cb/environment.yaml

# Conda environment:
#   source: env/deeptools.yml
#   prefix: /conda-envs/2f716f46231f821a5e905c39c6060cff
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - deeptools=3.5.1
RUN mkdir -p /conda-envs/2f716f46231f821a5e905c39c6060cff
COPY env/deeptools.yml /conda-envs/2f716f46231f821a5e905c39c6060cff/environment.yaml

# Conda environment:
#   source: env/picard.yml
#   prefix: /conda-envs/5802f2d84ae022c00e054e6c16564f06
#   channels:
#     - bioconda
#   dependencies:
#     - picard=2.27.4
RUN mkdir -p /conda-envs/5802f2d84ae022c00e054e6c16564f06
COPY env/picard.yml /conda-envs/5802f2d84ae022c00e054e6c16564f06/environment.yaml

# Conda environment:
#   source: env/qc.yml
#   prefix: /conda-envs/86f8c98916543a9e87d4bc2611d26536
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - fastqc=0.11.9
#     - multiqc=1.14
RUN mkdir -p /conda-envs/86f8c98916543a9e87d4bc2611d26536
COPY env/qc.yml /conda-envs/86f8c98916543a9e87d4bc2611d26536/environment.yaml

# Conda environment:
#   source: env/samtools.yml
#   prefix: /conda-envs/0b7012a9a9c4bff84185eb9d96cb0332
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools=1.16.1
#     - ucsc-gtftogenepred
#     - ucsc-genepredtobed
RUN mkdir -p /conda-envs/0b7012a9a9c4bff84185eb9d96cb0332
COPY env/samtools.yml /conda-envs/0b7012a9a9c4bff84185eb9d96cb0332/environment.yaml

# Conda environment:
#   source: env/trimmomatic.yml
#   prefix: /conda-envs/b93daf96b2454232db6380819bb61725
#   channels:
#     - bioconda
#   dependencies:
#     - trimmomatic=0.39
RUN mkdir -p /conda-envs/b93daf96b2454232db6380819bb61725
COPY env/trimmomatic.yml /conda-envs/b93daf96b2454232db6380819bb61725/environment.yaml

# Conda environment:
#   source: env/wget.yml
#   prefix: /conda-envs/8deaf44f9ffd29816443812db4b8bb83
#   channels:
#     - anaconda
#     - conda-forge
#   dependencies:
#     - wget=1.21.3
#     - gzip
RUN mkdir -p /conda-envs/8deaf44f9ffd29816443812db4b8bb83
COPY env/wget.yml /conda-envs/8deaf44f9ffd29816443812db4b8bb83/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/40cdd57d2470dfd817a34e1ec7edeaa4 --file /conda-envs/40cdd57d2470dfd817a34e1ec7edeaa4/environment.yaml && \
    mamba env create --prefix /conda-envs/8e96037ab9b9dd95318e6dde69e1b470 --file /conda-envs/8e96037ab9b9dd95318e6dde69e1b470/environment.yaml && \
    mamba env create --prefix /conda-envs/7548059a7c044c6fa179ed2c582570cb --file /conda-envs/7548059a7c044c6fa179ed2c582570cb/environment.yaml && \
    mamba env create --prefix /conda-envs/2f716f46231f821a5e905c39c6060cff --file /conda-envs/2f716f46231f821a5e905c39c6060cff/environment.yaml && \
    mamba env create --prefix /conda-envs/5802f2d84ae022c00e054e6c16564f06 --file /conda-envs/5802f2d84ae022c00e054e6c16564f06/environment.yaml && \
    mamba env create --prefix /conda-envs/86f8c98916543a9e87d4bc2611d26536 --file /conda-envs/86f8c98916543a9e87d4bc2611d26536/environment.yaml && \
    mamba env create --prefix /conda-envs/0b7012a9a9c4bff84185eb9d96cb0332 --file /conda-envs/0b7012a9a9c4bff84185eb9d96cb0332/environment.yaml && \
    mamba env create --prefix /conda-envs/b93daf96b2454232db6380819bb61725 --file /conda-envs/b93daf96b2454232db6380819bb61725/environment.yaml && \
    mamba env create --prefix /conda-envs/8deaf44f9ffd29816443812db4b8bb83 --file /conda-envs/8deaf44f9ffd29816443812db4b8bb83/environment.yaml && \
    mamba clean --all -y
