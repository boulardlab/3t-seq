FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="dc6ebc11c614c46262153d8bc966e11a8439195eb3a32e54d44f1bc4460ee1e1"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v2.6.0/utils/datavzrd/environment.yaml
#   prefix: /conda-envs/346c7adef34200145d01dda184ac25b8
#   channels:
#     - conda-forge
#     - nodefaults
#   dependencies:
#     - datavzrd =2.23.8
RUN mkdir -p /conda-envs/346c7adef34200145d01dda184ac25b8
ADD https://github.com/snakemake/snakemake-wrappers/raw/v2.6.0/utils/datavzrd/environment.yaml /conda-envs/346c7adef34200145d01dda184ac25b8/environment.yaml

# Conda environment:
#   source: workflow/env/R.yml
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
COPY workflow/env/R.yml /conda-envs/40cdd57d2470dfd817a34e1ec7edeaa4/environment.yaml

# Conda environment:
#   source: workflow/env/alignment.yml
#   prefix: /conda-envs/8e96037ab9b9dd95318e6dde69e1b470
#   channels:
#     - bioconda
#   dependencies:
#     - star=2.7.6a
#     - subread=2.0.1
RUN mkdir -p /conda-envs/8e96037ab9b9dd95318e6dde69e1b470
COPY workflow/env/alignment.yml /conda-envs/8e96037ab9b9dd95318e6dde69e1b470/environment.yaml

# Conda environment:
#   source: workflow/env/bash.yml
#   prefix: /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181
#   channels:
#     - conda-forge
#   dependencies:
#     - bash=5.1.16
RUN mkdir -p /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181
COPY workflow/env/bash.yml /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181/environment.yaml

# Conda environment:
#   source: workflow/env/bedtools.yml
#   prefix: /conda-envs/7548059a7c044c6fa179ed2c582570cb
#   channels:
#     - bioconda
#   dependencies:
#     - bedtools=2.30.0
RUN mkdir -p /conda-envs/7548059a7c044c6fa179ed2c582570cb
COPY workflow/env/bedtools.yml /conda-envs/7548059a7c044c6fa179ed2c582570cb/environment.yaml

# Conda environment:
#   source: workflow/env/deeptools.yml
#   prefix: /conda-envs/2f716f46231f821a5e905c39c6060cff
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - deeptools=3.5.1
RUN mkdir -p /conda-envs/2f716f46231f821a5e905c39c6060cff
COPY workflow/env/deeptools.yml /conda-envs/2f716f46231f821a5e905c39c6060cff/environment.yaml

# Conda environment:
#   source: workflow/env/pandas.yml
#   prefix: /conda-envs/6f13cf6b27cebfa40962324b0eefc1f4
#   channels:
#     - conda-forge
#     - anaconda
#   dependencies:
#     - pandas=2.2.1
#     - openssl=3.3.0
RUN mkdir -p /conda-envs/6f13cf6b27cebfa40962324b0eefc1f4
COPY workflow/env/pandas.yml /conda-envs/6f13cf6b27cebfa40962324b0eefc1f4/environment.yaml

# Conda environment:
#   source: workflow/env/picard.yml
#   prefix: /conda-envs/5802f2d84ae022c00e054e6c16564f06
#   channels:
#     - bioconda
#   dependencies:
#     - picard=2.27.4
RUN mkdir -p /conda-envs/5802f2d84ae022c00e054e6c16564f06
COPY workflow/env/picard.yml /conda-envs/5802f2d84ae022c00e054e6c16564f06/environment.yaml

# Conda environment:
#   source: workflow/env/qc.yml
#   prefix: /conda-envs/c9db25f4fb5a47a64d4be70e476fd582
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.11 # fix multiqc imp import error
#     - fastqc=0.11.9
#     - multiqc=1.14
RUN mkdir -p /conda-envs/c9db25f4fb5a47a64d4be70e476fd582
COPY workflow/env/qc.yml /conda-envs/c9db25f4fb5a47a64d4be70e476fd582/environment.yaml

# Conda environment:
#   source: workflow/env/samtools.yml
#   prefix: /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools=1.16.1
#     - ucsc-gtftogenepred=447
#     - ucsc-genepredtobed=447
RUN mkdir -p /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d
COPY workflow/env/samtools.yml /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d/environment.yaml

# Conda environment:
#   source: workflow/env/trimmomatic.yml
#   prefix: /conda-envs/b93daf96b2454232db6380819bb61725
#   channels:
#     - bioconda
#   dependencies:
#     - trimmomatic=0.39
RUN mkdir -p /conda-envs/b93daf96b2454232db6380819bb61725
COPY workflow/env/trimmomatic.yml /conda-envs/b93daf96b2454232db6380819bb61725/environment.yaml

# Conda environment:
#   source: workflow/env/wget.yml
#   prefix: /conda-envs/2473cc7cad71f946585d67bc0fc276d6
#   channels:
#     - conda-forge
#     - anaconda
#   dependencies:
#     - wget=1.21.4
#     - gzip
#     - openssl=3.3.0
RUN mkdir -p /conda-envs/2473cc7cad71f946585d67bc0fc276d6
COPY workflow/env/wget.yml /conda-envs/2473cc7cad71f946585d67bc0fc276d6/environment.yaml

# Conda environment:
#   source: workflow/env/yte.yml
#   prefix: /conda-envs/3674dd77e0957c12c05158a88113106a
#   channels:
#     - conda-forge
#   dependencies:
#     - yte=1.5.1
RUN mkdir -p /conda-envs/3674dd77e0957c12c05158a88113106a
COPY workflow/env/yte.yml /conda-envs/3674dd77e0957c12c05158a88113106a/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/346c7adef34200145d01dda184ac25b8 --file /conda-envs/346c7adef34200145d01dda184ac25b8/environment.yaml && \
    mamba env create --prefix /conda-envs/40cdd57d2470dfd817a34e1ec7edeaa4 --file /conda-envs/40cdd57d2470dfd817a34e1ec7edeaa4/environment.yaml && \
    mamba env create --prefix /conda-envs/8e96037ab9b9dd95318e6dde69e1b470 --file /conda-envs/8e96037ab9b9dd95318e6dde69e1b470/environment.yaml && \
    mamba env create --prefix /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181 --file /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181/environment.yaml && \
    mamba env create --prefix /conda-envs/7548059a7c044c6fa179ed2c582570cb --file /conda-envs/7548059a7c044c6fa179ed2c582570cb/environment.yaml && \
    mamba env create --prefix /conda-envs/2f716f46231f821a5e905c39c6060cff --file /conda-envs/2f716f46231f821a5e905c39c6060cff/environment.yaml && \
    mamba env create --prefix /conda-envs/6f13cf6b27cebfa40962324b0eefc1f4 --file /conda-envs/6f13cf6b27cebfa40962324b0eefc1f4/environment.yaml && \
    mamba env create --prefix /conda-envs/5802f2d84ae022c00e054e6c16564f06 --file /conda-envs/5802f2d84ae022c00e054e6c16564f06/environment.yaml && \
    mamba env create --prefix /conda-envs/c9db25f4fb5a47a64d4be70e476fd582 --file /conda-envs/c9db25f4fb5a47a64d4be70e476fd582/environment.yaml && \
    mamba env create --prefix /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d --file /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d/environment.yaml && \
    mamba env create --prefix /conda-envs/b93daf96b2454232db6380819bb61725 --file /conda-envs/b93daf96b2454232db6380819bb61725/environment.yaml && \
    mamba env create --prefix /conda-envs/2473cc7cad71f946585d67bc0fc276d6 --file /conda-envs/2473cc7cad71f946585d67bc0fc276d6/environment.yaml && \
    mamba env create --prefix /conda-envs/3674dd77e0957c12c05158a88113106a --file /conda-envs/3674dd77e0957c12c05158a88113106a/environment.yaml && \
    mamba clean --all -y \
    update-ca-certificates
