FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="9cee23af23e90e6e41eedaa595815bd36011ca925cd9323685a63749565c102d"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: ../../workflow/env/R.yml
#   prefix: /conda-envs/77b8061a0f6d791781952ff9f600bdcb
#   channels:
#     - bioconda
#     - conda-forge
#     - r
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
#     - bioconductor-deseq2=1.30.1
#     - bioconductor-rtracklayer
#     - bioconductor-topgo=2.42.0
#     - bioconductor-vsn
#     - bioconductor-apeglm
#     - bioconductor-reactomepa=1.34.0
RUN mkdir -p /conda-envs/77b8061a0f6d791781952ff9f600bdcb
COPY ../../workflow/env/R.yml /conda-envs/77b8061a0f6d791781952ff9f600bdcb/environment.yaml

# Conda environment:
#   source: ../../workflow/env/alignment.yml
#   prefix: /conda-envs/8e96037ab9b9dd95318e6dde69e1b470
#   channels:
#     - bioconda
#   dependencies:
#     - star=2.7.6a
#     - subread=2.0.1
RUN mkdir -p /conda-envs/8e96037ab9b9dd95318e6dde69e1b470
COPY ../../workflow/env/alignment.yml /conda-envs/8e96037ab9b9dd95318e6dde69e1b470/environment.yaml

# Conda environment:
#   source: ../../workflow/env/bash.yml
#   prefix: /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181
#   channels:
#     - conda-forge
#   dependencies:
#     - bash=5.1.16
RUN mkdir -p /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181
COPY ../../workflow/env/bash.yml /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181/environment.yaml

# Conda environment:
#   source: ../../workflow/env/bedtools.yml
#   prefix: /conda-envs/7548059a7c044c6fa179ed2c582570cb
#   channels:
#     - bioconda
#   dependencies:
#     - bedtools=2.30.0
RUN mkdir -p /conda-envs/7548059a7c044c6fa179ed2c582570cb
COPY ../../workflow/env/bedtools.yml /conda-envs/7548059a7c044c6fa179ed2c582570cb/environment.yaml

# Conda environment:
#   source: ../../workflow/env/deeptools.yml
#   prefix: /conda-envs/b38caa9f62475994293262b96f06dcc8
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - deeptools=3.5.1
#     - setuptools<81.0.0
#     - matplotlib<3.9
RUN mkdir -p /conda-envs/b38caa9f62475994293262b96f06dcc8
COPY ../../workflow/env/deeptools.yml /conda-envs/b38caa9f62475994293262b96f06dcc8/environment.yaml

# Conda environment:
#   source: ../../workflow/env/pandas.yml
#   prefix: /conda-envs/22b95354972e48bc1415a03941f4a9da
#   channels:
#     - conda-forge
#     - anaconda
#   dependencies:
#     - pandas=2.2.1
#     - openssl=3.3.0
#     - setuptools<81.0.0
RUN mkdir -p /conda-envs/22b95354972e48bc1415a03941f4a9da
COPY ../../workflow/env/pandas.yml /conda-envs/22b95354972e48bc1415a03941f4a9da/environment.yaml

# Conda environment:
#   source: ../../workflow/env/picard.yml
#   prefix: /conda-envs/5802f2d84ae022c00e054e6c16564f06
#   channels:
#     - bioconda
#   dependencies:
#     - picard=2.27.4
RUN mkdir -p /conda-envs/5802f2d84ae022c00e054e6c16564f06
COPY ../../workflow/env/picard.yml /conda-envs/5802f2d84ae022c00e054e6c16564f06/environment.yaml

# Conda environment:
#   source: ../../workflow/env/qc.yml
#   prefix: /conda-envs/d0b615e1a38d223cd5e2527a5f9259c6
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.11 # fix multiqc imp import error
#     - fastqc=0.11.9
#     - multiqc=1.14
#     - setuptools<81.0.0
RUN mkdir -p /conda-envs/d0b615e1a38d223cd5e2527a5f9259c6
COPY ../../workflow/env/qc.yml /conda-envs/d0b615e1a38d223cd5e2527a5f9259c6/environment.yaml

# Conda environment:
#   source: ../../workflow/env/refgenie.yml
#   prefix: /conda-envs/d8b0a908802597799d00afd556e97962
#   name: refgenie
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - refgenie=0.12.1
#     - setuptools<81
#     - samtools=1.16.1
RUN mkdir -p /conda-envs/d8b0a908802597799d00afd556e97962
COPY ../../workflow/env/refgenie.yml /conda-envs/d8b0a908802597799d00afd556e97962/environment.yaml

# Conda environment:
#   source: ../../workflow/env/samtools.yml
#   prefix: /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools=1.16.1
#     - ucsc-gtftogenepred=447
#     - ucsc-genepredtobed=447
RUN mkdir -p /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d
COPY ../../workflow/env/samtools.yml /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d/environment.yaml

# Conda environment:
#   source: ../../workflow/env/trimmomatic.yml
#   prefix: /conda-envs/b93daf96b2454232db6380819bb61725
#   channels:
#     - bioconda
#   dependencies:
#     - trimmomatic=0.39
RUN mkdir -p /conda-envs/b93daf96b2454232db6380819bb61725
COPY ../../workflow/env/trimmomatic.yml /conda-envs/b93daf96b2454232db6380819bb61725/environment.yaml

# Conda environment:
#   source: ../../workflow/env/wget.yml
#   prefix: /conda-envs/f00fd813319f551472eebd27835bf2b0
#   channels:
#     - conda-forge
#     - anaconda
#     - bioconda
#   dependencies:
#     - wget=1.21.4
#     - gzip
#     - openssl=3.3.0
#     - samtools=1.16.1
RUN mkdir -p /conda-envs/f00fd813319f551472eebd27835bf2b0
COPY ../../workflow/env/wget.yml /conda-envs/f00fd813319f551472eebd27835bf2b0/environment.yaml

# Conda environment:
#   source: ../../workflow/env/yte.yml
#   prefix: /conda-envs/d6fd740cd80b0b1a2fa95ea615907de9
#   channels:
#     - conda-forge
#   dependencies:
#     - yte=1.5.1
#     - setuptools<81.0.0
RUN mkdir -p /conda-envs/d6fd740cd80b0b1a2fa95ea615907de9
COPY ../../workflow/env/yte.yml /conda-envs/d6fd740cd80b0b1a2fa95ea615907de9/environment.yaml

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

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/77b8061a0f6d791781952ff9f600bdcb --file /conda-envs/77b8061a0f6d791781952ff9f600bdcb/environment.yaml && \
    mamba env create --prefix /conda-envs/8e96037ab9b9dd95318e6dde69e1b470 --file /conda-envs/8e96037ab9b9dd95318e6dde69e1b470/environment.yaml && \
    mamba env create --prefix /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181 --file /conda-envs/16e1dc5e3e5976d71e955eaf12ac9181/environment.yaml && \
    mamba env create --prefix /conda-envs/7548059a7c044c6fa179ed2c582570cb --file /conda-envs/7548059a7c044c6fa179ed2c582570cb/environment.yaml && \
    mamba env create --prefix /conda-envs/b38caa9f62475994293262b96f06dcc8 --file /conda-envs/b38caa9f62475994293262b96f06dcc8/environment.yaml && \
    mamba env create --prefix /conda-envs/22b95354972e48bc1415a03941f4a9da --file /conda-envs/22b95354972e48bc1415a03941f4a9da/environment.yaml && \
    mamba env create --prefix /conda-envs/5802f2d84ae022c00e054e6c16564f06 --file /conda-envs/5802f2d84ae022c00e054e6c16564f06/environment.yaml && \
    mamba env create --prefix /conda-envs/d0b615e1a38d223cd5e2527a5f9259c6 --file /conda-envs/d0b615e1a38d223cd5e2527a5f9259c6/environment.yaml && \
    mamba env create --prefix /conda-envs/d8b0a908802597799d00afd556e97962 --file /conda-envs/d8b0a908802597799d00afd556e97962/environment.yaml && \
    mamba env create --prefix /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d --file /conda-envs/9a7384898f8f9aa25cd1d29f531a7f7d/environment.yaml && \
    mamba env create --prefix /conda-envs/b93daf96b2454232db6380819bb61725 --file /conda-envs/b93daf96b2454232db6380819bb61725/environment.yaml && \
    mamba env create --prefix /conda-envs/f00fd813319f551472eebd27835bf2b0 --file /conda-envs/f00fd813319f551472eebd27835bf2b0/environment.yaml && \
    mamba env create --prefix /conda-envs/d6fd740cd80b0b1a2fa95ea615907de9 --file /conda-envs/d6fd740cd80b0b1a2fa95ea615907de9/environment.yaml && \
    mamba env create --prefix /conda-envs/346c7adef34200145d01dda184ac25b8 --file /conda-envs/346c7adef34200145d01dda184ac25b8/environment.yaml && \
    mamba clean --all -y
