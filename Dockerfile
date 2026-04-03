FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="2a4c9ee664ce41a562c0007df6c2c4782f10edc9f14fb7c10de48cf954e15448"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: ../../workflow/env/R.yml
#   prefix: /conda-envs/07ad9ede6e9fcc2a6f0fafe905771466
#   channels:
#     - bioconda
#     - conda-forge
#     - r
#   dependencies:
#     - r-base>=4.4
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
#     - bioconductor-deseq2=1.50.2
#     - bioconductor-rtracklayer
#     - bioconductor-topgo=2.62.0
#     - bioconductor-apeglm
#     - bioconductor-reactomepa=1.54.0
RUN mkdir -p /conda-envs/07ad9ede6e9fcc2a6f0fafe905771466
COPY ../../workflow/env/R.yml /conda-envs/07ad9ede6e9fcc2a6f0fafe905771466/environment.yaml

# Conda environment:
#   source: ../../workflow/env/alignment.yml
#   prefix: /conda-envs/5707810691930ad68b31a3345c6cfa9d
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - star=2.7.11b
#     - subread=2.1.1
RUN mkdir -p /conda-envs/5707810691930ad68b31a3345c6cfa9d
COPY ../../workflow/env/alignment.yml /conda-envs/5707810691930ad68b31a3345c6cfa9d/environment.yaml

# Conda environment:
#   source: ../../workflow/env/bash.yml
#   prefix: /conda-envs/666f7a7f9eab53c3aad55f72fd3e02b2
#   channels:
#     - conda-forge
#   dependencies:
#     - bash>=5.1.16
RUN mkdir -p /conda-envs/666f7a7f9eab53c3aad55f72fd3e02b2
COPY ../../workflow/env/bash.yml /conda-envs/666f7a7f9eab53c3aad55f72fd3e02b2/environment.yaml

# Conda environment:
#   source: ../../workflow/env/bedtools.yml
#   prefix: /conda-envs/aab52617c40794ab4b77450dd57249f0
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bedtools=2.31.1
RUN mkdir -p /conda-envs/aab52617c40794ab4b77450dd57249f0
COPY ../../workflow/env/bedtools.yml /conda-envs/aab52617c40794ab4b77450dd57249f0/environment.yaml

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
#   prefix: /conda-envs/2be69f68af45bfe533dab091f4daca3c
#   channels:
#     - conda-forge
#     - anaconda
#   dependencies:
#     - pandas=2.2.1
#     - openssl=3.3.0
#     - setuptools<81.0.0
#     - python>=3.9,<3.13
RUN mkdir -p /conda-envs/2be69f68af45bfe533dab091f4daca3c
COPY ../../workflow/env/pandas.yml /conda-envs/2be69f68af45bfe533dab091f4daca3c/environment.yaml

# Conda environment:
#   source: ../../workflow/env/picard.yml
#   prefix: /conda-envs/b736d0d3a3e6c301c2f92186255630b3
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - openjdk=17
#     - picard=2.27.4
RUN mkdir -p /conda-envs/b736d0d3a3e6c301c2f92186255630b3
COPY ../../workflow/env/picard.yml /conda-envs/b736d0d3a3e6c301c2f92186255630b3/environment.yaml

# Conda environment:
#   source: ../../workflow/env/qc.yml
#   prefix: /conda-envs/eefea21802b84f8a9b0891967d8a74b8
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - openjdk=17
#     - python=3.11 # fix multiqc imp import error
#     - fastqc=0.11.9
#     - multiqc=1.14
#     - setuptools<81.0.0
RUN mkdir -p /conda-envs/eefea21802b84f8a9b0891967d8a74b8
COPY ../../workflow/env/qc.yml /conda-envs/eefea21802b84f8a9b0891967d8a74b8/environment.yaml

# Conda environment:
#   source: ../../workflow/env/refgenie.yml
#   prefix: /conda-envs/f27a3e4081582b1d67581fa85cc39f2a
#   name: refgenie
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - refgenie=0.13.0
#     - refgenconf>=0.12.1
#     - yacman=0.9.5
#     - python>=3.9,<3.13
#     - setuptools<81
#     - samtools=1.21
RUN mkdir -p /conda-envs/f27a3e4081582b1d67581fa85cc39f2a
COPY ../../workflow/env/refgenie.yml /conda-envs/f27a3e4081582b1d67581fa85cc39f2a/environment.yaml

# Conda environment:
#   source: ../../workflow/env/samtools.yml
#   prefix: /conda-envs/1b261668d4e3aa64c8224616402eca72
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools=1.21
#     - ucsc-gtftogenepred=482
#     - ucsc-genepredtobed=482
RUN mkdir -p /conda-envs/1b261668d4e3aa64c8224616402eca72
COPY ../../workflow/env/samtools.yml /conda-envs/1b261668d4e3aa64c8224616402eca72/environment.yaml

# Conda environment:
#   source: ../../workflow/env/trimmomatic.yml
#   prefix: /conda-envs/50e30cdbbcff8910a88966012cee653e
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - openjdk=17
#     - trimmomatic=0.40
RUN mkdir -p /conda-envs/50e30cdbbcff8910a88966012cee653e
COPY ../../workflow/env/trimmomatic.yml /conda-envs/50e30cdbbcff8910a88966012cee653e/environment.yaml

# Conda environment:
#   source: ../../workflow/env/wget.yml
#   prefix: /conda-envs/e80319233658ca84f4254bb4039938d2
#   channels:
#     - conda-forge
#     - anaconda
#     - bioconda
#   dependencies:
#     - wget=1.21.4
#     - gzip
#     - openssl=3.3.0
RUN mkdir -p /conda-envs/e80319233658ca84f4254bb4039938d2
COPY ../../workflow/env/wget.yml /conda-envs/e80319233658ca84f4254bb4039938d2/environment.yaml

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
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v9.4.1/utils/datavzrd/environment.yaml
#   prefix: /conda-envs/f4325b5fac23e8aea0f27b96bf2a0bb7
#   channels:
#     - conda-forge
#     - nodefaults
#   dependencies:
#     - datavzrd =2.58.8
#     - yte =1.9.4
#     - numpy
#     - pandas
#     - polars
RUN mkdir -p /conda-envs/f4325b5fac23e8aea0f27b96bf2a0bb7
ADD https://github.com/snakemake/snakemake-wrappers/raw/v9.4.1/utils/datavzrd/environment.yaml /conda-envs/f4325b5fac23e8aea0f27b96bf2a0bb7/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/07ad9ede6e9fcc2a6f0fafe905771466 --file /conda-envs/07ad9ede6e9fcc2a6f0fafe905771466/environment.yaml && \
    mamba env create --prefix /conda-envs/5707810691930ad68b31a3345c6cfa9d --file /conda-envs/5707810691930ad68b31a3345c6cfa9d/environment.yaml && \
    mamba env create --prefix /conda-envs/666f7a7f9eab53c3aad55f72fd3e02b2 --file /conda-envs/666f7a7f9eab53c3aad55f72fd3e02b2/environment.yaml && \
    mamba env create --prefix /conda-envs/aab52617c40794ab4b77450dd57249f0 --file /conda-envs/aab52617c40794ab4b77450dd57249f0/environment.yaml && \
    mamba env create --prefix /conda-envs/b38caa9f62475994293262b96f06dcc8 --file /conda-envs/b38caa9f62475994293262b96f06dcc8/environment.yaml && \
    mamba env create --prefix /conda-envs/2be69f68af45bfe533dab091f4daca3c --file /conda-envs/2be69f68af45bfe533dab091f4daca3c/environment.yaml && \
    mamba env create --prefix /conda-envs/b736d0d3a3e6c301c2f92186255630b3 --file /conda-envs/b736d0d3a3e6c301c2f92186255630b3/environment.yaml && \
    mamba env create --prefix /conda-envs/eefea21802b84f8a9b0891967d8a74b8 --file /conda-envs/eefea21802b84f8a9b0891967d8a74b8/environment.yaml && \
    mamba env create --prefix /conda-envs/f27a3e4081582b1d67581fa85cc39f2a --file /conda-envs/f27a3e4081582b1d67581fa85cc39f2a/environment.yaml && \
    mamba env create --prefix /conda-envs/1b261668d4e3aa64c8224616402eca72 --file /conda-envs/1b261668d4e3aa64c8224616402eca72/environment.yaml && \
    mamba env create --prefix /conda-envs/50e30cdbbcff8910a88966012cee653e --file /conda-envs/50e30cdbbcff8910a88966012cee653e/environment.yaml && \
    mamba env create --prefix /conda-envs/e80319233658ca84f4254bb4039938d2 --file /conda-envs/e80319233658ca84f4254bb4039938d2/environment.yaml && \
    mamba env create --prefix /conda-envs/d6fd740cd80b0b1a2fa95ea615907de9 --file /conda-envs/d6fd740cd80b0b1a2fa95ea615907de9/environment.yaml && \
    mamba env create --prefix /conda-envs/f4325b5fac23e8aea0f27b96bf2a0bb7 --file /conda-envs/f4325b5fac23e8aea0f27b96bf2a0bb7/environment.yaml && \
    mamba clean --all -y
