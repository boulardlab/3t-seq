# 3t-seq Manual

Welcome to the **3t-seq** documentation.

3t-seq is a robust, automated Snakemake pipeline designed for end-to-end total RNA-seq analysis, with specialized modules for tRNA and transposable element (TE) processing.

## At a Glance

- **Easy to use**: Simplified configuration and environment management via [Pixi](https://pixi.sh).
- **Reproducible**: Container-native execution with Apptainer/Docker.
- **Comprehensive**: Handles everything from quality control to differential expression.
- **Scalable**: Runs on local machines or high-performance computing (HPC) clusters.

## Workflow Overview

```mermaid
graph TD
    %% Input Layer
    RAW{{Raw FASTQ}}
    REF[(References)]

    %% Reference Processing
    REF --> GENOME[Genome Indexing]
    REF --> TE_REF[TE Libraries]
    REF --> TRNA_REF[tRNA Resources]

    %% Initial QC
    RAW --> QC_RAW[FastQC Raw]

    %% Trimming
    QC_RAW --> TRIM_DERIVE[Adaptive Params]
    RAW --> TRIM[Trimmomatic]
    TRIM_DERIVE -.->|JSON| TRIM

    %% Trimmed QC
    TRIM --> QC_TRIM[FastQC Trimmed]

    %% Alignment Layer
    TRIM --> STAR[STAR Alignment]
    GENOME --> STAR

    %% Post-Alignment
    STAR --> MARKDUP[Picard MarkDuplicates]
    STAR --> BW[BigWig Coverage]

    %% Specialized Analysis
    RAW --> SALMON_TE[SalmonTE quant + test]

    TRIM --> STAR_TE[STAR-TE Alignment]
    GENOME --> STAR_TE
    STAR_TE --> FILTER[TE-only BAM filtering]
    TE_REF --> FILTER
    FILTER --> FC[featureCounts]
    TE_REF --> FC

    STAR --> TRNA[tRNA Quantification]
    TRNA_REF --> TRNA

    %% Differential expression
    STAR --> DESEQ2[DESeq2 Expression]
    FC --> DESEQ2
    TRNA --> DESEQ2

    %% Reporting
    QC_RAW --> MULTIQC[MultiQC]
    QC_TRIM --> MULTIQC
    TRIM --> MULTIQC
    STAR --> MULTIQC
    MARKDUP --> MULTIQC
    MULTIQC --> REPORT((MultiQC report per comparison))

    %% Styling
    style RAW fill:#f9f,stroke:#333,stroke-width:2px
    style REF fill:#bbf,stroke:#333,stroke-width:2px
    style REPORT fill:#bfb,stroke:#333,stroke-width:2px
    style TRIM_DERIVE stroke-dasharray: 5 5
```

`TRIM --> MULTIQC` stands for the Trimmomatic statistics; `STAR --> MULTIQC` for the STAR
logs, gene counts and FastQC of the STAR BAM; `MARKDUP --> MULTIQC` for the duplication
statistics and FastQC of the duplicate-marked BAM. SalmonTE uses the raw reads and the TE
index shipped in its container, and runs its own differential expression test. The
Snakemake HTML report is a separate step (see
[Running & Reporting](tutorial/running.md#2-detailed-reporting)).

## Navigation

- [Getting Started](getting-started/index.md): Installation and your first run.
- [Configuration](configuration/index.md): Detailed guide on how to set up your analysis.
- [Usage](usage/index.md): Commands and workflows.
- [Migrating from 1.x](migration/index.md): Config, sample sheet and environment changes in 2.0.
- [Deduplication](algorithms/deduplication.md): Understanding the deduplication strategy.
