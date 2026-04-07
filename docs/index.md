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
    STAR --> FILTER[BAM Filtering]
    FILTER --> MARKDUP[Picard MarkDuplicates]
    MARKDUP --> BW[BigWig Coverage]

    %% Specialized Analysis
    TRIM --> SALMON_TE[SalmonTE]
    TE_REF --> SALMON_TE

    MARKDUP --> STAR_TE[STAR-TE Analysis]
    TE_REF --> STAR_TE

    MARKDUP --> TRNA[tRNA Quantification]
    TRNA_REF --> TRNA

    %% Quantification
    STAR_TE --> DESEQ2[DESeq2 Expression]
    TRNA --> DESEQ2

    %% Reporting
    QC_RAW --> MULTIQC[MultiQC]
    QC_TRIM --> MULTIQC
    STAR --> MULTIQC
    STAR_TE --> MULTIQC
    MULTIQC --> REPORT((Final Report))

    %% Styling
    style RAW fill:#f9f,stroke:#333,stroke-width:2px
    style REF fill:#bbf,stroke:#333,stroke-width:2px
    style REPORT fill:#bfb,stroke:#333,stroke-width:2px
    style TRIM_DERIVE stroke-dasharray: 5 5
```

## Navigation

- [Getting Started](getting-started.md): Installation and your first run.
- [Configuration](configuration.md): Detailed guide on how to set up your analysis.
- [Usage](usage.md): Commands and workflows.
- [Deduplication](deduplication.md): Understanding the deduplication strategy.
