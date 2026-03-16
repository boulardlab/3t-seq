# Sequence Alignment and Trimming Deduplication

To optimize performance and reduce storage usage, the 3t-seq pipeline implements a deduplication mechanism for expensive processing steps, specifically read trimming and sequence alignment.

## Overview

When the same samples (same raw FASTQ files) or samples with identical processing parameters appear in multiple sequencing libraries (series) within a single run, the pipeline identifies these overlaps and executes the processing jobs only once.

## Mechanism: Content-Based Hashing

The deduplication is powered by a hashing system implemented in `workflow/rules/common.smk`. Before the workflow starts, the `populate_sample_registry` function calculates a unique SHA-256 hash for each task based on:

1.  **Input Files**: Absolute paths to the raw FASTQ files.
2.  **Processing Parameters**: Tool-specific settings (e.g., Trimmomatic flags, STAR alignment parameters).
3.  **Reference Genome**: The genome labels and paths to FASTA/GTF files.

### Hashing Levels

-   **Trimming Hash**: Derived from raw input files and Trimmomatic parameters.
-   **Alignment Hash**: Derived from the Trimming Hash (ensuring upstream changes propagate), STAR parameters, and genome reference.
-   **Downstream Hashes**: Similar logic applies to `starTE` and `Picard MarkDuplicates`.

## Directory Structure

Processed files are stored in a `_shared` directory within their respective result folders:

```text
results/
├── trim/
│   └── _shared/
│       └── <trim_hash>/
│           └── <sample>.fastq.gz
└── alignments/
    └── star/
        └── _shared/
            └── <align_hash>/
                └── <sample>.Aligned.sortedByCoord.out.bam
```

### Preservation of Per-Series Structure

To maintain compatibility with downstream analysis tools (like DESeq2) that expect a per-series directory structure, the pipeline automatically creates symbolic links from the per-series folders to the shared processed files:

`results/alignments/star/GSE123456/sample1.bam` -> `../_shared/<hash>/sample1.bam`

## Benefits

-   **Performance**: Drastically reduces compute time when analyzing overlapping cohorts or re-running samples with identical parameters across different series.
-   **Storage**: Prevents redundant storage of large BAM and FASTQ files.
-   **Reliability**: Ensures that identical inputs always produce identical, shared outputs, reducing variability in the analysis.

## Usage Transparency

This mechanism is entirely transparent to the user. No changes to the configuration files are required. The pipeline handles the routing, hashing, and link creation automatically.
