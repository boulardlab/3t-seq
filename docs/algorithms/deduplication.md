# Sequence Alignment and Trimming Deduplication

To optimize performance and reduce storage usage, the 3t-seq pipeline implements a deduplication mechanism for expensive processing steps, specifically read trimming and sequence alignment.

## Overview

The pipeline uses **content-based hashing** to uniquely identify every processing step. Think of a "hash" as a unique barcode or fingerprint that represents a specific combination of input data and settings. This ensures that expensive computational work (like mapping millions of reads to a genome) is only performed once. If two samples have the exact same barcode, the software knows they are identical and skips re-doing the work.

## The Hashing Mechanism

The pipeline computes a hierarchy of these **SHA256 barcodes/hashes** for every sample in the `SAMPLE_HASHES` registry before the workflow starts.

### 1. Trimming Hash

This is the root of the dependency tree. It is derived from:

- **Input Files**: Absolute paths to the raw FASTQ files.
- **Processing Parameters**: The specific **Trimmomatic** settings (e.g., adaptive vs. fixed string).
- **Protocol**: Whether the library is **Single-End** or **Paired-End**.

### 2. Alignment Hash

The alignment hash ensures that any upstream changes propagate correctly. It depends on:

- The **Trimming Hash**.
- The **STAR** parameters (default or overrides).
- The **Genome Configuration** (labels, FASTA, GTF).

### 3. Analysis Hashes

Higher-level analysis hashes (like `starTE` or `MarkDup`) are derived from the upstream **Alignment Hash** and module-specific settings:

- **starTE**: Alignment Hash + Strandedness + starTE-specific parameters.
- **MarkDup**: Alignment Hash + Picard settings.

## Shared Result Cache (`_shared/`)

Results for hashed tasks are stored in a centralized directory structure:

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

Instead of recalculating results, the pipeline checks if the corresponding barcode exists in the `_shared/` folder. If found, it creates a **symbolic link** from your library's folder to the shared file.

Think of a symbolic link as a shortcut on your desktop, or a "see page 42" reference in a lab notebook. It looks and acts like the actual file, but takes up zero extra hard-drive space.

`results/alignments/star/GSE123456/sample1.bam` -> `../_shared/<hash>/sample1.bam`

## Benefits

- **Performance**: Drastically reduces compute time when analyzing overlapping cohorts or re-running samples.
- **Storage**: Prevents redundant storage of large BAM and FASTQ files.
- **Reliability**: Ensures that identical inputs always produce identical, shared outputs, reducing variability.

## Usage Transparency

This mechanism is entirely transparent to the user. No changes to configuration files are required; the pipeline handles routing, hashing, and link creation automatically.
