# Sequence Alignment and Trimming Deduplication

To optimize performance and reduce storage usage, the 3t-seq pipeline implements a deduplication mechanism for expensive processing steps: read trimming, STAR alignment, STAR-TE alignment and Picard MarkDuplicates.

## Overview

The pipeline uses **content-based hashing** to uniquely identify every processing step. Think of a "hash" as a unique barcode or fingerprint that represents a specific combination of input data and settings. This ensures that expensive computational work (like mapping millions of reads to a genome) is only performed once. If two samples have the exact same barcode, the software knows they are identical and skips re-doing the work.

## The Hashing Mechanism

The pipeline computes a hierarchy of these hashes for every sample in the `SAMPLE_HASHES`
registry before the workflow starts. Each hash is the SHA256 of the inputs below
serialized as JSON with sorted keys, truncated to the first 16 hex characters
(`calculate_content_hash` in `workflow/rules/common.smk`).

### 1. Trimming Hash

This is the root of the dependency tree. It is derived from:

- **Input Files**: The paths of the raw FASTQ files as resolved from the sample sheet
  (relative paths made absolute against the working directory, extension added if it
  was omitted).
- **Processing Parameters**: The comparison's `trimmomatic` value as written in the
  config: `default` when omitted, the plain string, or the object (`trimmomatic: adaptive`
  counts as `{adaptive: true}`).
- **Protocol**: Whether the library is **Single-End** or **Paired-End**.
- **Adaptive pool** (adaptive trimming only): the sorted raw read paths of every sample
  in the comparison. [Adaptive trimming](adaptive-trimming.md) derives one parameter set
  from all of them, so a sample shared by two adaptive comparisons with different sample
  sets is trimmed once per comparison, each time with that comparison's parameters.
  Identical comparisons still share the result.

### 2. Alignment Hash

The alignment hash ensures that any upstream changes propagate correctly. It depends on:

- The **Trimming Hash**.
- The **STAR** parameters (default or overrides).
- The **Genome Configuration** (labels, FASTA, GTF).

### 3. Analysis Hashes

Two more hashes cover the STAR-TE alignments and duplicate marking:

- **starTE**: the **Trimming Hash**, the genome configuration (label, FASTA, GTF), the
  resolved `strandedness`, `starTE_random` and `starTE_multihit` settings, and a fixed
  version tag. It does not depend on the Alignment Hash, so changing `star` does not
  re-run STAR-TE.
- **MarkDup**: the **Alignment Hash** plus a fixed tag (`picard_v1_fixed`); Picard has no
  configurable settings.

## Shared Result Cache (`_shared/`)

Results for hashed tasks are stored in a centralized directory structure:

```text
<results_folder>/
├── trim/_shared/<trim_hash>/
│   ├── <sample>.fastq.gz                       (se)
│   ├── <sample>_1.fastq.gz, <sample>_2.fastq.gz (pe)
│   ├── <sample>_1.unpaired.fastq.gz, <sample>_2.unpaired.fastq.gz (pe)
│   ├── <sample>.summary.txt, <sample>.stats.txt
│   └── <sample>.adaptive_params.json           (adaptive trimming only; same content for every sample of the comparison)
└── alignments/
    ├── star/_shared/<align_hash>/
    │   └── <sample>.Aligned.sortedByCoord.out.bam, ReadsPerGene.out.tab, ...
    ├── star_markdup/_shared/<markdup_hash>/
    │   └── <sample>.markdup.bam, <sample>.markdup.stats.txt
    └── starTE/_shared/<starte_hash>/{random,multihit}/
        └── <sample>.Aligned.out.bam, <sample>.Log.final.out
```

### Preservation of Per-Series Structure

Each comparison folder holds **symbolic links** into `_shared/`. The `symlink_*` rules
always create them with `ln -sfr`; there is no explicit lookup. What avoids recomputation
is Snakemake itself: the `*_shared` rule that writes `_shared/<hash>/...` only runs when
that file is missing or out of date, so a second comparison with the same hash just gets
new links.

Think of a symbolic link as a shortcut on your desktop, or a "see page 42" reference in a lab notebook. It looks and acts like the actual file, but takes up zero extra hard-drive space.

```text
<results_folder>/alignments/star/<comparison>/<sample>.Aligned.sortedByCoord.out.bam
    -> ../_shared/<align_hash>/<sample>.Aligned.sortedByCoord.out.bam
```

## Benefits

- **Performance**: Drastically reduces compute time when analyzing overlapping cohorts or re-running samples.
- **Storage**: Prevents redundant storage of large BAM and FASTQ files.
- **Reliability**: Ensures that identical inputs always produce identical, shared outputs, reducing variability.

## Usage Transparency

This mechanism is entirely transparent to the user. No changes to configuration files are required; the pipeline handles routing, hashing, and link creation automatically.
