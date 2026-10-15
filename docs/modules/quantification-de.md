# Quantification & Differential Expression

The Quantification & DE module provides statistical assessment of gene, TE, and tRNA expression differences across experimental conditions.

## Overview

Differential expression analysis identifies transcripts that change significantly between biological groups, providing insights into biological processes and molecular mechanisms.

## Workflow

```mermaid
graph LR
    A[BAM Files] --> B[featureCounts]
    B --> C[DESeq2 Analysis]
    C --> D[Differential Expression Reports]
```

## Gene Quantification (featureCounts)

[featureCounts](https://subread.sourceforge.net/featureCounts.html) is used to quantify the abundance of gene transcripts from mapped BAM files.

- **Accurate**: Handles overlapping features.
- **Efficient**: High performance with large genomes.

## Differential Expression (DESeq2)

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a leading Biocunductor package for differential expression analysis.

### Features

- **Shrinkage Estimation**: Accurately estimates dispersion and fold changes.
- **Statistical Tests**: Robust identification of differentially expressed genes.
- **Normalization**: Corrects for differences in library size and sequencing depth.

### Result Visualization

DESeq2 results are visualized through:

- **Volcano Plots**: Highlight significantly changing genes.
- **MA Plots**: Visualize the relationship between expression level and fold change.
- **PCA Plots**: Assess the similarity between biological replicates.

## Parameters & Defaults

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `deseq2.test` | `Wald` | Statistical test (`Wald` or `LRT`). |
| `deseq2.variable` | `genotype` | The column in the sample sheet to use for analysis. |
| `deseq2.reference_level` | - | The baseline level of the `variable` (e.g., `WT`). |

---

## Results

| Location | Description |
| :--- | :--- |
| `results/analysis/tables/deseq2/<comparison>/` | Gene DESeq2 results (`results.csv`, `results.shrink.csv`). |
| `results/analysis/rdata/deseq2/<comparison>/` | Gene DESeq2 object (`dds.rds`). |
| `results/analysis/datavzrd/<comparison>/` | Interactive gene results report. |
| `results/qc/multiqc/` | Integrated reports from the entire pipeline. |

TE and tRNA results are written next to their quantification; see
[Results layout](../usage/index.md#results-layout).
