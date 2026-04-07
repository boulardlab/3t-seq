# Alignment

The alignment module maps indexed and trimmed reads to the reference genome, filters high-quality results, and identifies PCR duplicates.

## Workflow

```mermaid
graph LR
    A[Trimmed FASTQ] --> B[STAR Alignment]
    B --> C[BAM Filtering]
    C --> D[Mark Duplicates]
    D --> E[BigWig Tracks]
```

## STAR Alignment

[STAR](https://github.com/alexdobin/STAR) is recognized for its high speed and accuracy in RNA-seq alignment.

### Features

- Support for split-read mapping (splicing).
- Optimized for large genomes.
- Integrated sorting and indexing.

### Configuration

By default, STAR is configured with:

- `outFilterMultimapNmax: 20`
- `outSAMunmapped: Within`
- `quantMode: GeneCounts`

## Post-Alignment Filtering

After alignment, the pipeline filters the resulting BAM files to remove low-quality mappings and ensure data integrity.

### BAM Filtering

- Minimizes noise by discarding reads with low mapping quality scores.
- Standard filter: `q >= 255` (uniquely mapped reads only) or as configured.

### Duplicate Identification (Picard)

[Picard MarkDuplicates](https://broadinstitute.github.io/picard/) is used to identify and flag PCR duplicates, which is crucial for accurate quantification in low-input libraries.

## Parameters & Defaults

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `outFilterMultimapNmax` | `20` | Max alignments for a read to be considered. |
| `outSAMunmapped` | `Within` | Controls if unmapped reads appear in output. |
| `quantMode` | `GeneCounts` | Enables gene-level quantification during alignment. |

---

## Results

| Location | Description |
| :--- | :--- |
| `results/alignments/star/` | BAM files and STAR log files. |
| `results/alignments/star_markdup/` | BAM files with flagged PCR duplicates. |
| `results/qc/star/` | Alignment statistics and log files. |
| `results/alignments/bigwig/` | Genome browser tracks (BigWig). |
