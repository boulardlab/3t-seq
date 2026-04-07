# Configuration Reference

This page describes all the parameters available in the `config.yaml` file. The pipeline uses [Snakemake's schema validation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#validation) to ensure your configuration is correct before starting.

## Example Configuration

Here is a complete `config.yaml` demonstrating most pipeline features:

```yaml
# Global settings for the pipeline
globals:
  results_folder: "results/experiment_1/"

# Reference genome configuration
genome:
  label: "mm10" # Automated download of mm10
  # Uncomment below to use custom local references:
  # fasta_path: "references/genome.fa"
  # gtf_path: "references/annotation.gtf"
  # annotation_type: "ensembl"

# Shared defaults for all libraries
defaults:
  strandedness: 0 # 0: unstranded, 1: stranded, 2: reverse
  deseq2:
    variable: "treatment" # The column in your sample sheet to test
    test: "Wald"

# List of sequencing libraries to process
sequencing_libraries:
  - name: "GSE123456" # Unique identifier (e.g., GEO accession)
    protocol: "pe" # Paired-end
    sample_sheet: "data/gse123456_samples.csv"
    trimmomatic:
      adaptive: true # Enable automated parameter derivation

  - name: "GSE987654" # Another experiment in the same run
    protocol: "se" # Single-end
    sample_sheet: "data/gse987654_samples.csv"
    # Override defaults for this specific library:
    trimmomatic: "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 MINLEN:50"
    star: "--seedSearchStartLmax 30"

# Feature flags (optional)
disable_TE_analysis: false
disable_tRNA_analysis: false
```

## Sample Sheet Format

Each library requires a corresponding CSV sample sheet. This file maps your raw sequencing files to sample names and biological variables.

### Example: `data/gse123456_samples.csv`

```csv
name,filename_1,filename_2,condition
sample1,raw/s1_R1.fq.gz,raw/s1_R2.fq.gz,control
sample2,raw/s2_R1.fq.gz,raw/s2_R2.fq.gz,control
sample3,raw/s3_R1.fq.gz,raw/s3_R2.fq.gz,treated
sample4,raw/s4_R1.fq.gz,raw/s4_R2.fq.gz,treated
```

!!! tip "Column Names"
    - `name`: Must be a unique identifier for the sample.
    - `filename_1` / `filename_2`: Paths to the raw sequence files (FASTQ). These paths must be **relative to the project root**. This means if your config file is at `3t-seq/config.yaml`, and your data is at `3t-seq/data/reads.fq.gz`, your path here should just be `data/reads.fq.gz`.
    - Additional columns (like `condition` or `genotype`) can be used for downstream **DESeq2** analysis to compare groups.

---

## Globals

The `globals` section defines the root paths for your analysis.

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `results_folder` | `string` | **Required**. The main directory where all outputs (BAMs, counts, logs) will be saved. |

---

## Genome

The `genome` section configures your reference resources.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `label` | `string` | - | **Required**. Use `mm10` or `mm39` for automated downloads. |
| `fasta_path` | `string` | - | Local path to a custom genome FASTA file. |
| `gtf_path` | `string` | - | Local path to a custom genome GTF annotation. |
| `annotation_type` | `string` | `ensembl` | Format of the annotation (`mgi`, `gencode`, `ensembl`). |
| `selected_chromosomes` | `array` | `null` | Optional list of chromosomes to process (e.g., `["chr1", "chr2"]`). |

---

## Defaults

The `defaults` section provides global settings for all sequencing libraries. These can be overridden at the library level.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `strandedness` | `integer` | `0` | Library preparation geometry. **`0`**: Unstranded (can't tell which DNA strand was transcribed). **`1`**: Forward stranded. **`2`**: Reversely stranded (most modern dUTP-based RNA-seq kits, like Illumina TruSeq Stranded, are `2`). |
| `starTE_random` | `object` | - | Settings for starTE random mode (multi-mapper assignment). |
| `starTE_multihit` | `object` | - | Settings for starTE multihit mode (fractional counting). |
| `deseq2` | `object` | - | Global DESeq2 parameters (test type, design variable). |

---

## Sequencing Libraries

A list of libraries to be processed by the pipeline.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `name` | `string` | - | **Required**. Unique identifier for the library. |
| `sample_sheet` | `string` | - | **Required**. Path to the CSV sample sheet. |
| `protocol` | `string` | `se` | Sequencing protocol (`se` for single-end, `pe` for paired-end). |
| `trimmomatic` | `object` | - | Custom trimming parameters (can be `adaptive: true`). |
| `strandedness` | `integer` | `0` | Library-specific override for strandedness. |

---

## Feature Flags

Use these to disable specific analysis modules.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `disable_TE_analysis` | `boolean` | `false` | Set to `true` to skip all TE-related steps. |
| `disable_salmonTE_analysis` | `boolean` | `false` | Skip SalmonTE specific quantification. |
| `disable_tRNA_analysis` | `boolean` | `false` | Skip tRNA-related steps. |

---

## tRNA Quantification

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `method` | `string` | `standard` | Quantification method (`standard` or `mim-tRNA-seq`). |
| `mimseq_params` | `object` | - | Specific parameters for the mim-tRNA-seq method. |
