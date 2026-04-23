# Preparing Your Own Data

This page explains how to create the two files you need to run 3t-seq on your own data:
a **sample sheet** (CSV) and a **config file** (YAML). Both are demonstrated using the
bundled GSE130735 test as a concrete reference point.

---

## 1. The Sample Sheet

The sample sheet is a CSV file that tells the pipeline where your FASTQ files are and
what metadata each sample carries.

### Fixed column names

The pipeline looks for these exact column names:

| Column | Required when | Description |
| :--- | :--- | :--- |
| `name` | Always | Unique sample identifier. Used to name output files. |
| `filename` | Single-end only | Path to the FASTQ file. |
| `filename_1` | Paired-end only | Path to the R1 / mate 1 FASTQ file. |
| `filename_2` | Paired-end only | Path to the R2 / mate 2 FASTQ file. |

These four names are reserved — do not use them for metadata.

### Metadata columns

Every other column is a **metadata column** and you choose the names freely.
You can add as many as you like (`genotype`, `treatment`, `sex`, `batch`, `timepoint`, …).

The pipeline carries all metadata columns through to DESeq2, but **only the column named in
`deseq2.variable` in your config is used to build the differential expression contrast**.
All other metadata columns are informational; they appear in DESeq2's `colData` but are
not used in the model unless you add them explicitly.

!!! note "Schema default"
    The default value of `deseq2.variable` in the schema is `genotype`. If you do not set
    `deseq2.variable` explicitly, your sample sheet must have a column named `genotype`.

### Path resolution

- **Absolute paths** are used as-is.
- **Relative paths** are resolved relative to the Snakemake working directory (the value
  of `--directory`, or the current directory if not set).
- The `.fastq.gz` / `.fq.gz` extension is optional — the pipeline tests common extensions
  automatically if the exact path is not found.

### Examples

#### Single-end

```csv
name,filename,genotype
WT_Rep1,reads/wt_rep1.fq.gz,WT
WT_Rep2,reads/wt_rep2.fq.gz,WT
KO_Rep1,reads/ko_rep1.fq.gz,KO
KO_Rep2,reads/ko_rep2.fq.gz,KO
```

#### Paired-end (the GSE130735 test)

```csv
name,filename_1,filename_2,genotype
SRX5795112_SRR9016958,GSE130735-subset/SRX5795112_SRR9016958_1.fq.gz,GSE130735-subset/SRX5795112_SRR9016958_2.fq.gz,WT
SRX5795113_SRR9016959,GSE130735-subset/SRX5795113_SRR9016959_1.fq.gz,GSE130735-subset/SRX5795113_SRR9016959_2.fq.gz,WT
SRX5795117_SRR9016963,GSE130735-subset/SRX5795117_SRR9016963_1.fq.gz,GSE130735-subset/SRX5795117_SRR9016963_2.fq.gz,KO
SRX5795118_SRR9016964,GSE130735-subset/SRX5795118_SRR9016964_1.fq.gz,GSE130735-subset/SRX5795118_SRR9016964_2.fq.gz,KO
```

#### Multiple metadata columns

You can record additional information without affecting the analysis:

```csv
name,filename_1,filename_2,genotype,sex,batch
s01,/data/s01_R1.fq.gz,/data/s01_R2.fq.gz,WT,M,b1
s02,/data/s02_R1.fq.gz,/data/s02_R2.fq.gz,WT,F,b1
s03,/data/s03_R1.fq.gz,/data/s03_R2.fq.gz,KO,M,b2
s04,/data/s04_R1.fq.gz,/data/s04_R2.fq.gz,KO,F,b2
```

With this sample sheet and `deseq2.variable: genotype` in the config, the contrast is
WT vs KO. The `sex` and `batch` columns are available in the DESeq2 object's `colData`
but are not included in the model.

---

## 2. The Config File

The config file (`config.yaml`) controls every aspect of the run. The minimum viable
config has three sections:

```yaml
globals:
  results_folder: results/my_analysis/

genome:
  label: mm10

comparisons:
  - name: my_experiment
    protocol: pe
    sample_sheet: samples.csv
    deseq2:
      variable: genotype
      reference_level: WT
```

### The `comparisons` list

Each entry in `comparisons` represents one independent analysis: one sample sheet, one
set of alignment parameters, and one DESeq2 contrast.

| Key | Required | Default | Description |
| :--- | :--- | :--- | :--- |
| `name` | Yes | — | Label for this comparison. Used to name output subdirectories. |
| `protocol` | Yes | — | `pe` (paired-end) or `se` (single-end). Must match the sample sheet columns. |
| `sample_sheet` | Yes | — | Path to the CSV file. |
| `deseq2.test` | No | `Wald` | Statistical test: `Wald` or `LRT`. |
| `deseq2.variable` | No | `genotype` | Sample-sheet column used to define groups. |
| `deseq2.reference_level` | No | *(none)* | Baseline level of `variable`. Defines the denominator of the fold-change. Strongly recommended. |
| `strandedness` | No | `0` | `0` unstranded, `1` forward-stranded, `2` reverse-stranded. |
| `trimmomatic` | No | TruSeq3 defaults | Fixed string or `{adaptive: true}` for automatic parameter derivation. |
| `star` | No | `""` | Extra flags passed to STAR. |
| `bamCoverage` | No | `""` | Extra flags for deeptools bamCoverage. |

### The `genome` section

| Key | Required | Default | Description |
| :--- | :--- | :--- | :--- |
| `label` | Yes | — | `mm10` or `mm39`. Downloads references automatically via Refgenie. |
| `fasta_path` | No | derived from `label` | Absolute path to a local FASTA. Requires `gtf_path` and `annotation_type`. |
| `gtf_path` | No | derived from `label` | Absolute path to a local GTF. Requires `fasta_path` and `annotation_type`. |
| `annotation_type` | No | `ensembl` | `ensembl`, `gencode`, or `mgi`. Required when using custom FASTA/GTF. |
| `selected_chromosomes` | No | all | Restrict to a subset of chromosomes (useful for pilot runs). |

### Using your own reference files

Remove `label` and provide explicit paths:

```yaml
genome:
  fasta_path: /data/refs/custom.fa
  gtf_path: /data/refs/custom.gtf
  annotation_type: ensembl
```

!!! warning
    `label` and `fasta_path`/`gtf_path` are mutually exclusive. The config validator
    will raise an error if you mix them.

---

## 3. Putting It Together

A complete config for a typical WT vs KO experiment:

```yaml
globals:
  results_folder: results/wt_vs_ko/

genome:
  label: mm10

comparisons:
  - name: lung_WT_vs_KO
    protocol: pe
    sample_sheet: samples/lung.csv
    deseq2:
      variable: genotype
      reference_level: WT
```

Run it:

```bash
pixi run snakemake \
    --profile profiles/laptop \
    --configfile config/my_config.yaml
```

---

## Next Steps

- [**Advanced Profiles**](profiles.md): encapsulate resource settings and executor choice.
- [**Running & Reporting**](running.md): HPC execution and HTML report generation.
- [**Configuration Reference**](../configuration/index.md): exhaustive list of all parameters.
