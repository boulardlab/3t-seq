# Tutorial: Running the Test Dataset

This tutorial walks through running the bundled integration test from scratch. By the end
you will have produced real output files and understood every config option and sample-sheet
column — making it straightforward to adapt the pipeline to your own data.

The test dataset is a chr19 subset of [GSE130735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130735),
a paired-end mouse lung RNA-seq experiment comparing wild-type (WT) and knockout (KO) animals.

---

## 1. Set Up the Environment

If you haven't done this yet, follow [Getting Started](../getting-started/index.md) to install
Pixi, git-lfs, and the pipeline dependencies.

```bash
git clone https://github.com/boulardlab/3t-seq.git
cd 3t-seq
pixi run -e dev setup   # pull LFS test data
pixi install            # install all tools
```

!!! tip "Using the `dev` branch"
      When cloning the repository, the `main` branch is checked out by default.
      Many of the feature described here are available only in the `dev` branch.

      To switch to the `dev` branch, use:

    ```bash
    cd 3t-seq
    git checkout origin/dev
    ```

---

## 2. The Test Data

The integration test lives under `.tests/integration/`.

```text
.tests/integration/
├── sample-sheet.csv          ← the 4-sample manifest
├── configs/
│   ├── local-references.yaml ← config used for the test
│   └── ...
├── profiles/
│   └── laptop/               ← resource settings for local execution
└── GSE130735-subset/         ← chr19-subset FASTQ files (via git-lfs)
```

### The sample sheet

```csv title=".tests/integration/sample-sheet.csv"
name,filename_1,filename_2,genotype
SRX5795112_SRR9016958,GSE130735-subset/SRX5795112_SRR9016958_1.fq.gz,GSE130735-subset/SRX5795112_SRR9016958_2.fq.gz,WT
SRX5795113_SRR9016959,GSE130735-subset/SRX5795113_SRR9016959_1.fq.gz,GSE130735-subset/SRX5795113_SRR9016959_2.fq.gz,WT
SRX5795117_SRR9016963,GSE130735-subset/SRX5795117_SRR9016963_1.fq.gz,GSE130735-subset/SRX5795117_SRR9016963_2.fq.gz,KO
SRX5795118_SRR9016964,GSE130735-subset/SRX5795118_SRR9016964_1.fq.gz,GSE130735-subset/SRX5795118_SRR9016964_2.fq.gz,KO
```

Four columns:

- **`name`** — a unique sample identifier (required, fixed column name).
- **`filename_1`** / **`filename_2`** — paths to the paired-end FASTQ files (required for
  paired-end, fixed column names). Paths are relative to the directory from which Snakemake
  is launched.
- **`genotype`** — a metadata column you choose freely. Here it encodes the biological
  condition used for the differential expression contrast (WT vs KO). The column name must
  match `deseq2.variable` in the config (see below).

### The config

```yaml title=".tests/integration/configs/local-references.yaml"
comparisons:
  - name: GSE130735-subset          # (1)
    protocol: pe                    # (2)
    sample_sheet: sample-sheet.csv  # (3)
    trimmomatic: "ILLUMINACLIP:..."  # (4)
    star: "--seedSearchStartLmax 30 ..."
    bamCoverage: "--binSize 10 --normalizeUsing None"
    deseq2:
      test: Wald                    # (5)
      variable: genotype            # (6)
      reference_level: WT           # (7)

disable_TE_analysis: false
disable_salmonTE_analysis: false
disable_tRNA_analysis: false

globals:
  results_folder: results/          # (8)

genome:
  label: mm10                       # (9)
  fasta_path: references/GRCm38.primary_assembly.genome.chr19.fa  # (10)
  gtf_path: references/MGI.chr19.gff3
  annotation_type: mgi
  selected_chromosomes:
    - chr19                         # (11)
```

1. A label for this comparison — used to name output subdirectories.
2. `pe` = paired-end; `se` = single-end. Must match the columns in your sample sheet.
3. Path to the sample sheet, relative to the Snakemake working directory.
4. Fixed Trimmomatic parameters. Omit this key to use the defaults, or set `adaptive: true`
   to let the pipeline derive parameters automatically from FastQC output.
5. DESeq2 test type (`Wald` or `LRT`).
6. Which sample-sheet column to use for the contrast. Must match a column name in the CSV.
7. The level of `variable` used as the baseline (denominator) of the fold-change. Here WT
   is the reference, so positive LFC means higher in KO.
8. Where all output files will be written.
9. Genome label. The pipeline downloads references automatically for `mm10` and `mm39`.
10. For the test we override with local chr19-only files so the run completes quickly.
11. Restrict analysis to chr19 — drastically speeds up the integration test.

---

## 3. Run the Pipeline

From the repository root:

```bash
pixi run test local-references laptop
```

This expands to:

```bash
pixi run snakemake \
    --profile .tests/integration/profiles/laptop \
    --configfile .tests/integration/configs/local-references.yaml \
    --directory .tests/integration
```

- `--profile` sets resource limits and executor (local, 4 cores by default).
- `--configfile` points to the config above.
- `--directory` sets the working directory, so relative paths in the sample sheet resolve
  correctly.

Snakemake will print a job plan, then start executing. Expected runtime on a laptop: 10–20
minutes.

!!! tip "Dry run first"
    See exactly what will be executed without running anything:
    ```bash
    pixi run dry-run local-references laptop
    ```

---

## 4. Explore the Results

When the run finishes, outputs are in `.tests/integration/results/`:

```text
results/
├── qc/
│   └── multiqc/GSE130735-subset/multiqc_report.html   ← QC overview
├── alignments/
│   └── star/GSE130735-subset/                          ← BAM files
├── analysis/
│   ├── tables/GSE130735-subset/                        ← count matrices
│   └── pictures/GSE130735-subset/                      ← DESeq2 plots
└── tRNA_coverage/
    └── GSE130735-subset/                               ← tRNA quantification
```

Open `results/qc/multiqc/GSE130735-subset/multiqc_report.html` in your browser for a
summary of read counts, mapping rates, and duplication metrics.

---

## 5. Generate a Snakemake Report

```bash
pixi run test-report local-references laptop
```

This creates a `report.zip` in `.tests/integration/`. Unzip it and open `report.html` for
an interactive provenance report that includes tool versions, rule parameters, and QC plots.

---

## Next Steps

Now that you've run the test, adapt the pipeline to your own data:

- [**Preparing Data & Samples**](setup.md): how to write your own sample sheet and config.
- [**Advanced Profiles**](profiles.md): resource settings and HPC/Slurm execution.
- [**Running & Reporting**](running.md): execution patterns and monitoring.
