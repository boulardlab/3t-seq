# Migrating from 1.x to 2.0

3t-seq 2.0 is a breaking release. A 1.x config file (last release: `v1.1.0`) will not
run unchanged. This page lists every change you need to make to an existing project,
in the order you are likely to hit them.

!!! summary "TL;DR"
    1. Install [Pixi](https://pixi.sh) and run everything through `pixi run`.
    2. Rename `sequencing_libraries` → `comparisons` and add `protocol: pe` or
       `protocol: se` to every comparison.
    3. Delete everything under `globals` except `results_folder`.
    4. Move the top-level `deseq2` block to `defaults.deseq2` (or into each comparison).
    5. Replace `genome.fasta_url` / `gtf_url` / `gtrnadb_url` with `genome.label: mm10` or
       `mm39` (references fetched automatically), or with local `fasta_path` / `gtf_path` /
       `annotation_type`. Any genome other than mm10/mm39 also needs
       `disable_TE_analysis: true` and `disable_tRNA_analysis: true`.
    6. Turn the sample sheet `filename*` columns into real paths to your FASTQ files.
    7. Write results to a fresh `results_folder`; 1.x outputs are not reused.

---

## 1. Software environment

| | 1.x | 2.0 |
| :--- | :--- | :--- |
| Install | hand-made conda env with Snakemake `<8` and Singularity | `pixi install` from the repository root |
| Snakemake | 7.x | 8.9.0 (pinned in `pixi.toml`) |
| Run | `snakemake --profile profile/default` | `pixi run snakemake --profile <profile> --configfile <config.yaml>` |
| Slurm | built-in `slurm: True` | `snakemake-executor-plugin-slurm` (shipped in the Pixi env) |

```bash
git clone https://github.com/boulardlab/3t-seq.git
cd 3t-seq
pixi install
pixi run snakemake --profile profile/default --configfile config/my_config.yaml
```

The workflow no longer declares a default config file, so always pass `--configfile`.

### Profiles

Snakemake 8 changed several profile keys. Update any 1.x profile you kept:

| 1.x profile key | 2.0 (Snakemake 8) |
| :--- | :--- |
| `slurm: True` | `executor: slurm` |
| `use-conda: True` + `use-singularity: True` | `software-deployment-method: [conda, apptainer]` (the old flags still work but are deprecated) |

`profile/slurm-example/` is a ready-to-copy Snakemake 8 Slurm profile; the 1.x
`profile/slurm-GSE130735/` was removed (rules now size `mem_mb`, `runtime` and `threads`
themselves). See [Advanced Profiles](../tutorial/profiles.md) for complete examples.

---

## 2. Config file: old → new keys

Every key that was removed (all `globals.*` folders, `genome.fasta_url`, `fasta_link`,
`gtf_url`, `gaf_url`, `rmsk_path`, `gtrnadb_url`, `gtrnadb_bed`, and the top-level
`deseq2` block) is a hard error, not a warning. All of them are reported together in one
`Invalid config file:` error at start-up, each with a hint on what replaces it; see
[§5](#5-check-the-migrated-config) for an example.

### Top level

| 1.x | 2.0 | Notes |
| :--- | :--- | :--- |
| `sequencing_libraries` | `comparisons` | The old key is still accepted and renamed at start-up with a `CONFIG DEPRECATION` warning. It will be removed in a future release. A config with both `sequencing_libraries` and `comparisons` is rejected. |
| `deseq2` (top level) | `defaults.deseq2` | A top-level `deseq2` block is **rejected** at start-up. Move it to `defaults.deseq2` or into `comparisons[].deseq2`. |
| `deseq2.working_directory` | — | Removed. |
| — | `defaults` | New. Global fallbacks for `strandedness`, `starTE_random`, `starTE_multihit`, `deseq2`. |
| — | `disable_salmonTE_analysis` | New. Skip SalmonTE while keeping STAR-TE. |
| — | `tRNA_quantification` | New. `method: standard` (default) or `mim-tRNA-seq`. |
| `disable_TE_analysis`, `disable_tRNA_analysis` | unchanged | |

### `globals`

All derived folders were removed. The schema rejects any key other than `results_folder`
(including `data_folder`, `container_folder` and `tmp_folder`).

| 1.x | 2.0 |
| :--- | :--- |
| `globals.results_folder` | `globals.results_folder` (only key left) |
| `globals.reads_folder` | removed; read paths come from the sample sheet (see [§3](#3-sample-sheet)) |
| `globals.qc_folder` | `<results_folder>/qc` |
| `globals.log_folder` | `<results_folder>/log/<timestamp>` |
| `globals.references_folder` | `<results_folder>/references` |
| `globals.tmp_folder` | `<results_folder>/tmp` |
| `globals.analysis_folder` | `<results_folder>/analysis` |

### Per-comparison keys (formerly per-library)

| 1.x | 2.0 | Notes |
| :--- | :--- | :--- |
| `name` | `name` | Required. |
| `sample_sheet` | `sample_sheet` | Required. Relative paths resolve against the Snakemake working directory. |
| — | `protocol` | **Required: `pe` or `se`.** 1.x detected paired-end from a `filename_2` column; 2.0 reads it from the config, and a comparison without `protocol` fails validation. |
| `trimmomatic` (required) | `trimmomatic` (optional) | String, or `{adaptive: true, extra_params: "..."}` for [adaptive trimming](../algorithms/adaptive-trimming.md). Default: `ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36` (`TruSeq3-SE.fa` for single-end). |
| `star` (required) | `star` (optional) | Default: `--seedSearchStartLmax 30 --outFilterMismatchNoverReadLmax 0.04 --winAnchorMultimapNmax 40`. |
| `bamCoverage` (required) | `bamCoverage` (optional) | Default: `--binSize 50 --normalizeUsing None`. |
| `deseq2.{test,variable}` | `deseq2.{test,variable,reference_level}` | Optional per-comparison override of `defaults.deseq2`. |
| — | `strandedness` | New. `0` unstranded (default), `1` stranded, `2` reversely stranded. |
| — | `starTE_random`, `starTE_multihit` | New. Optional STAR-TE overrides. |

!!! warning "Quoted tool strings"
    The 1.x `config/README.md` example wrapped tool strings in an extra pair of quotes
    inside a YAML block scalar (`trimmomatic: >- "ILLUMINACLIP:..."`). Those quotes become
    part of the string. In 2.0 either drop the parameters to use the defaults, or write
    them as a plain YAML string.

### `genome`

| 1.x | 2.0 | Notes |
| :--- | :--- | :--- |
| `label` (free text, e.g. `hg38`) | `label` | Required. Without local files it must be `mm10` or `mm39` (references fetched with Refgenie). With `fasta_path` + `gtf_path` any label of 1–64 letters, digits, `.`, `_`, `-` is accepted. |
| `fasta_url` | `fasta_path` | Local file only. Optional; overrides the Refgenie FASTA. |
| `gtf_url` | `gtf_path` | Local file only. Optional; overrides the Refgenie GTF. |
| `annotation_type` (incl. `other`) | `annotation_type` | One of `ensembl` (default), `gencode`, `mgi`; `other` is no longer accepted, pick the closest format. Required together with `fasta_path` / `gtf_path`. |
| `gtrnadb_url`, `gtrnadb_bed` | — | Removed (error if present). GtRNAdb is fetched automatically for `mm10`/`mm39`. |
| `rmsk_path` | — | Removed (error if present). RepeatMasker annotation is fetched automatically for `mm10`/`mm39`. |
| `fasta_link`, `gaf_url` | — | Removed (error if present). |
| — | `selected_chromosomes` | New. Restrict the analysis to a list of chromosomes (pilot runs). |

`fasta_path`, `gtf_path` and `annotation_type` must be given together, and `label` is
still required when you use them. RepeatMasker and GtRNAdb annotations are only fetched
for mm10 and mm39, so a 1.x project on another genome (e.g. `hg38` with its own FASTA/GTF)
must turn off STAR-TE and tRNA analysis, or the pipeline stops at start-up:

```yaml
genome:
  label: hg38
  fasta_path: /data/refs/hg38.fa.gz
  gtf_path: /data/refs/hg38.gtf.gz
  annotation_type: gencode

disable_TE_analysis: true
disable_tRNA_analysis: true
```

SalmonTE still runs for `hg*` labels (human TE index); see
[SalmonTE supported genomes](../configuration/index.md#salmonte-supported-genomes).

### Before / after

=== "1.x"

    ```yaml
    sequencing_libraries:
      - name: GSE130735
        sample_sheet: sample-sheet.csv
        trimmomatic: "ILLUMINACLIP:TruSeq3-PE.fa:1:0:15:2 SLIDINGWINDOW:20:22 MINLEN:75"
        star: "--seedSearchStartLmax 30 --outFilterMismatchNoverReadLmax 0.04 --winAnchorMultimapNmax 40"
        bamCoverage: "--binSize 50 --normalizeUsing None"

    disable_TE_analysis: false
    disable_tRNA_analysis: false

    globals:
      reads_folder: reads
      results_folder: results/
      qc_folder: results/qc
      log_folder: results/log
      references_folder: results/references
      tmp_folder: /tmp
      analysis_folder: results/analysis

    genome:
      label: mm10
      annotation_type: ensembl
      fasta_url: https://example.org/mm10.fa.gz
      gtf_url: https://example.org/mm10.gtf.gz
      gtrnadb_url: http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz

    deseq2:
      working_directory: ../../..
      test: Wald
      variable: genotype
      reference_level: wt
    ```

=== "2.0"

    ```yaml
    comparisons:
      - name: GSE130735
        protocol: pe
        sample_sheet: sample-sheet.csv
        # trimmomatic / star / bamCoverage are optional; omit them to use the defaults
        trimmomatic: "ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:1:0:15:2 SLIDINGWINDOW:20:22 MINLEN:75"

    disable_TE_analysis: false
    disable_salmonTE_analysis: false
    disable_tRNA_analysis: false

    globals:
      results_folder: results/

    genome:
      label: mm10

    defaults:
      deseq2:
        test: Wald
        variable: genotype
        reference_level: wt
    ```

---

## 3. Sample sheet

The column names are the same as in 1.x; what changed is what the values mean.

| | 1.x | 2.0 |
| :--- | :--- | :--- |
| Sample ID | `name` | `name` |
| Read files | `filename_1` / `filename_2` (paired) or `filename` (single): bare file stems, extensions forbidden | same columns, but each value is a **path** to the FASTQ file |
| Where reads live | `<reads_folder>/<library name>/<stem>.fastq.gz` | absolute path, or relative to the Snakemake working directory (`--directory`, default: current directory) |
| Extension | must be omitted | optional; `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq` are tried automatically |
| Paired vs single | detected from the presence of `filename_2` | set by `protocol: pe` or `protocol: se` in the config (required) |

To convert a 1.x sheet, prefix every stem with the old `<reads_folder>/<library name>/`:

```csv
# 1.x  (reads_folder: reads, library name: GSE130735)
name,filename_1,filename_2,genotype
SRX5795112_SRR9016958,SRX5795112_SRR9016958_1,SRX5795112_SRR9016958_2,WT

# 2.0
name,filename_1,filename_2,genotype
SRX5795112_SRR9016958,reads/GSE130735/SRX5795112_SRR9016958_1,reads/GSE130735/SRX5795112_SRR9016958_2,WT
```

In 2.0 every sample sheet and every read file it lists is checked before any job starts.
Missing sample sheets, missing columns (for the comparison's `protocol`) and read files
that cannot be found are reported together in a single `Sample sheet validation failed:`
error; see [Troubleshooting](../usage/index.md#sample-sheet-files-not-found).

Every other column is free-form metadata. The column named in `deseq2.variable`
(default `genotype`) defines the contrast. See
[Preparing Your Own Data](../tutorial/setup.md#1-the-sample-sheet) for details.

---

## 4. Results and re-runs

- The results tree is now derived from `results_folder` (see [§2](#globals)), and trimming,
  alignment and STAR-TE outputs are stored under content-hashed `_shared/` folders so
  identical samples in different comparisons are processed once
  ([Deduplication](../algorithms/deduplication.md)).
- Snakemake will not pick up 1.x outputs. Point `results_folder` at a new directory
  instead of re-using a 1.x results tree.
- Reference files and STAR indexes are rebuilt under `<results_folder>/references`.

---

## 5. Check the migrated config

Validate the config and look at the job plan before launching a full run:

```bash
pixi run snakemake --configfile config/my_config.yaml --profile profile/default -n
```

Every 1.x key that was removed is reported in one `Invalid config file:` error, with a
hint per key, for example:

```text
Invalid config file:
  deseq2: The top-level 'deseq2' block was moved: set test/variable/reference_level under 'defaults.deseq2' (or per comparison under 'comparisons[].deseq2'). 'working_directory' and 'notebook_path' were removed.
  genome.gtrnadb_url: Custom GtRNAdb annotations are no longer supported; they are downloaded automatically for mm10 or mm39.
  globals.qc_folder: Output folders are now derived from 'globals.results_folder'; remove this key.
See the 2.0 migration guide: docs/migration/ (https://github.com/boulardlab/3t-seq/tree/main/docs/migration)
```

See [Troubleshooting](../usage/index.md#troubleshooting) for other errors.
