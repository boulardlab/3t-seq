# Configuration Reference

This page is an exhaustive reference for every parameter in `config.yaml`.
The pipeline validates your config against a JSON Schema
(`workflow/schemas/config.schema.yaml`) every time the workflow starts, before any job
runs, so typos and missing required fields are caught early. All problems are reported
together in one `Invalid config file:` error. The sample sheets and the FASTQ files they
list are checked right after.

To check a config without running anything, do a dry-run with it:

```bash
pixi run snakemake --profile profile/default --configfile my_config.yaml -n
```

`pixi run lint` is a different check: it runs `snakemake --lint` on the workflow code with
one of the bundled test configs and does not look at your config.

---

## Minimal working config

```yaml
globals:
  results_folder: results/my_analysis/

genome:
  label: mm10

comparisons:
  - name: my_experiment
    protocol: pe
    sample_sheet: samples.csv
```

This processes a paired-end dataset using mm10 references downloaded automatically,
with all defaults applied.

---

## `globals`

| Key | Required | Default | Description |
| :--- | :--- | :--- | :--- |
| `results_folder` | **Yes** | — | Root directory where all outputs (BAMs, counts, logs) are written. Use a unique path per experiment. |

---

## `genome`

Controls which reference genome is used.

=== "Automatic (Refgenie)"

    Provide a `label` and omit `fasta_path` / `gtf_path`. References are fetched automatically.

    ```yaml
    genome:
      label: mm10
    ```

    Supported labels: `mm10`, `mm39`. Any other label without `fasta_path`/`gtf_path` is
    rejected with `genome.label: '<label>' has no automatically fetched references
    (supported: mm10, mm39). For any other genome set 'genome.fasta_path',
    'genome.gtf_path' and 'genome.annotation_type'.`

    | Key | Required | Default | Description |
    | :--- | :--- | :--- | :--- |
    | `label` | **Yes** | — | Genome version. |
    | `annotation_type` | No | `ensembl` | GTF format: `ensembl`, `gencode`, or `mgi`. |
    | `selected_chromosomes` | No | all | Restrict analysis to a chromosome list. Useful for pilot runs. |

=== "Custom local files"

    Add explicit paths next to `label`. `fasta_path`, `gtf_path` and `annotation_type`
    are required together. With custom files `label` can be any name of 1–64 letters,
    digits, `.`, `_` or `-` (starting with a letter or digit); it is used in file names.

    ```yaml
    genome:
      label: mm10
      fasta_path: /data/refs/my_genome.fa
      gtf_path: /data/refs/my_genome.gtf
      annotation_type: ensembl
    ```

    | Key | Required | Default | Description |
    | :--- | :--- | :--- | :--- |
    | `label` | **Yes** | — | Any safe name; `mm10`/`mm39` if the files are one of those assemblies. |
    | `fasta_path` | **Yes** (with gtf) | — | Absolute path to genome FASTA. |
    | `gtf_path` | **Yes** (with fasta) | — | Absolute path to gene annotation GTF. |
    | `annotation_type` | **Yes** (with fasta) | — | `ensembl`, `gencode`, or `mgi`. |
    | `selected_chromosomes` | No | all | Restrict to a chromosome subset. |

    !!! warning
        The validator rejects a config that sets only some of `fasta_path`, `gtf_path`
        and `annotation_type`.

    RepeatMasker and GtRNAdb annotations are only fetched for `mm10` and `mm39`. If your
    files are one of those assemblies, use that label and every analysis runs. With any
    other label you must set `disable_TE_analysis: true` and `disable_tRNA_analysis: true`,
    otherwise the pipeline stops at start-up:

    ```text
    genome.label 'hg38' is not one of mm10, mm39, so these enabled analyses cannot run:
      - STAR-TE analysis needs a RepeatMasker annotation, which is only fetched for mm10 or mm39. Set 'disable_TE_analysis: true'.
      - tRNA analysis needs GtRNAdb annotations, which are only fetched for mm10 or mm39. Set 'disable_tRNA_analysis: true'.
    If your custom references are an mm10/mm39 assembly, use that label instead.
    ```

    SalmonTE is not affected by this check: it depends on the label prefix and turns
    itself off with a warning (see [SalmonTE supported genomes](#salmonte-supported-genomes)).

---

## `comparisons`

A list of independent analyses. Each entry specifies a sample sheet and its own alignment
and differential expression settings. You can have multiple comparisons in one config —
the pipeline runs them in parallel.

```yaml
comparisons:
  - name: experiment_A
    protocol: pe
    sample_sheet: samples_A.csv
    deseq2:
      variable: genotype
      reference_level: WT
  - name: experiment_B
    protocol: se
    sample_sheet: samples_B.csv
    strandedness: 2
    deseq2:
      variable: treatment
      reference_level: control
```

### Required keys per comparison

| Key | Required | Description |
| :--- | :--- | :--- |
| `name` | **Yes** | Label for this comparison. Used to name output subdirectories. |
| `protocol` | **Yes** | `pe` for paired-end data (sample sheet columns `filename_1`/`filename_2`) or `se` for single-end data (`filename` column). No other value is accepted, and there is no default. |
| `sample_sheet` | **Yes** | Path to the CSV file (relative to `--directory`, or absolute). |

### DESeq2 settings

These can be set per-comparison, or as global defaults under `defaults.deseq2`.
A per-comparison value always takes precedence over the global default.

| Key | Default | Description |
| :--- | :--- | :--- |
| `deseq2.test` | `Wald` | Statistical test: `Wald` or `LRT`. |
| `deseq2.variable` | `genotype` | Sample-sheet column used to define contrast groups. Must match an existing column name. |
| `deseq2.reference_level` | *(none)* | The baseline level of `variable` (fold-change denominator). Strongly recommended — without it DESeq2 picks an arbitrary reference. |

### Alignment overrides

| Key | Default | Description |
| :--- | :--- | :--- |
| `strandedness` | `0` | Library strandedness: `0` unstranded, `1` forward, `2` reverse. |
| `trimmomatic` | TruSeq3 defaults | Fixed Trimmomatic string, or `{extra_params: "...", adaptive: true}`; see [Trimming](../modules/preprocessing.md#trimming-trimmomatic). |
| `star` | `--seedSearchStartLmax 30 --outFilterMismatchNoverReadLmax 0.04 --winAnchorMultimapNmax 40` | Extra CLI flags passed to the STAR aligner. |
| `bamCoverage` | `--binSize 50 --normalizeUsing None` | Extra flags for `deeptools bamCoverage`. |

---

## `defaults`

Global fallbacks applied to every comparison when no per-comparison override is set.
Useful to avoid repeating the same setting in each comparison.

```yaml
defaults:
  strandedness: 2
  deseq2:
    test: Wald
    variable: genotype
    reference_level: WT
```

| Key | Default | Description |
| :--- | :--- | :--- |
| `strandedness` | `0` | Applied to all comparisons that don't set their own `strandedness`. |
| `deseq2.test` | `Wald` | Default test type. |
| `deseq2.variable` | `genotype` | Default contrast variable. |
| `deseq2.reference_level` | *(none)* | Default reference level. |

---

## Module flags

Disable optional analysis modules to save computation time.

```yaml
disable_TE_analysis: false        # set true to skip STAR-TE and SalmonTE
disable_salmonTE_analysis: false  # set true to skip SalmonTE only
disable_tRNA_analysis: false      # set true to skip tRNA quantification
```

### SalmonTE supported genomes

SalmonTE ships TE indexes for a few species only. The pipeline picks the index from the
first two letters of `genome.label`:

| `genome.label` prefix | Species | SalmonTE `--reference_genome` |
| :--- | :--- | :--- |
| `mm` (e.g. `mm10`, `mm39`) | Mouse | `mm` |
| `hg` (e.g. `hg38`) | Human | `hs` |
| `dr` | Zebrafish | `dr` |
| `dm` | Fruit fly | `dm` |

With any other label, SalmonTE and its outputs are skipped and the pipeline logs one
warning at start-up (`SalmonTE DISABLED: ...`); the rest of the pipeline runs normally.
Set `disable_salmonTE_analysis: true` to skip SalmonTE without the warning.

---

## Expert parameters

These settings are rarely needed. Refer to the [Analysis Modules](../modules/te-analysis.md)
documentation for context.

??? info "STAR-TE multi-mapping modes"

    ```yaml
    defaults:
      starTE_random:
        outFilterMultimapNmax: 5000    # default
        winAnchorMultimapNmax: 5000    # default
        alignTranscriptsPerWindowNmax: 300  # default
      starTE_multihit:
        outFilterMultimapNmax: 1       # default
        winAnchorMultimapNmax: 5000    # default
        alignTranscriptsPerWindowNmax: 3000 # default
    ```

    Per-comparison overrides use the same keys under the comparison entry.

??? info "tRNA quantification"

    ```yaml
    tRNA_quantification:
      method: standard          # or "mim-tRNA-seq"
      mimseq_params:
        max_mismatches: 2       # default
        cluster_identity: 0.8   # default
        min_cov: 10             # default
    ```

    Use `mim-tRNA-seq` only if your libraries were prepared with the mim-tRNA-seq protocol.
