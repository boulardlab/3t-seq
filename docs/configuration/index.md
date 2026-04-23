# Configuration Reference

This page is an exhaustive reference for every parameter in `config.yaml`.
The pipeline validates your config against a [JSON Schema](../../workflow/schemas/config.schema.yaml)
before any jobs run — so typos and missing required fields are caught early.

Run the validator manually at any time with:

```bash
pixi run lint local-references laptop
```

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

    Supported labels: `mm10`, `mm39`.

    | Key | Required | Default | Description |
    | :--- | :--- | :--- | :--- |
    | `label` | **Yes** | — | Genome version. |
    | `annotation_type` | No | `ensembl` | GTF format: `ensembl`, `gencode`, or `mgi`. |
    | `selected_chromosomes` | No | all | Restrict analysis to a chromosome list. Useful for pilot runs. |

=== "Custom local files"

    Remove `label` and provide explicit paths. All three keys are required together.

    ```yaml
    genome:
      fasta_path: /data/refs/my_genome.fa
      gtf_path: /data/refs/my_genome.gtf
      annotation_type: ensembl
    ```

    | Key | Required | Default | Description |
    | :--- | :--- | :--- | :--- |
    | `fasta_path` | **Yes** (with gtf) | — | Absolute path to genome FASTA. |
    | `gtf_path` | **Yes** (with fasta) | — | Absolute path to gene annotation GTF. |
    | `annotation_type` | **Yes** (with fasta) | — | `ensembl`, `gencode`, or `mgi`. |
    | `selected_chromosomes` | No | all | Restrict to a chromosome subset. |

    !!! warning
        `label` and `fasta_path`/`gtf_path` are mutually exclusive. The validator will
        reject a config that mixes them.

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
| `protocol` | **Yes** | `pe` (paired-end) or `se` (single-end). Must match the sample sheet columns. |
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
| `trimmomatic` | TruSeq3 defaults | Fixed Trimmomatic string, or `{adaptive: true}` for automatic parameter derivation. |
| `star` | `""` | Extra CLI flags passed to the STAR aligner. |
| `bamCoverage` | `""` | Extra flags for `deeptools bamCoverage`. |

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
