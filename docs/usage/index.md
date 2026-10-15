# Usage

## Pixi tasks

Pixi wraps the most common Snakemake invocations as named tasks so you don't have to
remember the flags. Every task accepts two positional arguments:
`<config-name>` (the config file stem under `.tests/integration/configs/`) and
`<profile-name>` (a profile folder name).

| Task | Command | What it does |
| :--- | :--- | :--- |
| **setup** | `pixi run -e dev setup` | Installs git-lfs and pulls LFS-tracked test data |
| **lint** | `pixi run lint <config> <profile>` | Validates config YAML against the schema |
| **dry-run** | `pixi run dry-run <config> <profile>` | Shows jobs that would run without executing them |
| **test** | `pixi run test <config> <profile>` | Runs the pipeline on the integration test data |
| **test-report** | `pixi run test-report <config> <profile>` | Generates a Snakemake HTML report after a test run |

### Examples

```bash
# Validate your config before committing to a long run
pixi run lint local-references laptop

# See the job plan without running anything
pixi run dry-run local-references laptop

# Run the integration test on a laptop
pixi run test local-references laptop

# Run using remote references (downloads mm10 from Refgenie)
pixi run test remote-references laptop
```

---

## Direct Snakemake invocation

For production runs or custom flags, call Snakemake directly through Pixi:

```bash
pixi run snakemake \
    --profile profile/slurm \
    --configfile config/my_config.yaml
```

Override a single config value without editing the file:

```bash
pixi run snakemake \
    --profile profile/slurm \
    --configfile config/my_config.yaml \
    --config globals:results_folder="results/rerun_v2"
```

---

## Results layout

Everything is written below `globals.results_folder`; the sub-folders are fixed and
cannot be configured. `<comparison>` is the `name` of an entry under `comparisons`,
`<sample>` the `name` column of its sample sheet.

```text
<results_folder>/
├── qc/
│   ├── fastqc-raw/<comparison>/          FastQC on raw reads
│   ├── fastqc-trimmed/<comparison>/      FastQC on trimmed reads
│   ├── fastqc-markdup/<comparison>/      FastQC on deduplicated alignments
│   ├── star/<comparison>/                FastQC on STAR alignments
│   └── multiqc/<comparison>/multiqc_report.html
├── trim/
│   ├── _shared/<hash>/                   trimmed reads, shared between comparisons
│   └── <comparison>/                     links to _shared/, trimming stats
├── alignments/
│   ├── star/
│   │   ├── _shared/<hash>/               single-copy gene alignments
│   │   └── <comparison>/                 BAM, BAI, bigWig, ReadsPerGene, STAR logs
│   ├── star_markdup/
│   │   ├── _shared/<hash>/
│   │   └── <comparison>/                 duplicate-marked BAM and stats
│   └── starTE/
│       ├── _shared/<hash>/{random,multihit}/
│       └── <comparison>/
│           ├── random/ multihit/         TE alignments
│           ├── filter/                   TE-only BAMs
│           ├── featureCount/             TE counts (random, multihit)
│           ├── DESeq2/                   TE differential expression (lfc.csv, dds_random.rds)
│           └── random/datavzrd/          interactive TE results report
├── salmonTE/
│   ├── quant/<comparison>/               SalmonTE quantification
│   └── de_analysis/<comparison>/         SalmonTE differential expression
├── tRNA_coverage/<comparison>/           tRNA coverage, counts, DESeq2 results, datavzrd report
├── analysis/
│   ├── tables/deseq2/<comparison>/       gene DESeq2 results (results.csv, results.shrink.csv)
│   ├── rdata/deseq2/<comparison>/        DESeq2 object (dds.rds)
│   └── datavzrd/<comparison>/            interactive gene results report
├── references/
│   ├── <label>/                          processed FASTA and GTF
│   ├── STAR/                             STAR index
│   ├── rmsk/                             RepeatMasker annotation
│   ├── gtrnadb/                          GtRNAdb tRNA annotation
│   └── refgenie/                         Refgenie genome config
├── log/<YYYY-MM-DD>/<tool>/              job logs, one folder per run date
└── tmp/                                  scratch space for STAR
```

`salmonTE/`, `alignments/starTE/` and `tRNA_coverage/` are only created when the
corresponding analysis is enabled (see `disable_*_analysis` in the
[configuration reference](../configuration/index.md)). The `_shared/` folders are
explained in [Deduplication](../algorithms/deduplication.md).

---

## Generating a report

```bash
pixi run snakemake \
    --report report.html \
    --profile profile/slurm \
    --configfile config/my_config.yaml
```

Use the **same profile and configfile** as the run that produced the results, otherwise
Snakemake cannot match outputs to rules for the provenance data.

The integration test shortcut:

```bash
pixi run test-report local-references laptop
```

---

## Troubleshooting

### Config validation fails

```text
WorkflowError: Config file is not valid
```

Run `pixi run lint <config> <profile>` to get a human-readable list of what's wrong.
Common causes:

- Missing required key (`globals.results_folder`, `genome.label`, or `comparisons`).
- `fasta_path` provided without `gtf_path` (or vice versa).
- `annotation_type` not set when using custom FASTA/GTF.
- A comparison entry missing `name` or `sample_sheet`.
- A 1.x key still present (e.g. `globals.reads_folder`, `globals.qc_folder`). See
  [Migrating from 1.x](../migration/index.md).

### Sample sheet files not found

```text
Validation failed: The following necessary input read files could not be found
```

- Check that `protocol` in your config matches the columns in your CSV (`protocol: pe`
  requires `filename_1`/`filename_2`; without it the comparison is treated as single-end
  and `filename` is required).
- Relative paths in the sample sheet are resolved relative to `--directory` (the working
  directory passed to Snakemake), not relative to the config file or sample sheet location.
- Confirm the FASTQ files actually exist at those paths on disk.

### STAR runs out of memory

STAR needs ~32 GB of RAM to load the mm10 genome index. On a laptop the `laptop` profile
limits memory, which may cause STAR to fail on large inputs.

- For full genomes, use an HPC cluster with the `hpc`/`slurm` profile.
- For local testing, use `selected_chromosomes: [chr19]` in the `genome` section to
  restrict to a single chromosome.

### Pipeline appears stuck

Check the log files — every rule writes to `results/log/<date>/`. On a cluster you can
also inspect the Slurm queue:

```bash
squeue -u $USER
```

Snakemake's `--rerun-incomplete` flag re-runs any jobs that terminated abnormally:

```bash
pixi run snakemake \
    --profile profile/slurm \
    --configfile config/my_config.yaml \
    --rerun-incomplete
```

### Deprecation warning about `sequencing_libraries`

If you see:

```text
CONFIG DEPRECATION: The 'sequencing_libraries' key has been renamed to 'comparisons'.
```

Open your config file and rename the top-level key from `sequencing_libraries` to
`comparisons`. The old key still works but will be removed in a future release. Other
1.x keys need changes too; see [Migrating from 1.x](../migration/index.md).
