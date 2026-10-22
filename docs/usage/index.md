# Usage

## Pixi tasks

Pixi wraps the Snakemake invocations used on the bundled integration test as named
tasks, so you don't have to remember the flags. They all run on `.tests/integration/`;
for your own data call Snakemake directly (see below).

`lint`, `dry-run`, `test` and `test-report` take two optional positional arguments:

- `<config>`: a config stem under `.tests/integration/configs/`: `remote-references`,
  `local-references` or `adaptive-trim`.
- `<profile>`: a test profile under `.tests/integration/profiles/`: `laptop`, `hpc` or
  `ci` (see [Test profiles](../tutorial/profiles.md#4-test-profiles)).

Defaults: `remote-references laptop` for `lint` and `dry-run`, `local-references laptop`
for `test` and `test-report`. Extra Snakemake options go after `--`, e.g.
`pixi run -e dev test local-references laptop -- --cores 4`.

| Task | Command | What it does |
| :--- | :--- | :--- |
| **setup** | `pixi run -e dev setup` | Runs `git lfs install && git lfs pull` (git-lfs comes with the `dev` environment) to fetch the test data |
| **lint** | `pixi run lint <config> <profile>` | Runs `snakemake --lint` (code-quality checks on the workflow) with a test config |
| **dry-run** | `pixi run dry-run <config> <profile>` | Shows jobs that would run without executing them |
| **test** | `pixi run -e dev test <config> <profile>` | Runs `setup`, then the pipeline on the integration test data |
| **test-report** | `pixi run -e dev test-report <config> <profile>` | Runs `test` (a no-op if results are up to date), then writes the Snakemake report to `.tests/integration/report.zip` |

### Examples

```bash
# Lint the workflow code (developers)
pixi run lint local-references laptop

# See the job plan without running anything
pixi run dry-run local-references laptop

# Run the integration test on a laptop
pixi run -e dev test local-references laptop

# Run using remote references (downloads mm10 from Refgenie)
pixi run -e dev test remote-references laptop
```

---

## Direct Snakemake invocation

For production runs or custom flags, call Snakemake directly through Pixi. The examples
below use `profile/slurm`, your copy of `profile/slurm-example`; use `profile/default` to
run locally (see [Advanced Profiles](../tutorial/profiles.md)). `config/my_config.yaml`
stands for your own config file; `config/examples/` has templates to start from.

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
│   ├── rdata/deseq2/ann.rds              gene annotation used by DESeq2 (from the GTF)
│   ├── datavzrd/<comparison>/            interactive gene results report
│   └── pictures/ notebooks/              created at start-up, currently left empty
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
    --report report.zip \
    --profile profile/slurm \
    --configfile config/my_config.yaml
```

Unzip `report.zip` and open `report.html`; see
[Running & Reporting](../tutorial/running.md#2-detailed-reporting) for what it contains.
Use the **same profile and configfile** as the run that produced the results, otherwise
Snakemake cannot match outputs to rules for the provenance data.

The integration test shortcut, which writes `.tests/integration/report.zip`:

```bash
pixi run -e dev test-report local-references laptop
```

---

## Troubleshooting

### Config validation fails

```text
WorkflowError:
Invalid config file:
  comparisons[0] (comparison 'lung'): missing 'protocol'. Set 'protocol: pe' (paired-end, sample sheet columns filename_1/filename_2) or 'protocol: se' (single-end, column filename). 1.x inferred it from the sample sheet; 2.0 does not.
  globals.qc_folder: Output folders are now derived from 'globals.results_folder'; remove this key.
See the 2.0 migration guide: docs/migration/ (https://github.com/boulardlab/3t-seq/tree/main/docs/migration)
```

Every schema problem is listed, one per line, with the path of the offending key. Keys
removed since 1.x get a hint on what replaces them; any other unknown key lists the keys
allowed at that place. Re-check after editing with a dry-run
(`pixi run snakemake --profile profile/default --configfile config/my_config.yaml -n`);
`pixi run lint` does not read your config. Common causes:

- Missing required key (`globals.results_folder`, `genome.label`, or `comparisons`).
- A comparison entry missing `name`, `sample_sheet` or `protocol`.
- No `deseq2.reference_level` for a comparison (neither in the comparison nor in
  `defaults.deseq2`).
- `bamCoverage` with `--normalizeUsing RPGC` but no `--effectiveGenomeSize`, on a genome
  other than `mm10`/`mm39`.
- `fasta_path` provided without `gtf_path` (or vice versa), or without `annotation_type`.
- `genome.label` other than `mm10`/`mm39` without `fasta_path` and `gtf_path`.
- A 1.x key still present (e.g. `globals.reads_folder`, `globals.qc_folder`). See
  [Migrating from 1.x](../migration/index.md).

### Sample sheet files not found

```text
WorkflowError:
Sample sheet validation failed:
Sample sheet missing: /path/to/samples.csv
Sample sheet /path/to/lung.csv lacks column(s) filename_1, filename_2 required by 'protocol: pe' of comparison 'lung' in the config file.
Serie 'liver' / Sample 'WT_1': Could not resolve filename ('reads/wt_1')
```

All sample sheets are checked when the workflow starts, and every problem is reported in
one error. Each line is one of:

- `Sample sheet missing: …`: the `sample_sheet` path does not exist.
- `Sample sheet … lacks column(s) …`: the columns do not match `protocol`
  (`pe` needs `filename_1` and `filename_2`, `se` needs `filename`; both need `name`).
- `Serie '…' / Sample '…': Could not resolve …`: the FASTQ file was not found, with or
  without the extensions `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`.

To fix them:

- Relative paths in the sample sheet are resolved relative to `--directory` (the working
  directory passed to Snakemake), not relative to the config file or sample sheet location.
- Confirm the FASTQ files actually exist at those paths on disk.

### STAR runs out of memory

The STAR rules request 32 GB for the genome index plus 4 GB per GB of compressed FASTQ
input (at least 5 GB), so at least 37 GB per alignment job. On a cluster these values are
passed to Slurm. Locally they are only enforced if you pass `--resources mem_mb=...`;
otherwise STAR fails when the machine runs out of memory.

- For full genomes, run on a cluster with a Slurm profile based on `profile/slurm-example`.
- For local testing, use `selected_chromosomes: [chr19]` in the `genome` section to
  restrict to a single chromosome.

### Pipeline appears stuck

Check the log files: every rule writes to `<results_folder>/log/<YYYY-MM-DD>/`, dated by
the day the run started. On a cluster you can also inspect the Slurm queue:

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
CONFIG DEPRECATION: 'sequencing_libraries' has been renamed to 'comparisons'. Please update your config file; 'sequencing_libraries' will be removed in a future release. See the migration guide: docs/migration/ (https://github.com/boulardlab/3t-seq/tree/main/docs/migration)
```

Open your config file and rename the top-level key from `sequencing_libraries` to
`comparisons`. The old key still works but will be removed in a future release. Other
1.x keys need changes too; see [Migrating from 1.x](../migration/index.md).

A config with both keys is rejected:

```text
Config has both 'sequencing_libraries' (1.x name) and 'comparisons'. Remove 'sequencing_libraries'. See the migration guide: docs/migration/ (https://github.com/boulardlab/3t-seq/tree/main/docs/migration)
```
