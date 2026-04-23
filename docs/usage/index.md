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
    --profile profiles/slurm \
    --configfile config/my_config.yaml
```

Override a single config value without editing the file:

```bash
pixi run snakemake \
    --profile profiles/slurm \
    --configfile config/my_config.yaml \
    --config globals:results_folder="results/rerun_v2"
```

---

## Generating a report

```bash
pixi run snakemake \
    --report report.html \
    --profile profiles/slurm \
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
- A comparison entry missing `name`, `protocol`, or `sample_sheet`.

### Sample sheet files not found

```text
Validation failed: The following necessary input read files could not be found
```

- Check that `protocol` in your config matches the columns in your CSV (`pe` requires
  `filename_1`/`filename_2`; `se` requires `filename`).
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
    --profile profiles/slurm \
    --configfile config/my_config.yaml \
    --rerun-incomplete
```

### Deprecation warning about `sequencing_libraries`

If you see:

```text
CONFIG DEPRECATION: The 'sequencing_libraries' key has been renamed to 'comparisons'.
```

Open your config file and rename the top-level key from `sequencing_libraries` to
`comparisons`. The old key still works but will be removed in a future release.
