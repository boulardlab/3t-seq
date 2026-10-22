# Running & Reporting

Once your data is prepared and you've selected a profile, you're ready to perform the full analysis for your experiment.

---

## 1. Execution Patterns

A standard Snakemake invocation for 3t-seq follows a predictable pattern. Using a Slurm
profile as an example (`profile/slurm` is your copy of `profile/slurm-example`, see
[Advanced Profiles](profiles.md#profileslurm-example-slurm-clusters); use
`profile/default` to run locally):

```bash
pixi run snakemake \
    --profile profile/slurm \
    --configfile my_config.yaml
```

The pipeline will then automatically:

1. **Validate** your configuration and sample sheets.
2. **Download** genomic references (if not already local).
3. **Deploy** the 3t-seq and SalmonTE Apptainer containers (the Conda environments are
   pre-built inside the 3t-seq image).
4. **Execute** the rules in the optimal order for your cluster.

### Overriding Parameters

Sometimes you may want to override a parameter without changing your `config.yaml`. The
most reliable way is a second, small config file: Snakemake merges every `--configfile`
in order, nested keys included, and later files win.

```yaml
# rerun.yaml
globals:
  results_folder: results/rerun_test
```

```bash
pixi run snakemake --profile profile/slurm --configfile my_config.yaml rerun.yaml
```

`--config` only accepts top-level keys (`--config globals:results_folder=...` is rejected
with `Config entry must start with a valid identifier`). A nested value has to be written
as a YAML mapping, e.g. `--config 'globals={results_folder: results/rerun_test}'`.

---

## 2. Detailed Reporting

Snakemake can generate a self-contained report that collects the QC and results files
of the run, together with the workflow graph, run times and the code of every rule.

### Generating the Report

To ensure the report accurately reflects the run that generated the results, **use the exact same profile and configuration flags** in your reporting command:

```bash
pixi run snakemake \
    --report report.zip \
    --profile profile/slurm \
    --configfile my_config.yaml
```

The report is written as a zip archive (the datavzrd reports are folders); unzip it and
open `report.html`. For the integration test, `pixi run -e dev test-report` does the same
and writes `.tests/integration/report.zip`.

### What's in the Report?

The report shows the files the rules mark with `report()`, grouped by category and
comparison:

- **FastQC**: per-sample reports for raw reads, trimmed reads, STAR alignments and
  duplicate-marked alignments.
- **MultiQC**: one `multiqc_report.html` per comparison.
- **Differential expression**: the interactive datavzrd reports for single-copy genes,
  STAR-TE (random mode) and tRNAs.
- **SalmonTE**: the plots of its differential expression test.

Snakemake adds the workflow graph, per-rule run times, and the code and software
environment of each rule.

---

## 3. Monitoring Your Run

You can monitor the progress of your jobs in real-time:
- **Terminal Output**: Snakemake provides a live count of finished/remaining jobs.
- **Log Files**: Every rule execution creates a corresponding log file under
  `<results_folder>/log/<YYYY-MM-DD>/`, dated by the day the run started.

---

## Next Steps

For a complete reference of every available parameter, visit the detailed configuration guide:

- [**Configuration Reference**](../configuration/index.md): Exhaustive list of all `config.yaml` parameters and their use-cases.
