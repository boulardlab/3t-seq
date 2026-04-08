# Running & Reporting

Once your data is prepared and you've selected a profile, you're ready to perform the full analysis for your experiment.

---

## 1. Execution Patterns

A standard Snakemake invocation for 3t-seq follows a predictable pattern. Using the `hpc` profile as an example:

```bash
pixi run snakemake \
    --profile profiles/slurm \
    --configfile my_config.yaml
```

The pipeline will then automatically:

1. **Validate** your configuration and sample sheets.
2. **Download** genomic references (if not already local).
3. **Deploy** Conda environments or Singularity containers.
4. **Execute** the rules in the optimal order for your cluster.

### Overriding Parameters

Sometimes you may want to override a single parameter without changing your `config.yaml`. You can do this using the `--config` flag:

```bash
pixi run snakemake --profile profiles/slurm --configfile my_config.yaml --config globals:results_folder="results/rerun_test"
```

---

## 2. Detailed Reporting

Snakemake can generate a comprehensive, self-contained HTML report that includes quality metrics, statistics, and interactive plots for every stage of your pipeline.

### Generating the Report

To ensure the report accurately reflects the run that generated the results, **use the exact same profile and configuration flags** in your reporting command:

```bash
pixi run snakemake \
    --report report.html \
    --profile profiles/slurm \
    --configfile my_config.yaml
```

### What's in the Report?

The generated `report.html` is a stand-alone file that you can share with colleagues. It includes:
- **Statistics**: Read counts, mapping percentages, and duplication rates.
- **QC Plots**: FastQC, MultiQC, and library-specific metrics.
- **Provenance**: The exact versions of tools, environments, and configuration parameters used for the run.
- **Workflow Steps**: An interactive diagram of the analysis steps.

---

## 3. Monitoring Your Run

You can monitor the progress of your jobs in real-time:
- **Terminal Output**: Snakemake provides a live count of finished/remaining jobs.
- **Log Files**: Every rule execution creates a corresponding log file in your `results/logs/` directory.

---

## Next Steps

For a complete reference of every available parameter, visit the detailed configuration guide:

- [**Configuration Reference**](../configuration/index.md): Exhaustive list of all `config.yaml` parameters and their use-cases.
