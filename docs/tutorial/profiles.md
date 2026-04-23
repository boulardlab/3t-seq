# Advanced Profiles

In **3t-seq**, we separate your **biological setup** (what files to process) from your **computational strategy** (how to run them). The primary way to manage this is through the `--profile` flag.

---

## 1. The `--profile` Flag as a Configuration Hub

The `--profile` flag is much more than a simple resource setting. It acts as a **comprehensive encapsulation hub** for all your command-line parameters.

Think of a profile as a pre-packaged "Run Configuration" that points Snakemake to:

- **Resource Limits**: How many CPUs and how much RAM to use.
- **Execution Strategy**: Whether to run locally, via Slurm, or SGE.
- **Environment Management**: Automatic use of Conda or Apptainer (Singularity).
- **Default Configs**: Bundling specific `--config` or `--configfile` flags into the profile folder itself.

---

## 2. Built-in Profiles: Laptop vs. HPC Slurm

We provide two primary profiles as starting points:

### The `laptop` Profile

Designed for small-scale pilot runs or local debugging.

- **Resource Limits**: Pins `--cores` to a safe number (e.g., 4 or 8) so your machine remains responsive.
- **Isolation**: Automatically handles **bind mounts** for containers, allowing them to write results to your local folder.
- **Usage**:

  ```bash
  pixi run snakemake --profile .tests/integration/profiles/laptop --configfile my_config.yaml
  ```

### The `hpc` (Slurm) Profile

Designed for scaling up to hundreds of samples on a High-Performance Computing cluster.

- **Orchestration**: Communicates with the **Slurm scheduler** to submit jobs with appropriate resources.
- **Concurrency**: Allows many jobs to run in parallel (e.g., `--jobs 100`).
- **Resilience**: Configured with automatic re-tries for intermittent cluster failures.
- **Usage**:

  ```bash
  pixi run snakemake --profile profiles/slurm --configfile my_config.yaml
  ```

---

## 3. Profiles as CLI Encapsulation

**Reproducibility** is at the heart of 3t-seq. Typing out long, complex commands with dozens of `--config` flags is error-prone and hard to share.

Instead, we recommend **encapsulating these settings into a profile**. You can create a folder (e.g., `profiles/my_experiment/`) and add a `config.yaml` inside it:

```yaml
# profiles/my_experiment/config.yaml
executor: slurm
jobs: 50
default-resources:
  slurm_account: "my_lab_account"
  runtime: 240
```

Now, instead of a long CLI string, you simply run:

```bash
pixi run snakemake --profile profiles/my_experiment --configfile my_samples.yaml
```

---

## Next Steps

Now that you understand the power of profiles, learn how to execute the pipeline and generate meaningful reports:

- [**Running & Reporting**](running.md): Execution patterns and how to generate the master report.
