# Quickstart Tutorial

This tutorial provides a step-by-step guide to running the **3t-seq pipeline** for the first time. We will use a small subset of a mouse lung dataset (**GSE130735**) and demonstrate both local and cluster-based execution.

---

## 1. Environment Setup

The pipeline requires **Pixi** for environment management. If you haven't installed it yet, follow the instructions on [pixi.sh](https://pixi.sh).

1. Clone the repository:

   ```bash
   git clone https://github.com/boulardlab/3t-seq.git
   cd 3t-seq
   ```

2. Install dependencies:

   ```bash
   pixi install
   ```

---

## 2. Preparing Your Run

To run the pipeline on your own data, you need two files: a **Configuration File** and a **Sample Sheet**.

### Writing a Config File (`my_config.yaml`)

Create a file named `my_config.yaml` in your project root. This file tells the pipeline where to save results and which genome to use.

```yaml
globals:
  results_folder: "results/my_first_run/"

genome:
  label: "mm10" # Automated download

sequencing_libraries:
  - name: "GSE130735"
    protocol: "pe"
    sample_sheet: "my_samples.csv"
```

### Writing a Sample Sheet (`my_samples.csv`)

The sample sheet maps your FASTQ files to biological conditions.

```csv
name,filename_1,filename_2,condition
WT_Rep1,.tests/integration/GSE130735-subset/SRX5795112_SRR9016958_1.fq.gz,.tests/integration/GSE130735-subset/SRX5795112_SRR9016958_2.fq.gz,WT
WT_Rep2,.tests/integration/GSE130735-subset/SRX5795113_SRR9016959_1.fq.gz,.tests/integration/GSE130735-subset/SRX5795113_SRR9016959_2.fq.gz,WT
KO_Rep1,.tests/integration/GSE130735-subset/SRX5795117_SRR9016963_1.fq.gz,.tests/integration/GSE130735-subset/SRX5795117_SRR9016963_2.fq.gz,KO
KO_Rep2,.tests/integration/GSE130735-subset/SRX5795118_SRR9016964_1.fq.gz,.tests/integration/GSE130735-subset/SRX5795118_SRR9016964_2.fq.gz,KO
```

---

## 3. Understanding Workflow Profiles

In 3t-seq, your **biological configuration** (like sample names and genome choices, which you put in `my_config.yaml`) is completely separate from your **computational configuration** (hardware settings like how many CPUs to use or how to talk to a cluster).

We manage these hardware settings using **Profiles**.

### Local Execution (`laptop`)

The `laptop` profile is optimized for running on your personal computer. It limits the number of CPUs used so your computer doesn't freeze. It also automatically sets up **bind mounts**—think of a bind mount as a secure window that allows the isolated software container to temporarily look into and write files to your project folder.

Usage:

```bash
pixi run snakemake --profile .tests/integration/profiles/laptop --configfile my_config.yaml
```

### Scaling to an HPC Cluster (`slurm`)

For large datasets (anything more than a few pilot samples), your laptop won't be powerful enough. You should run the pipeline on a **High-Performance Computing (HPC) cluster**.
An HPC cluster is essentially a supercomputer constructed by networking many smaller computers (called **nodes**) together.

To manage all the users sharing this supercomputer, clusters use a **scheduler**, most commonly **Slurm**. Imagine Slurm as a restaurant manager:

- You submit a **job** (your order).
- Slurm looks at how many **cores/CPUs** (cooks) and **memory/RAM** (kitchen space) your job needs.
- When enough resources become available, Slurm assigns your job to specific nodes to be executed.

You can create a Slurm profile by creating a folder (e.g., `profiles/slurm/`) and adding a `config.yaml` file. This tells Snakemake how to format its "orders" to the manager:

```yaml
# profiles/slurm/config.yaml
executor: slurm
jobs: 100 # Maximum number of concurrent jobs Snakemake can submit
software-deployment-method: [conda, apptainer]
default-resources:
  slurm_account: my_account # The billing account for your lab
  slurm_partition: interactive # The queue you are submitting to
  runtime: 240 # Maximum minutes a job is allowed to run before being killed
  mem_mb: 16000 # Memory requested in Megabytes (16 GB)
```

Run with:

```bash
pixi run snakemake --profile profiles/slurm --configfile my_config.yaml
```

---

## 4. Execution & Monitoring

A typical Snakemake invocation with custom data looks like this:

```bash
pixi run snakemake \
    --profile .tests/integration/profiles/laptop \
    --configfile my_config.yaml \
    --cores 4 \
    --use-conda \
    --use-singularity
```

### Generating the Report

Once the pipeline completes, you can generate a comprehensive, self-contained HTML report:

```bash
pixi run snakemake --report report.html
```

This report includes quality metrics, statistics, and interactive plots for every step of the pipeline.

---

## 5. Exploring Results

After the run, your `results/my_first_run/` folder will contain:

- **`qc/multiqc/`**: Integrated quality control reports.
- **`alignments/`**: Filtered and sorted BAM files.
- **`analysis/tables/`**: Quantification tables for Genes, TEs, and tRNAs.
- **`analysis/pictures/`**: Differential expression plots (Volcano, PCA).
