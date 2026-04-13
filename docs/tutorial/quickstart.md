# Quickstart Tutorial

This tutorial provides a step-by-step guide to running the **3t-seq pipeline** for the first time. We will use a small subset of a mouse lung dataset (**GSE130735**) to demonstrate the "Happy Path" of the pipeline.

---

## 1. Environment Setup

3t-seq uses **Pixi** for automated environment management.

### a. Install Pixi

If you haven't installed Pixi yet, run the following command (for Linux/macOS):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Or follow the instructions on [pixi.sh](https://pixi.sh).

### b. Clone and Initialize

1. Clone the repository:

   ```bash
   git clone https://github.com/boulardlab/3t-seq.git
   cd 3t-seq
   ```

2. Initialize Git LFS and pull test data:

   ```bash
   pixi run -e dev setup
   ```

   !!! tip "EMBL Cluster Users"
       On the EMBL cluster, you can load Git LFS with: `module load git-lfs`.

   !!! important
       Testing the pipeline requires Large File Storage (LFS) assets. The command above ensures that `git-lfs` is installed in your local environment and all binary assets (FASTQ, FASTA, etc.) are downloaded.

3. Install all dependencies:

   ```bash
   pixi install
   ```

---

## 2. A 5-Minute Run (Test Data)

We provide a built-in integration test dataset and configuration so you can verify your installation immediately.

### Run the Pipeline

Execute the following command to run the pipeline on the mouse lung subset. This uses the `laptop` profile, which is optimized for local execution.

```bash
pixi run snakemake \
    --profile .tests/integration/profiles/laptop \
    --configfile .tests/integration/configs/local-references.yaml
```

!!! tip
    This run will use reference files located in the `.tests/integration/references/` directory.

---

## 3. Exploring the Results

Once the run completes, you can find the outputs in the `results/` directory as specified in the configuration.

- **Quality Control**: Open `results/qc/multiqc/multiqc_report.html` for an overview of the run metrics.
- **Quantification Tables**: Gene, TE, and tRNA counts are available in `results/analysis/tables/`.
- **Alignments**: Sorted BAM files can be found in `results/alignments/`.

---

## Next Steps

Now that you've successfully run the pipeline on test data, move to the next sections to learn how to prepare your own data and configure the pipeline for your experiments:

1. [**Preparing Data & Samples**](setup.md): How to organize your FASTQ files and write a Sample Sheet.
2. [**Advanced Profiles**](profiles.md): Using the `--profile` flag to manage resources and configurations.
3. [**Running & Reporting**](running.md): Scaling up to HPC clusters and generating detailed reports.
