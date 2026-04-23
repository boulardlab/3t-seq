# Getting Started

## Prerequisites

### Pixi

3t-seq uses [Pixi](https://pixi.sh) to install and manage all software dependencies.
Without it you would need to manually install and version-lock Snakemake, conda/mamba,
STAR, Trimmomatic, featureCounts, DESeq2, samtools, deeptools, and a dozen other tools —
and ensure they all agree on Python and R versions. Pixi replaces that entire process with a
single command.

Install it once on your machine (Linux / macOS):

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Or follow the [official instructions](https://pixi.sh/latest/#installation).

### Git LFS

3t-seq uses [Git Large File Storage (LFS)](https://git-lfs.com) to store the integration
test data (chr19-subset FASTQs and reference files). Without git-lfs those files are
downloaded as small pointer stubs, not real data. The setup step will then fail, and the
test run will silently process empty inputs.

Install it before cloning:

=== "macOS (Homebrew)"
    ```bash
    brew install git-lfs
    ```

=== "Ubuntu / Debian"
    ```bash
    sudo apt install git-lfs
    ```

=== "EMBL cluster"
    ```bash
    module load git-lfs
    ```

Then activate it globally (once per machine):

```bash
git lfs install
```

---

## Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/boulardlab/3t-seq.git
    cd 3t-seq
    ```

    !!! tip "Development branch"
        For the latest (potentially unstable) features, clone `dev` instead:
        ```bash
        git clone -b dev https://github.com/boulardlab/3t-seq.git
        ```

2. Pull the LFS assets and install all dependencies:

    ```bash
    pixi run -e dev setup   # pulls git-lfs test data
    pixi install            # installs all software
    ```

---

## Your First Run

Run the integration test to verify the installation:

```bash
pixi run test local-references laptop
```

This processes a 4-sample chr19 subset of the [GSE130735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130735)
mouse lung dataset (WT vs KO, paired-end) and checks that all pipeline stages complete
successfully. Expected runtime: ~10–20 minutes on a laptop.

!!! tip "Step-by-step walkthrough"
    The [Tutorial](../tutorial/quickstart.md) walks through this exact test run in detail —
    explaining every config option, sample sheet column, and output file.

---

## Apple Silicon & macOS notes

3t-seq supports macOS including Apple Silicon (M1/M2/M3). A few things to be aware of:

- **GNU coreutils shims**: macOS ships BSD versions of `ln`, `sed`, etc. that differ from
  the Linux versions bioinformatics tools expect. When you run the pipeline via Pixi a
  setup script creates shims in `.pixi/macos-shims`. Install `coreutils` via Homebrew so
  the shims have something to point to:

  ```bash
  brew install coreutils
  ```

- **Singularity / Apptainer**: does not run natively on macOS. Container-based execution
  requires a Linux VM (e.g., [Lima](https://github.com/lima-vm/lima) or
  [Colima](https://github.com/abiosoft/colima)). Most containers are built for `linux/amd64`;
  running them on Apple Silicon adds architecture emulation overhead and may cause instability
  with Java-based tools.
