# Getting Started

## Prerequisites

- **Pixi**: We use [Pixi](https://pixi.sh) for environment and dependency management. Think of Pixi as a toolkit manager. Instead of you manually installing 50 different bioinformatics programs (and hoping they don't conflict with each other or your computer), Pixi creates an isolated, perfectly organized "toolbox" just for this pipeline.

## Installation

1. Get the source code:
    - **Download**: [Latest release](https://github.com/boulardlab/3t-seq/releases) (extract the archive).
    - **Clone**: `git clone https://github.com/boulardlab/3t-seq.git`

2. Install dependencies:

   ```bash
   cd 3t-seq
   pixi install
   ```

## Your First Run

To test the installation, you can run the integration test:

```bash
pixi run test remote-references laptop
```

This will run the pipeline on a small subset of data included in the repository.

!!! tip "Looking for a deeper guide?"
    Check out our new [Step-by-Step Tutorial](tutorial/quickstart.md) for a comprehensive guide on creating custom configurations, managing sample sheets, and scaling execution to HPC clusters via Slurm.

---

## Apple Silicon & macOS Support

3t-seq includes built-in support for macOS, particularly for Apple Silicon (`arm64`, like M1/M2/M3 chips). However, setting up the computational environment on a Mac can be tricky. Here are some things to keep in mind.

!!! warning "GNU Coreutils Shims (Translating Commands)"
    On macOS, standard system commands (like `ln` for making file links, or `sed` for replacing text) are slightly different from the Linux versions that most bioinformatics tools expect.

    When you run the pipeline via Pixi on macOS, a setup script automatically creates "shims" in `.pixi/macos-shims`. You can think of a **shim** as a small translator: when a bioinformatics tool asks for the Linux command, the shim steps in and translates the request so your Mac understands it.

    Ensure you have `coreutils` installed via Homebrew (`brew install coreutils`) so these translators have the right dictionary!

!!! warning "Using Singularity/Apptainer on macOS"
    **Containers** (like Singularity or Apptainer) are virtual software boxes. They pack up an entire program and all its requirements into one file, so it runs identically on your laptop or a supercomputer.

    However, Singularity does not run natively on macOS. If you plan to use containers on Apple Silicon, you must run them through a Linux virtual machine (like [Lima](https://github.com/lima-vm/lima) or [Colima](https://github.com/abiosoft/colima)). Think of this as running a computer-inside-your-computer.

    Most bioinformatics containers are built for Intel/AMD chips (`linux/amd64`). Running these on an Apple Silicon chip requires **architecture emulation** (on-the-fly translation from Intel instructions to Apple instructions). This can lead to:

    - **Performance overhead**: Emulation is significantly slower than running software built natively for your Mac.
    - **Stability issues**: Some older tools (especially Java-based ones) may crash during this translation process.
