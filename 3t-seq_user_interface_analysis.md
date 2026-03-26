# 3t-seq Pipeline User Interface Analysis

## Overview
This document analyzes the user interface of the `3t-seq` pipeline from the perspective of an end-user. In the context of a Snakemake workflow, the "user interface" primarily consists of the configuration files (`config.yaml`), the sample sheet (`.csv`/`.tsv` files), and the command-line execution instructions.

## Critical Points & Potential Pain Points
Based on the current structure (`config/README.md` and `config/config.yaml`), here are the critical touchpoints where a user interacts with the pipeline, along with associated challenges:

### 1. The `globals` Path Definitions
**Current State:** The user must explicitly define seven different paths (`reads_folder`, `results_folder`, `qc_folder`, `log_folder`, `references_folder`, `tmp_folder`, `analysis_folder`). 
**Critical Point:** In practice (as seen in `config.yaml`), most of these folders point to subdirectories within `results_folder`. Forcing the user to manually define each of these paths in the configuration file is highly prone to copy-paste errors, absolute vs relative path confusion, and creates unnecessary friction during setup.

### 2. Embedded Tool Parameters
**Current State:** Users are required to provide command-line string arguments for tools like `trimmomatic`, `star`, and `bamCoverage` directly in the `sequencing_libraries` definition block.
**Critical Point:** This is a major usability hurdle. Expecting a biologist or typical user to write raw tool flags like `"ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/... SLIDINGWINDOW:20:22 MAXINFO:4:20 LEADING:3 TRAILING:3 MINLEN:36"` makes the configuration file verbose, brittle, and intimidating. A single typo in these strings will cause the pipeline to crash deep into execution.

### 3. File Path Resolution and Naming Conventions
**Current State:** Reads must be located in a directory matching the `name` of the library inside the `reads_folder`. Additionally, reads must follow strict naming conventions for paired-end identifiers (e.g. `_1`/`_2`, `_R1`/`_R2` before the extension).
**Critical Point:** This implicit pairing mechanism limits flexibility. If a sequencing core returns files that don't match this exact folder structure or naming convention, the user must use symlinks or aggressively rename their files before running the pipeline. Furthermore, sometimes absolute paths are used in sample sheets, which conflicts with the `reads_folder` + `name` implicit resolution described in `config/README.md`.

### 4. Reference Genomes Configuration
**Current State:** Users must specify URLs or absolute paths to `fasta`, `gtf` and `gtrnadb` sources manually.
**Critical Point:** Fetching the right genome builds is a perennial struggle in bioinformatics. Forcing the user to define these URLs manually assumes they know exactly where to find high-quality annotations that match the expected formatting of the rest of the tools in the suite.

## Proposed Simplifications

### A. Auto-derive Results Subdirectories
Reduce the `globals` section to just a `results_folder` and a `reads_folder` (or remove `reads_folder` entirely if input paths are defined via the sample sheet). The pipeline should automatically create `qc/`, `logs/`, `references/`, `tmp/`, and `analysis/` folders inside the `results_folder` as standard.

*Example target config:*
```yaml
globals:
  results_folder: results/ # Pipeline automatically handles internal structures
```

### B. Sensible Defaults for Tool Parameters
Abstract away explicit command-line strings for `trimmomatic`, `star`, and `bamCoverage`. The pipeline should use pre-tested, sensible defaults internally. If advanced users want to override them, they can use an `extra_args` or `params` block, but these should not be strictly mandatory for the average run.

*Example target config:*
```yaml
sequencing_libraries:
  - name: GSE13073
    sample_sheet: sample-sheet.csv
    # trimmomatic, star, and bamCoverage are handled natively with good defaults
```

### C. Enhanced Sample Sheet Parsing
Instead of relying on rigid folder structures (`reads_folder/LIBRARY_NAME`) and strict suffix matching (`_1.fastq.gz`), rely entirely on the sample sheet to map `sample_id` to its raw read files. The sample sheet should allow absolute paths. If `filename_2` is provided, it's PE; if not, it's SE. This makes it trivial to run the pipeline on data sitting anywhere on a server without renaming or moving the data.

### D. Built-in Genome Profiles (iGenomes / Refgenie)
Provide a standard identifier shortcut for common species (e.g. `genome: mm10` or `genome: hg38`) and use an internal lookup table or tools like `refgenie` to automatically pull down the standard matched FASTA, GTF, and tRNAdb resources. Keep the manual URL definitions *only* for custom or bespoke reference genomes.
