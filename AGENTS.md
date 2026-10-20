# Agent Instructions for 3t-seq

This document provides comprehensive context and guidelines for AI agents working on the 3t-seq repository. Adherence to these guidelines is mandatory for all coding, modification, or testing tasks.

## Project Overview
**3t-seq** is a Snakemake workflow for integrated analysis of single-copy genes, transposable elements (TEs), and tRNAs from total RNA-seq data. The core objective is to ensure reproducible and high-quality bioinformatics analysis.

## Tech Stack
- **Workflow Engine:** Snakemake (v8.9.0)
- **Environment Management:** [Pixi](https://pixi.sh)
- **Languages:** Python, R, Bash
- **Containerization:** Apptainer/Singularity (configured via Snakemake profiles). Pixi does not install it: use a system or module Apptainer >= 1.2.5, never the conda-forge one (its `squashfuse_ll` hangs on SIF reads, #39).

## Always prefer using `pixi run <task>` to ensure the correct environment and dependencies are used. The project uses a parameterized task system and multiple environments.

### Pixi Environments
- **`default`**: The base environment containing Snakemake and core dependencies. Used for production runs.
- **`dev`**: Extends `default` with development tools (`snakefmt`, `pre-commit`, `pytest`). Use this for linting and formatting.
- **`docs`**: Extends `dev` with documentation tools (`mkdocs`, plugins). Use this for building and serving the manual.

To specify an environment, use the `-e` (or `--environment`) flag:
```bash
pixi run -e dev format
pixi run -e docs docs-serve
```

| Task | Command | Description |
| :--- | :--- | :--- |
| **Setup** | `pixi run -e dev setup` | Initializes Git LFS and pulls large test assets (FASTQs, etc.). |
| **Linting** | `pixi run lint <config> <profile>` | Runs `snakemake --lint`. Defaults: `remote-references`, `laptop`. |
| **Formatting** | `pixi run -e dev format` | Formats Snakemake files using `snakefmt`. |
| **Format check** | `pixi run -e dev format-check` | `snakefmt --check --diff`; what `lint.yml` runs. Does not modify files. |
| **Dry Run** | `pixi run dry-run <config> <profile>` | Snakemake dry-run (`-n`). CI `smoke` runs `pixi run -e dev dry-run <config> ci -- --sdm conda` for all three configs. |
| **Test** | `pixi run -e dev test <config> <profile>` | Executes the integration test. Needs `-e dev`: it depends on `setup`, which only exists in `dev`/`docs`. CI `e2e` runs `local-references ci`. |
| **Report** | `pixi run -e dev test-report <config> <profile>` | Generates a Snakemake report (depends on `test`). |
| **Containers** | `pixi run containerize` | Regenerates the `Dockerfile` and updates the 3t-seq image pin in `workflow/Snakefile`. Run after any conda env change and commit both files. |
| **Pre-build** | `pixi run build-conda-envs` | Pre-creates all Conda environments. |

**Available Arguments (Positional):**
- `config`: `remote-references`, `local-references`, `adaptive-trim`
- `profile`: `laptop`, `hpc`, `ci` (GitHub runner: 4 cores, 2 jobs)

Extra Snakemake arguments go after `--`, e.g. `pixi run -e dev test local-references laptop -- --cores 4`; without it pixi fails with "received more arguments than expected".

## Repository Structure
- `workflow/`: Contains the core logic.
    - `Snakefile`: Main entry point.
    - `rules/`: Modular Snakemake rules.
    - `scripts/`: R and Python scripts called by rules.
    - `envs/`: Conda environment definitions for specific rules.
- `config/`: Default configuration files and schema definitions.
- `.tests/`: Integration and unit tests.
    - `integration/`: Small dataset for end-to-end testing.
- `pixi.toml`: Project-level dependency and task management.

## Development Guidelines for Agents

### 1. Execution Flow
- **Environment First:** Always prefix commands with `pixi run` unless explicitly modifying a Snakefile or Bash script that does not rely on a specific environment setup.
- **Formatting:** Run `pixi run -e dev format` before committing. CI (`lint.yml`) runs `pixi run -e dev format-check`, which fails on any `snakefmt` diff instead of fixing it.
- **Validation:** Before committing changes to workflow logic or core scripts, run what the required CI job `smoke` runs: `pixi run -e dev lint <config> laptop` and `pixi run -e dev dry-run <config> ci -- --sdm conda` for `remote-references`, `local-references` and `adaptive-trim`.
- **Documentation:** Keep documentation in sync with feature changes. Any agent modifying source code MUST check the `docs/` directory and update relevant parts. Validate changes with `pixi run -e docs docs-build`.

### 2. Code Style Guidelines (Python/R/Bash)
- **Python:**
    - **Typing:** Use comprehensive type hinting for all function signatures and complex data structures.
    - **Naming:** Follow snake_case for functions and variables, PascalCase for classes. Modules should be descriptive.
    - **Imports:** Use absolute imports where possible. Group imports logically (standard library, third-party, local).
    - **Error Handling:** Use specific exceptions rather than generic ones. Handle I/O errors explicitly. Log errors using the standard logging mechanism defined in `scripts/`.
    - **Efficiency:** Prioritize time/space complexity for data processing pipelines.
- **R:**
    - **Style:** Follow standard R style conventions for vectorization and function definitions.
    - **Dependencies:** Explicitly define R dependencies in `envs/` or relevant Snakefile sections.
- **Bash:**
    - **Quoting:** Always quote file paths and arguments containing spaces.
    - **Scripting:** Scripts must be idempotent where possible. Use clear, descriptive variable names.
- **General:**
    - **Security:** Never handle secrets or keys directly in code. Use environment variables or secure configuration files.
    - **Documentation:** Functions and modules must have clear docstrings explaining parameters, return values, and potential exceptions.

### 3. Standardized Patterns & Best Practices

#### Configuration Management
- **Single Source of Truth:** Always define default values in `workflow/schemas/config.schema.yaml`.
- **Enforcement:** `workflow/Snakefile` calls `validate_config()` (defined in `workflow/rules/common.smk`), which validates against the schema, reports all errors, then applies the defaults via `apply_schema_defaults()` and runs `check_genome_support()`. Do not call `apply_schema_defaults()` separately. Most schema objects set `additionalProperties: false`, so every new key must be declared in the schema.
- **Resolution:** Use `get_resolved_param()` for parameters that can be defined at both the 'library' and 'global default' levels.

#### File Naming & Formats
- **Tabular Data:** Use `.csv` instead of `.txt` for output tables (e.g., DESeq2 results). This ensures tools like `datavzrd` and `pandas` can automatically detect the format without manual intervention.
- **Consistency:** Ensure file extensions match the content (e.g., `.fastq.gz` for compressed reads).

#### Cross-Platform Compatibility
- **macOS Coreutils:** Pixi puts `$PIXI_PROJECT_ROOT/.pixi/macos-shims` first on `PATH` (osx-arm64). Provide GNU tools there as plain-named symlinks to Homebrew `coreutils`/`gnu-sed` `gnubin` (see `docs/getting-started/index.md`, "Apple Silicon & macOS notes"). Do not manually override system binaries.
- **Temporary Directories:** In rule `shell:` blocks, use `{resources.tmpdir}` (e.g. `mktemp -d -p "{resources.tmpdir}"`), not `$TMPDIR`: Snakemake runs shells with `set -u` and `$TMPDIR` is often not propagated into containers ("TMPDIR: unbound variable"). Snakemake always sets `resources.tmpdir`; profiles may override it via `default-resources`.

#### Containerization & CI
The full reference (branch/PR model, each workflow and its failure modes, image tags, the `main` ruleset, release procedure, troubleshooting) is the Developer Guide in `docs/developer/`; start at `docs/developer/index.md`. Keep this section a summary and update the guide when CI or container behaviour changes.

- **Multi-arch Support:** `build-3t-seq.yml` and `build-salmonte.yml` build `linux/amd64` and `linux/arm64`. Every dependency in `workflow/envs/*.yaml` must be solvable for both, or the image build fails. CI runs the pipeline (`e2e`) on `linux/amd64` only.
- **SHA Tagging:** Images pushed to GHCR should be tagged with the commit SHA in addition to the branch/tag name.
- **Tests in CI:** `test.yml` job `smoke` (lint + conda-only dry-run of all three configs, ~30 s) is a required check on `main`. Job `e2e` (full `local-references` run with the `ci` profile, 90 min limit, Apptainer from `ppa:apptainer/ppa`) runs on release-please PRs, nightly and via `gh workflow run test.yml --ref <branch>`; it writes a per-rule timing table to the job summary. Required checks on `main`: `check`, `dockerfile-sync`, `smoke`.
- **Pinned Images:** Every `ghcr.io/boulardlab/*` image under `workflow/` must use an immutable tag (never `:latest` or untagged) that exists in GHCR. The `container-check.yml` gate enforces this on PRs to `main`/`dev` and on pushes to `main`.
- **3t-seq Image (env hash):** Snakemake looks conda envs up inside the `containerized:` image by content hash and fails if one is missing, so the image must match the envs of every commit.
    - The image is tagged `env-<conda_env_hash>`, using the label `snakemake --containerize` writes into the `Dockerfile`. `workflow/Snakefile` pins that tag, managed by `.github/scripts/container_pin.py`.
    - After changing anything in a conda env, run `pixi run containerize` and commit `Dockerfile` + `workflow/Snakefile` together. CI (`container-check.yml`, job `dockerfile-sync`) fails if either is stale.
    - `build-3t-seq.yml` builds each `env-<hash>` once (PRs from this repo, and pushes to `main`/`dev`) and never overwrites it. If the hash already exists, branch/`sha-*`/`latest` tags are pointed at it without rebuilding. `latest` moves only on pushes to `main`. A `v*` tag adds an `:X.Y.Z` alias.
    - PRs from forks cannot push images. A maintainer must build the hash: push the branch to this repo, then run `gh workflow run build-3t-seq.yml --ref <branch>` (a push alone does not trigger the build; push triggers are `main`/`dev` only) or open a same-repo PR from it.
- **SalmonTE Image (release version):** `workflow/rules/salmonTE.smk` pins `salmonte:X.Y.Z` (`x-release-please-version` markers). `build-salmonte.yml` builds `:X.Y.Z` from the release-please PR before it merges, or manually with the `version` input to bootstrap. It refuses to overwrite once git tag `vX.Y.Z` exists.
- **Releases:** release-please (`release-please-config.json`, `.release-please-manifest.json`) bumps `pixi.toml` and the SalmonTE pins. It needs the `RELEASE_PLEASE_TOKEN` secret (a PAT or App token), because PRs and tags created with `GITHUB_TOKEN` do not trigger other workflows.

### 5. Agent Specific Rules
- **Cursor Rules:** [Check for .cursor/rules/ or .cursorrules for specific agent instructions.]
- **Copilot Rules:** [Check for .github/copilot-instructions.md for specific Copilot guidelines.]

## Common Gotchas
- **Slurm/HPC:** Be mindful of resource allocations (`threads`, `mem_mb`) in rules, as the pipeline is designed for HPC execution.
- **Profiles:** Users typically run the pipeline with `--profile profile/default`.