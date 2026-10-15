# Agent Instructions for 3t-seq

This document provides comprehensive context and guidelines for AI agents working on the 3t-seq repository. Adherence to these guidelines is mandatory for all coding, modification, or testing tasks.

## Project Overview
**3t-seq** is a Snakemake workflow for integrated analysis of single-copy genes, transposable elements (TEs), and tRNAs from total RNA-seq data. The core objective is to ensure reproducible and high-quality bioinformatics analysis.

## Tech Stack
- **Workflow Engine:** Snakemake (v8.9.0)
- **Environment Management:** [Pixi](https://pixi.sh)
- **Languages:** Python, R, Bash
- **Containerization:** Apptainer/Singularity (configured via Snakemake profiles)

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
| **Setup** | `pixi run setup` | Initializes Git LFS and pulls large test assets (FASTQs, etc.). |
| **Linting** | `pixi run lint <config> <profile>` | Runs `snakemake --lint`. Defaults: `remote-references`, `laptop`. |
| **Formatting** | `pixi run format` | Formats Snakemake files using `snakefmt`. |
| **Dry Run** | `pixi run dry-run <config> <profile>` | Snakemake dry-run (`-n`). |
| **Test** | `pixi run test <config> <profile>` | Executes the integration test (depends on `setup`). |
| **Report** | `pixi run test-report <config> <profile>` | Generates a Snakemake report (depends on `test`). |
| **Containers** | `pixi run containerize` | Regenerates the `Dockerfile` and updates the 3t-seq image pin in `workflow/Snakefile`. Run after any conda env change and commit both files. |
| **Pre-build** | `pixi run build-conda-envs` | Pre-creates all Conda environments. |

**Available Arguments (Positional):**
- `config`: `remote-references`, `local-references`, `adaptive-trim`
- `profile`: `laptop`, `hpc`

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
- **Formatting Prerequisite:** Before any code change, run `pixi run format` to ensure code style consistency.
- **Validation:** Run `pixi run lint` before attempting to commit any changes that modify workflow logic or core scripts.
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
- **Enforcement:** Use `apply_schema_defaults(config, "workflow/schemas/config.schema.yaml")` in `workflow/rules/common.smk` to ensure these defaults are correctly populated at runtime.
- **Resolution:** Use `get_resolved_param()` for parameters that can be defined at both the 'library' and 'global default' levels.

#### File Naming & Formats
- **Tabular Data:** Use `.csv` instead of `.txt` for output tables (e.g., DESeq2 results). This ensures tools like `datavzrd` and `pandas` can automatically detect the format without manual intervention.
- **Consistency:** Ensure file extensions match the content (e.g., `.fastq.gz` for compressed reads).

#### Cross-Platform Compatibility
- **macOS Coreutils:** On macOS, the environment is shimmed to use GNU coreutils (e.g., `gln` as `ln`). This is handled automatically by the Pixi activation script (`scripts/macos-coreutils-setup.sh`). Do not manually override system binaries.
- **Temporary Directories:** Respect the `$TMPDIR` environment variable. Use centralized logic for temporary directory assignment to ensure compatibility with HPC nodes and local machines.

#### Containerization & CI
- **Multi-arch Support:** Docker builds must support both `linux/amd64` and `linux/arm64`. Verify that all dependencies in `workflow/envs/*.yaml` are solvable for both architectures.
- **SHA Tagging:** Images pushed to GHCR should be tagged with the commit SHA in addition to the branch/tag name.
- **Pinned Images:** Every `ghcr.io/boulardlab/*` image under `workflow/` must use an immutable tag (never `:latest` or untagged) that exists in GHCR. The `container-check.yml` gate enforces this on PRs to `main`/`dev` and on pushes to `main`.
- **3t-seq Image (env hash):** Snakemake looks conda envs up inside the `containerized:` image by content hash and fails if one is missing, so the image must match the envs of every commit.
    - The image is tagged `env-<conda_env_hash>`, using the label `snakemake --containerize` writes into the `Dockerfile`. `workflow/Snakefile` pins that tag, managed by `.github/scripts/container_pin.py`.
    - After changing anything in a conda env, run `pixi run containerize` and commit `Dockerfile` + `workflow/Snakefile` together. CI (`container-check.yml`, job `dockerfile-sync`) fails if either is stale.
    - `build-3t-seq.yml` builds each `env-<hash>` once (PRs from this repo, and pushes to `main`/`dev`) and never overwrites it. If the hash already exists, branch/`sha-*`/`latest` tags are pointed at it without rebuilding. `latest` moves only on pushes to `main`. A `v*` tag adds an `:X.Y.Z` alias.
    - PRs from forks cannot push images. A maintainer must build the hash (push the branch to this repo or run the workflow manually).
- **SalmonTE Image (release version):** `workflow/rules/salmonTE.smk` pins `salmonte:X.Y.Z` (`x-release-please-version` markers). `build-salmonte.yml` builds `:X.Y.Z` from the release-please PR before it merges, or manually with the `version` input to bootstrap. It refuses to overwrite once git tag `vX.Y.Z` exists.
- **Releases:** release-please (`release-please-config.json`, `.release-please-manifest.json`) bumps `pixi.toml` and the SalmonTE pins. It needs the `RELEASE_PLEASE_TOKEN` secret (a PAT or App token), because PRs and tags created with `GITHUB_TOKEN` do not trigger other workflows.

### 5. Agent Specific Rules
- **Cursor Rules:** [Check for .cursor/rules/ or .cursorrules for specific agent instructions.]
- **Copilot Rules:** [Check for .github/copilot-instructions.md for specific Copilot guidelines.]

## Common Gotchas
- **Slurm/HPC:** Be mindful of resource allocations (`threads`, `mem_mb`) in rules, as the pipeline is designed for HPC execution.
- **Profiles:** Users typically run the pipeline with `--profile profile/default`.