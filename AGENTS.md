# Agent Instructions for 3t-seq

This document provides context and guidelines for AI agents working on the 3t-seq repository.

## Project Overview
**3t-seq** is a Snakemake workflow for integrated analysis of single-copy genes, transposable elements (TEs), and tRNAs from total RNA-seq data.

## Tech Stack
- **Workflow Engine:** Snakemake (v8.9.0)
- **Environment Management:** [Pixi](https://pixi.sh)
- **Languages:** Python, R, Bash
- **Containerization:** Apptainer/Singularity (configured via Snakemake profiles)

## Key Workflows & Commands
Always prefer using `pixi run` to ensure the correct environment is used.

| Task | Command | Description |
|------|---------|-------------|
| **Linting** | `pixi run lint` | Runs `snakemake --lint` on the integration test config. |
| **Formatting** | `pixi run format` | Formats Snakemake files using `snakefmt`. |
| **Dry Run** | `pixi run run-test` | Performs a Snakemake dry-run using the integration test dataset. |
| **Hooks** | `pixi run hooks` | Runs `pre-commit` hooks on all files. |
| **Full Test** | `snakemake --directory .tests/integration ...` | Runs the actual integration test (see README for details). |

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
1. **Use Pixi:** Always use `pixi run <task>` or prefix commands with `pixi run`. For parameters run as `pixi run <task> <param1> <param2>`, for example `pixi run run-test remote-references laptop`.
2. **Follow Formatting:** Run `pixi run format` before submitting changes to Snakemake files.
3. **Lint Early:** Run `pixi run lint` to catch workflow errors before execution.
4. **Testing:** Refer to [Plan-Testing-Strategy.md](file:///Users/francesco/scratch/3t-seq/main.worktrees/dev/worktree-2026-03-31T14-14-10/Plan-Testing-Strategy.md) for the long-term testing goals, including unit test generation.
5. **Pathing:** When referencing files in the repo, use relative paths from the root.

## Common Gotchas
- **Slurm/HPC:** The pipeline is designed to run on HPCs. Be mindful of resource allocations (`threads`, `mem_mb`) in rules.
- **Profiles:** Users typically run the pipeline with `--profile profile/default`.
