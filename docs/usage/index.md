# Usage

The 3t-seq pipeline leverages Pixi to simplify execution and ensure environment consistency.

## Running the Pipeline

To execute the pipeline, use the `pixi run test` command or call `snakemake` directly via `pixi run`.

### Using Pixi Tasks

We have predefined several tasks to help you manage the pipeline:

- **Environment Setup**:
  Initializes Git LFS and pulls large data files needed for integration tests and some analysis steps.

  ```bash
  pixi run setup
  ```

- **Linting (Spell-checking)**:

  Linting checks your configuration files for typos or formatting errors before you run anything, preventing the pipeline from crashing hours later.

  ```bash
  pixi run lint remote-references laptop
  ```

- **Dry-run (Rehearsal)**:
  A dry-run tells the pipeline to plan out all the steps and show you exactly what files it *would* create and what jobs it *would* run, without actually doing any of the heavy computations.

  ```bash
  pixi run dry-run remote-references laptop
  ```

- **Execution (Run)**:
  This launches the actual data processing.

  ```bash
  pixi run test remote-references laptop
  ```

### Advanced Usage

You can pass extra arguments to Snakemake through Pixi:

```bash
pixi run snakemake --cores 8 --use-conda --directory .tests/integration
```

## Reports

After the pipeline completes, you can generate a comprehensive HTML report:

```bash
pixi run test-report remote-references laptop
```

This will create a `report.zip` file which you can extract to view results in your browser.
