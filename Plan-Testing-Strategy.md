# Implementation Plan: Concrete Testing Strategy

## Goal Description
Robust pipelines require automated testing to prevent regressions when dependencies update or rules are refactored. Currently, 3t-seq provides a `tests/` directory for manual execution, but lacks an automated continuous integration (CI) mechanism that asserts the *correctness* of the output reliably. This plan details the implementation of a concrete, automated testing strategy, leveraging Snakemake's native testing APIs, comprehensive test datasets, linting, and integration into a CI platform like GitHub Actions.

## Proposed Changes

### 1. Standardized Expected Outputs via Snakemake API
Instead of custom Bash scripts and MD5 manual hashing, we will utilize Snakemake's native `--generate-unit-tests` feature to automatically generate pytest-based unit tests for each pipeline rule. This method formally guarantees that refactoring steps do not silently alter the pipeline's computational results, by asserting that the outputs of individual rules mach pre-recorded expected outputs.

- We will generate test scaffolds stored under `.tests/unit/`.
- This ensures test generation is driven entirely by the Snakemake graph, significantly reducing maintenance overhead compared to custom bash harness validation.

### 2. Test Dataset Generation
To ensure reproducibility without inflating repository size or CI runtime, we will generate and version-control a minimal, representative test dataset.

- **Data Source:** We will take one representative 3t-seq raw FASTQ file pair from existing data (e.g., an existing standard wild-type sample).
- **Subsampling Method:** We will randomly subsample ~10,000 to 50,000 reads using `seqtk sample` (e.g., `conda run -n snakemake-8.9.0 seqtk sample -s100 real_sample.fq.gz 10000 > test_sample.fq`).
- **Storage:** The resulting downsampled fastq files will be committed to `tests/data/raw/` in the repository alongside a minimal genome index (`tests/data/genome/` or by downloading a canonical small reference like yeast or a single human chromosome during the test run). We will explicitly define the genomic features we expect to intersect these reads.

### 3. Supporting Pipeline Development
A testing harness isn't useful if it's painful to maintain. When pipeline parameters change (e.g., a new hardcoded strand flag or a threshold is altered) which *intentionally* changes the correct expected output, developers need a way to easily update the expected outputs.

**Updating test cases:**
Developers can regenerate the full expected test suite outputs seamlessly by re-running the Snakemake test generator over the test dataset after validating the new output:
```bash
# Delete old test definitions representing old behavior
rm -rf .tests/unit/

# Generate new expected outputs and pytest asserts
conda run -n snakemake-8.9.0 snakemake --directory tests --use-conda --generate-unit-tests
```
This single command explicitly replaces manual MD5 file adjustments and supports seamless pipeline evolution.

### 4. Usability for Local, Faster Development
The testing infrastructure is built to easily run locally so developers can catch bugs prior to pushing commits.

Instead of needing a custom shell script like `tests/run_tests.sh`, developers can execute tests rapidly using `pytest`:
```bash
# Run all unit tests locally using the snakemake-8.9.0 environment
conda run -n snakemake-8.9.0 pytest .tests/unit
```
To further smooth the local developer experience, we will add a trivial `Makefile` to the root:
```makefile
.PHONY: test update-tests

test:
	conda run -n snakemake-8.9.0 pytest .tests/unit

update-tests:
	rm -rf .tests/unit/
	conda run -n snakemake-8.9.0 snakemake --directory tests --use-conda --generate-unit-tests
```
This enables developers to type `make test` and instantly validate changes against the expected data locally.

### 5. Continuous Integration Setup
We will automate the `pytest` suite to run on every code change submitted to the repository via GitHub Actions.

#### [NEW] `.github/workflows/ci.yml` (or equivalent CI definition)
- Define a GitHub Actions workflow that:
  - Checks out the repository.
  - Installs Micromamba and caches the `snakemake-8.9.0` environment.
  - Generates or executes the tests via `pytest .tests/unit` natively.

**Example GitHub Actions Workflow snippet (`ci.yml`):**
```yaml
      - name: Setup environment
        uses: mamba-org/setup-micromamba@v1
        with:
          environment-file: environment.yaml # or wherever snakemake-8.9.0 is defined
          environment-name: snakemake-8.9.0
          cache-environment: true

      - name: Run PyTest on Snakemake unit tests
        run: micromamba run -n snakemake-8.9.0 pytest .tests/unit
```

### 6. Code Quality and Linting
To maintain readability and catch syntax issues early, we will incorporate linting into the test pipeline.

- Include a step to run `conda run -n snakemake-8.9.0 snakemake --lint` to ensure the `Snakefile` adheres to formatting standards.
- Incorporate R and Python linters for custom scripts inside `workflow/scripts` (e.g., `lintr` for R scripts).

## Verification Plan

### Automated Tests
- Push a branch containing the `.github/workflows/ci.yml`.
- Ensure the GitHub Action triggers immediately and successfully runs `pytest .tests/unit` to green completion.
- Review the CI action logs to confirm that environments are correctly cached.

### Manual Verification
- **Simulate a Regression:** Introduce an intentional bug locally (e.g., modify parameters in the `build_tRNA_coverage_matrix` R script).
- Run `make test` (or `conda run -n snakemake-8.9.0 pytest .tests/unit`) locally.
- Verify that the `pytest` output successfully flags the mismatched file assertion, proving the efficacy of the test harness.
- **Simulate an Update:** Revert the bug, update the parametrization intentionality, and run `make update-tests` to ensure the test scaffolds properly refresh and record the newly anticipated output state.
