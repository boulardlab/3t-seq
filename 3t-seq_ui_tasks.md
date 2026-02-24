# 3t-seq UI Simplification Tasks

Based on the user interface analysis, here is a precise, actionable list of technical tasks to simplify the 3t-seq pipeline configuration and improve user adoption.

## 0. Preparation
- [ ] **Task 0.1:** Create a new git branch from the `dev` branch to develop these features.

## 1. Auto-derive Results Subdirectories
**Objective:** Remove the need for the user to manually specify internal pipeline folders.
- [ ] **Task 1.1:** Modify the Snakefile's configuration parsing (or `common.smk` if used) to automatically define `qc_folder`, `log_folder`, `references_folder`, `tmp_folder`, and `analysis_folder` based on the provided `results_folder`.
- [ ] **Task 1.2:** Remove `qc_folder`, `log_folder`, `references_folder`, `tmp_folder`, and `analysis_folder` from the configuration JSON/YAML schema entirely (create the schema if it doesn't already exist to enforce this).
- [ ] **Task 1.3:** Align the example `config.yaml` files in the `config/` directory to the newly defined pipeline configuration requirements (i.e., by removing these global folders).
- [ ] **Task 1.4:** Update `config/README.md` to explain that these folders are now handled automatically by the pipeline.

## 2. Implement Sensible Defaults for Tool Parameters
**Objective:** Abstract away complex command-line strings for `trimmomatic`, `star`, and `bamCoverage`, making them optional.
- [ ] **Task 2.1:** Define a standard set of default parameters for `trimmomatic`, `star`, and `bamCoverage` directly within the Snakefile rules. For `trimmomatic`, smartly infer the default parameters from the sequencing library properties (e.g., single-end vs paired-end) instead of relying on static strings. Do not use an external `defaults.yaml` file.
- [ ] **Task 2.2:** Update the configuration schema for `sequencing_libraries` to make `trimmomatic`, `star`, and `bamCoverage` fields optional rather than required.
- [ ] **Task 2.3:** Update the Snakefile rules to use `config["sequencing_libraries"][X].get("trimmomatic", default_trimmomatic_params)` to fall back to defaults if the user omits them.
- [ ] **Task 2.4:** Remove these verbose parameters from the default `config.yaml` templates and update the `config/README.md` to demonstrate how they can be overridden if needed.

## 3. Enhance Sample Sheet Parsing
**Objective:** Allow flexible file naming and absolute paths, removing the strict folder structure and suffix requirements.
- [ ] **Task 3.1:** Locate the Python function responsible for parsing the `sample_sheet` and determining input reads (usually in `Snakefile` or a helper script).
- [ ] **Task 3.2:** Modify the path resolution logic: If the paths in `filename`, `filename_1`, or `filename_2` are absolute (e.g., start with `/`), use them directly instead of prepending `reads_folder/library_name`.
- [ ] **Task 3.3:** Remove the strict suffix validation logic (`_1`, `_R1`) from the pipeline's wildcards/inputs if explicit file paths are defined in the sample sheet. Ensure this does not break filename parsing in `common.smk` (e.g., skip wildcard-based regular expression inference or adjust the regex to support arbitrary base names) so that the pipeline trusts the precise paths from the DataFrame.
- [ ] **Task 3.4:** Add a clear validation loop during DAG generation: Use `os.path.exists()` to verify all absolute/relative paths extracted from the sample sheet, raising an immediate, informative `RuntimeError` or `WorkflowError` if files cannot be found before any rules execute.
- [ ] **Task 3.5:** Update documentation to explain the new, flexible sample sheet formats (supporting arbitrary filenames and absolute paths).

## 4. Integrate Built-in Genome Profiles via Galaxy Data Cache
**Objective:** Leverage the Galaxy Reference Data project (`refgenie` or CVMFS) to automatically resolve reference genomes (e.g., `mm10`, `hg38`) without manual URL entry.
- [ ] **Task 4.1:** Investigate and integrate the Galaxy reference data cache (https://galaxyproject.org/admin/reference-data-repo/) as the backend source for standard genome FASTA, GTF, and GtRNAdb files.
- [ ] **Task 4.2:** Update configuration parsing: If the user provides a recognized `label` (e.g., `mm10`), dynamically pull the reference paths/URLs using the Galaxy data repository instead of custom dictionaries.
- [ ] **Task 4.3:** Update the `genome` object in the configuration schema to make `fasta_url`, `gtf_url`, and `gtrnadb_url` optional when using standardized genomes.
- [ ] **Task 4.4:** Refactor the `references` downloading or symlinking rules to natively interface with the Galaxy data cache.
- [ ] **Task 4.5:** Document the integration of the Galaxy reference data repository in `config/README.md` and instruct users on how to use standard genome labels.
