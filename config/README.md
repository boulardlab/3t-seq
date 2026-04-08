See below for an example config file with explanation of each option and description of common use-cases. For a complete, exhaustive guide on all available parameters, please refer to the [Full Configuration Documentation](../docs/configuration/index.md).

## Configuration Validation

3t-seq includes a built-in schema validator to ensure your `config.yaml` is correctly formatted before running the full pipeline. You can run the validator using:

```bash
pixi run lint
```

## A complete example

```yaml
# config/config.yaml

# A list of sequencing libraries to process.
# Each item in the list represents a distinct sequencing run or library batch, 
# defined by a unique 'name' and a 'sample_sheet' detailing the samples inside it.
# Tool parameter overrides (trimmomatic, star, bamCoverage, starTE_random, starTE_multihit, strandedness, deseq2) are optional per-library 
# and will fall back to pipeline defaults defined in the schema if omitted.
sequencing_libraries:
  - name: GSE13073
    sample_sheet: sample-sheet.csv
    # Optional library-level overrides:
    deseq2:
      variable: cell_type     # This overrides the global default 'genotype' for this library
```

## How 3t-seq resolves reads paths

The pipeline relies primarily on the paths provided in your sample sheet to locate input read files.
For each sample, you can provide paths in the `filename` (for single-end) or `filename_1` and `filename_2` (for paired-end) columns.

These paths can be provided in several flexible formats:
1. **Absolute Paths:** If an absolute path (e.g. `/data/my_experiment/reads/sample1_R1.fq.gz`) is provided, 3t-seq will use it directly.
2. **Relative Paths:** If a relative path is provided, 3t-seq will resolve it relative to the current working directory.

*Note:* You do not need to provide the file extension (`.fastq.gz`, `.fq.gz`, etc.) if you don't want to. 3t-seq will automatically test standard fastq extensions if the exact path specified doesn't exist.

Because 3t-seq now relies on explicit paths in the sample sheet, **you do not need to rename your raw data files or move them into rigid directory structures**.

### Sample Sheet Format

The sample sheet is a CSV file. It must contain a `name` column representing the sample ID, and either:
- A `filename` column for single-end libraries.
- `filename_1` and `filename_2` columns for paired-end libraries.

*Example (Absolute Paths):*
```csv
name,filename_1,filename_2,genotype
sampleA,/data/raw/batch1/A_read1.fq.gz,/data/raw/batch1/A_read2.fq.gz,WT
sampleB,/data/raw/batch1/B_read1.fq.gz,/data/raw/batch1/B_read2.fq.gz,KO
```

## Built-in Reference Genomes via Refgenie

3t-seq now natively integrates with **Refgenie** to automatically retrieve and manage standard reference genomes (like `mm10` or `hg38`). 

To use this feature, simply provide the `label` in your configuration's `genome` section and omit `fasta_path` and `gtf_path`. The pipeline will automatically fetch the FASTA and GTF files for that genome and symlink them into your results directory. 

*Note: You do not need to configure or set the `$REFGENIE` environment variable on your system; the pipeline handles initialization underneath the hood at runtime.*

## How to use custom local reference files

If you are working with a non-standard genome or prefer to use your own reference files, you can override the built-in profiles by removing the `label` property and providing absolute paths to your local `fasta_path` and `gtf_path`:

```yaml
genome:
  # Remove 'label' to prevent Refgenie lookups
  # label: mm10
  
  # Provide an absolute path to your fasta:
  fasta_path: /path/to/references/custom-mm10.fa
  
  # Provide an absolute path to your gtf:
  gtf_path: /path/to/references/custom-mm10.gtf
  
  # Explicitly specify the annotation type (required for custom references)
  annotation_type: ensembl
  
  # Specify the target species for SalmonTE quantification (e.g. mm, hs, dr, dm):
  salmonte_species: mm
```

These two approaches are mutually exclusive; the configuration validator will alert you if you mix a `label` with explicit `fasta_path` or `gtf_path` parameters.