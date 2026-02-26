# Configuration instructions

See below for an example config file with explanation of each option and description of common use-cases.

## A complete example

```yaml
# config/config.yaml

# A list of datasets 
# Every dataset is defined by a name and a path to a sample sheet.
# The trimmomatic, star, and bamCoverage options are optional and will fall back to sensible defaults if omitted.
sequencing_libraries:
  - name: GSE13073
    sample_sheet: sample-sheet.csv
    # Optional: overrides default parameters for trimmomatic
    trimmomatic: >-
      "ILLUMINACLIP:TruSeq3-PE.fa:1:0:15:2
       SLIDINGWINDOW:20:22
       MAXINFO:20:0.6
       LEADING:22
       TRAILING:20
       MINLEN:75"
    # Optional: overrides default parameters for star
    star: >-
      "--seedSearchStartLmax 30
       --outFilterMismatchNoverReadLmax 0.04
       --winAnchorMultimapNmax 40"
    # Optional: overrides default parameters for bamCoverage
    bamCoverage: "--binSize 50 --normalizeUsing None"

#   - name: ...
#     sample_sheet: ...

# Disable all functionalities related to TE analysis
disable_TE_analysis: false

# Disable tRNA analysis
disable_tRNA_analysis: false

globals:
  # path to reads folder 
  # NB: ./GSE13073 is expected to exist
  reads_folder: .

  # path to results folder
  results_folder: results/

# genome informations
genome:
  # genome label
  label: mm10

  # annotation type
  # can be ensembl, mgi, gencode
  annotation_type: ensembl

  # URL or path to genome sequence
  fasta_url: <Genome fasta URL>
  
  # URL or path to genome annotation file
  gtf_url: <Genome annotation URL>

  # URL to gtRNAdb zip file
  gtrnadb_url: <GtRNADb bundle URL>

# Differential expression analysis parameters
deseq2:
  # wd
  working_directory: ../../..  
  
  # DESeq2 test name, can be Wald or LRT
  test: Wald
  
  # name of the column from sample sheet with experimental variable
  variable: genotype

  # base level from variable column
  reference_level: wt
```

## How 3t-seq resolves reads paths

The pipeline relies primarily on the paths provided in your sample sheet to locate input read files.
For each sample, you can provide paths in the `filename` (for single-end) or `filename_1` and `filename_2` (for paired-end) columns.

These paths can be provided in several flexible formats:
1. **Absolute Paths:** If an absolute path (e.g. `/data/my_experiment/reads/sample1_R1.fq.gz`) is provided, 3t-seq will use it directly.
2. **Relative Paths:** If a relative path is provided, 3t-seq will first attempt to resolve it relative to the current working directory, and then fall back to looking inside `reads_folder / sequence_library_name` for legacy compatibility.

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

## How to use local reference files

You can override the reference files by providing absolute paths in the genome URLs:
```yaml
genome:
  # [...]
  # Provide an absolute path:
  fasta_url: /path/to/references/custom-mm10.fa.gz
  
  # Provide an absolute path:
  gtf_url: /path/to/references/custom-mm10.gtf.gz
```

This allows users to host their own reference files locally and set `genome` informations accordingly

```yaml
genome:
  # [...]
  # This will evaluate to /path/to/references/custom-mm10.fa.gz
  fasta_url: custom-mm10.fa.gz
  
  # This will evaluate to /path/to/references/custom-mm10.gtf.gz
  gtf_url: custom-mm10.gtf.gz
```