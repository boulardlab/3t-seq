# Configuration instructions

See below for an example config file with explanation of each option and description of common use-cases.

## A complete example

```yaml
# config/config.yaml

# A list of datasets 
# Every dataset is defined by a name, a path to a sample sheet, trimmomatic, star and bamCoverage options.
# All these options are mandatory.
sequencing_libraries:
  - name: GSE13073
    sample_sheet: sample-sheet.csv
    trimmomatic: >-
      "ILLUMINACLIP:TruSeq3-PE.fa:1:0:15:2
       SLIDINGWINDOW:20:22
       MAXINFO:20:0.6
       LEADING:22
       TRAILING:20
       MINLEN:75"
    star: >-
      "--seedSearchStartLmax 30
       --outFilterMismatchNoverReadLmax 0.04
       --winAnchorMultimapNmax 40"
    bamCoverage: "--binSize 50 --normalizeUsing None"

#   - name: ...
#     sample_sheet: ...
#     trimmomatic: ...
#     star: ...
#     bamCoverage: ...

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

  # path to qc
  qc_folder: results/qc

  # path to log
  log_folder: results/log

  # path to references
  references_folder: results/references

  # temp folder
  tmp_folder: /tmp

  # path to analysis
  analysis_folder: results/analysis

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

The pipeline resolves reads paths starting from two bits of information:
1. `reads_folder` in the `globals` sections
2. The name of a library in `sequencing_libraries` list of objects

In the example configuration above, `reads_folder: .` and `sequencing_libraries[0].name: GSE13073`. These resolve to `./GSE13073`. It is crucial that **this folder exists before starting the pipeline**. This is because in this folder, the pipeline will look for input files.

Another example:
```yaml
sequencing_libraries:
  - name: first-batch
    sample_sheet: sample-sheet-first-batch.csv
    # [...]
  
  - name: second-batch
    sample_sheet: sample-sheet-second-batch.csv
    # [...]

globals:
  reads_folder: reads
  # [...]
```

In this scenario, 3t-seq will look for the `reads` folder and inside of it will look for two folders names `first-batch` and `second-batch`: `reads/first-batch` and `reads/second-batch`. If any of the two is not detected, the pipeline will fail.

### Naming convetion

Reads files need to have one of the following extensions: `fq`, `fq.gz`, `fastq`, `fastq.gz`. For a given sequencing library, the pipeline expects files to have the same extension.

For paired-end reads, the two mates should have one of the following idenfiers **before** the extension: `(_1, _2)`, `(_R1, _R2)`, `(_1_sequence, _2_sequence)`.


## How to use local reference files

The `references_folder` can be outside of `results_folder`. For instance:
```yaml
globals:
  # [...]
  # path to results folder
  results_folder: results/

  # path to references
  references_folder: /path/to/references
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