# Configuration instructions

See below for an example config file with explanation of each option.

```yaml
# config/config.yaml

# A list of datasets 
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

# 
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
