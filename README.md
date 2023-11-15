# RNA-seq Analysis Pipeline

## Overview

This is a Snakemake pipeline for the integrated analysis of single copy genes, transposable elements and tRNAs. It performs standard quality control checks and genome alignment in three different ways specialized either for single copy genes or transposable elements. It then quantifies gene expression depending on how the alignement step was performed. Finally it performs differential gene expression analysis yielding lists of genes significantly deregulated between two given conditions.

![Overall 3t-seq workflow](docs/figures/3t-wf.png)

## Requirements


- [Conda](https://conda.io/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Apptainer](https://apptainer.org/docs/user/latest/)

## Usage

### 1. Download the latest version of the pipeline

For example, to download version 1.0.0:

```bash
curl -LJO "https://github.com/boulardlab/3t-seq/archive/v1.0.0.zip"
unzip "3t-seq-1.0.0.zip"
```

### 2. Create a conda environment

```bash
conda create -n snakemake-latest -c bioconda -c conda-forge snakemake singularity
conda activate snakemake-latest
```

### 3. Configure samples and parameters

Edit the `config.yaml` file to specify your sample information and analysis parameters.

### 4. Execute the pipeline

```bash
snakemake --profile profile/default
```

### 5. View results

After the pipeline completes, you can find the results in the `results/` directory.

## Configuration

Adjust parameters in the `config.yaml` file to match your experimental setup.

Here is an example config file:

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
# 	  sample_sheet: ...
# 	  trimmomatic: ...
# 	  star: ...
# 	  bamCoverage: ...

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

### Sample sheet preparation

The sample sheet is a csv file that describe samples metadata:

- The `sample` column reports a human readable name for each sample.
- For pe libraries, `filename_1` and `filename_2` columns report file names for each of the two
sequencing reads mates. For se libraries, `filename` is sufficient. **The pipeline will use these columns to determine if a given dataset was sequenced with pe or se method**.
- The `genotype` column reports the variable of interest. The name of this column is flexible and can be anything as long as you specify what's this name in the config file (in the `deseq2` section).


## Directory Structure

The pipeline will generate an ouput folder tree like so

![Folder tree generated from 3t-seq](docs/figures/folder-tree-tikz.png)

## References

Tabaro F, Boulard M, *3t-seq: automatic gene expression analysis of single copy genes, transposable elements and tRNAs from total RNA-seq data*, Submitted

## License

This project is licensed under the [MIT License](LICENSE).