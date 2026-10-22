# Reference Management

The Reference Management module ensures that the pipeline has access to all necessary genomic resources, including FASTA files, gene annotations, and specialized repeat/tRNA databases.

All paths below are relative to `globals.results_folder`.

## Automated Downloads

For `mm10` and `mm39` the genome FASTA and GTF are fetched with
[Refgenie](http://refgenie.databio.org/) when `fasta_path`/`gtf_path` are not set.

1. **Initialize**: a local Refgenie config is created at `references/refgenie/genome_config.yaml`.
2. **Pull**: the pipeline pulls the `fasta` asset and the `ensembl_gtf` asset (falling
   back to `gencode_gtf`), then writes `references/<label>/fasta.fa` and
   `references/<label>/annotation.gtf`.

The TE and tRNA annotations do not come from Refgenie. They are fetched for `mm10` and
`mm39` whether or not you use local FASTA/GTF files:

- **RepeatMasker**: `workflow/scripts/get_rmsk.py` downloads the `rmsk` track from the
  UCSC Genome Browser REST API (`api.genome.ucsc.edu`), chromosome by chromosome (only
  `selected_chromosomes` if set), and writes `references/rmsk/<label>.gtf` and `.bed`.
- **GtRNAdb**: the tRNA set is downloaded from `gtrnadb.ucsc.edu` (e.g.
  `mm10-tRNAs.tar.gz`) into `references/gtrnadb/<label>/`. With `selected_chromosomes`
  it is first filtered to those chromosomes.

You cannot supply your own RepeatMasker or GtRNAdb files. The 1.x keys for that
(`rmsk_path`, `gtrnadb_url`, `gtrnadb_bed`) and for remote references (`fasta_url`,
`fasta_link`, `gtf_url`, `gaf_url`) are rejected at start-up; see
[Migrating from 1.x](../migration/index.md#genome).

## Manual Overrides

To use local files, set `fasta_path`, `gtf_path` and `annotation_type` together (the
schema rejects a config with only some of them). Files may be gzip-compressed (`.gz`).
The pipeline writes decompressed, formatted copies (restricted to `selected_chromosomes`
if set) under `references/`, named after the input file without `.gz`, and builds the
STAR index from those.

!!! important "Label Requirement"
    `label` is required even when overriding paths; it is used in reference file names.
    With custom files any label of 1–64 letters, digits, `.`, `_` or `-` is accepted. Only
    `mm10` and `mm39` have **RepeatMasker** and **GtRNAdb** annotations fetched
    automatically: keep one of those labels if your files are that assembly. Any other
    label requires `disable_TE_analysis: true` and `disable_tRNA_analysis: true`, and
    SalmonTE is skipped unless the label starts with `mm`, `hg`, `dr` or `dm`
    (see [SalmonTE supported genomes](../configuration/index.md#salmonte-supported-genomes)).

```yaml
genome:
  label: "mm10" # Required; use mm10/mm39 if your FASTA/GTF are that assembly
  fasta_path: "/path/to/my_genome.fa"
  gtf_path: "/path/to/my_annotation.gtf"
  annotation_type: "ensembl" # Must be one of: ensembl, gencode, mgi
```

## Chromosome Subsetting

The pipeline can filter your FASTA and GTF files to include only a specific set of chromosomes.

```yaml
genome:
  selected_chromosomes: ["chr19", "chrX"]
```

!!! warning "Mostly for Testing"
    This feature is primarily intended for **testing purposes** or for focusing on specific small genomic regions. Subsetting the genome can significantly simplify the analysis for debugging but may lead to mapping biases if used in production without caution.

## STAR index

The index is built once per results folder from the formatted FASTA and GTF with
`--sjdbOverhang 100` and `--genomeSAindexNbases min(14, log2(genome length)/2 - 1)`, as
the STAR manual recommends. The second value is 14 for mammalian genomes and smaller for
subsets such as a single chromosome, so `selected_chromosomes` runs do not need a
full-size index. The rule requests 8 threads and 256 GB of memory.

## Results

| Location | Description |
| :--- | :--- |
| `references/` | Root directory for all genomic resources, under `globals.results_folder`. |
| `references/<label>/fasta.fa`, `annotation.gtf` | Formatted genome and annotation (automatic downloads). |
| `references/<file name>` | Formatted copies of `fasta_path` / `gtf_path` (local files). |
| `references/refgenie/genome_config.yaml` | Refgenie configuration. |
| `references/STAR/` | The generated STAR genome index. |
| `references/rmsk/` | RepeatMasker annotations (`<label>.gtf`, `<label>.bed`). |
| `references/gtrnadb/<label>/` | GtRNAdb tRNA sequences and annotations. |
