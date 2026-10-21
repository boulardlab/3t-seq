# Reference Management

The Reference Management module ensures that the pipeline has access to all necessary genomic resources, including FASTA files, gene annotations, and specialized repeat/tRNA databases.

## Automated Downloads (Refgenie)

For standard genomes like `mm10` and `mm39`, 3t-seq uses [Refgenie](http://refgenie.databio.org/) to automatically fetch high-quality references.

### Workflow

1. **Initialize**: A local Refgenie repository is created in `results/references/refgenie/`.
2. **Pull**: The pipeline pulls the `fasta` and `ensembl_gtf` (or `gencode_gtf`) assets.
3. **Specialized Assets**:
    - **RepeatMasker**: Fetched and converted to GTF/BED for TE analysis.
    - **GtRNAdb**: Custom tRNA sequences and annotations are downloaded directly from the GtRNAdb servers.

## Manual Overrides

If you are working with a custom genome assembly or want to use specific local files, you can override the automated downloads.

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

## Results

| Location | Description |
| :--- | :--- |
| `results/references/` | Root directory for all genomic resources. |
| `results/references/STAR/` | The generated STAR genome index. |
| `results/references/rmsk/` | RepeatMasker annotations (GTF/BED). |
| `results/references/gtrnadb/` | tRNA sequences and annotations. |
