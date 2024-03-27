# List of rules in alphabetical order

1. `build_trna_coverage_matrix`: Merge per-sample tRNA read counts into a single read counts matrix
1. `coverage`: Calculate single copy genes coverage and outputs Bigwig files
1. `coverage_trna`: Calculate tRNA coverage
1. `datavzrd_single_copy_genes`: generates interactive Volcano and MA plots for single copy genes downstream of DESeq2
1. `datavzrd_starTE_random`: generates interactive Volcano and MA plots for TE (STAR-TE random) downstream of DESeq2
1. `datavzrd_trna`: generates interactive Volcano and MA plots for tRNAs downstream of DESeq2
1. `deseq2`: Run DESeq2 analysis on single copy genes read counts
1. `deseq2_starTE_random`: Run DESeq2 analysis on TE read counts from STAR-TE random quantification
1. `deseq2_tRNA`: Run DESeq2 analysis on tRNA read counts
1. `download_genome_annotation_file`: Download GTF annotations
1. `download_genome_fasta_file`: Download genome Fasta sequence
1. `download_gtRNAdb`: Download gtRNAdb annotation
1. `download_repeatmasker_annotation_file`: Download RepeatMasker and produces a GTF out of it
1. `edit_condition_file`: prepares a so-called condition file for SalmonTE analysis
1. `fastqc_markdup`: Run FastQC on deduplicated alignments
1. `fastqc_raw`: Run FastQC on raw reads
1. `fastqc_star`: Run FastQC on BAM files downstream of STAR alignment
1. `fastqc_trim`: Run FastQC on trimmed reads (single-end libraries)
1. `fastqc_trim_pe`: Run FastQC on trimmed reads (paired-end libraries)
1. `featureCounts_multihit`: Run featureCounts downstream of STAR-TE multihit
1. `featureCounts_random`: Run featureCounts downstream of STAR-TE random
1. `filter_bam`
1. `get_runinfo`
1. `index_bam`
1. `mk_genome_tsv`
1. `multiqc_markdup`
1. `multiqc_raw`: Aggregates FastQC results of raw reads
1. `multiqc_star`: Aggregates FastQC results of STAR alignment (single copy genes)
1. `multiqc_trim`: Aggregates FastQC results of trimming procedure
1. `parse_runinfo`
1. `picard_markdup`: Marks duplicated reads from alignment using Picard MarkDuplicates
1. `salmonTE_quant`: Run SalmonTE quantification routing
1. `salmonTE_test`: Perform differential expression testing from SalmonTE quantifications
1. `star`: Run STAR alignment (single copy genes)
1. `star_genome_preparation`: Generates genome indexes
1. `starTE_multihit`: Run STAR-TE multihit against RepeatMasker annotation
1. `starTE_random`: Run STAR-TE random against RepeatMasker annotation
1. `subset_gtf`
1. `trimmomatic_pe`: Perform reads trimming with Trimmomatic (paired-end libraries)
1. `trimmomatic_se`: Perform reads trimming with Trimmomatic (single-end libraries)
1. `validate_genome_and_annotation`: Check that chromosome names match between Fasta sequence and GTF annotation
1. `verify_star`: Validate STAR outputs
1. `yte_single_copy_genes`: Prepare YAML template for Datavzrd - single copy genes
1. `yte_starTE_random`: Prepare YAML template for Datavzrd - TE from STAR-TE random
1. `yte_trna`: Prepare YAML template for Datavzrd - tRNA
