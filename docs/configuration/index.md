# Configuration Reference

This page describes all parameters in the `config.yaml` file. The pipeline uses [Snakemake's schema validation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#validation) to ensure your configuration is correct before starting.

---

## 1. Core Configuration

The most commonly modified settings for any 3t-seq run. Use this annotated YAML as a guide for your basic setup.

```yaml
globals:
  results_folder: "results/my_analysis/" # (1)

genome:
  label: "mm10" # (2)

sequencing_libraries:
  - name: "GSE123456" # (3)
    protocol: "pe" # (4)
    sample_sheet: "samples.csv" # (5)

defaults: # (6)
  strandedness: 0
```

1. **results_folder**: **Required**. Root directory where all BAMs, counts, and logs will be saved. Always use a unique name for each experiment.
2. **label**: **Required**. Genome version (e.g., `mm10`, `mm39`). This triggers automated downloads of all necessary references.
3. **name**: **Required**. Unique identifier for the library or series (e.g., GEO accession).
4. **protocol**: Sequencing geometry (`se` for single-end, `pe` for paired-end).
5. **sample_sheet**: Path to the CSV mapping raw files to biological sample names.
6. **defaults**: Global parameters that apply to all libraries unless overridden individually.

---

## 2. Library & Alignment Settings

Use these tabs to configure specialized alignment logic and library-specific overrides.

=== "Reference genome"
    Advanced control for custom genome and which chromosomes to process.

    ```yaml
    genome:
      label: "mm10"
      fasta_path: "refs/genome.fa"
      selected_chromosomes: ["chr1", "chr2"]
    ```

    `fasta_path`
    : **Type**: `string` | **Default**: Derived from `label`
    : Local path to a custom genome FASTA file. Required if using a non-standard reference.

    `gtf_path`
    : **Type**: `string` | **Default**: Derived from `label`
    : Local path to a custom genome GTF annotation.

    `annotation_type`
    : **Type**: `string` | **Default**: `ensembl`
    : Format of your GTF (`ensembl`, `gencode`, or `mgi`).

    `selected_chromosomes`
    : **Type**: `array` | **Default**: `null` (All)
    : Optional list of chromosomes (e.g., `["chr1", "chr2"]`). Useful for focusing analysis or speeding up pilot runs.

=== "Sequencing libraries"
    Custom parameters for each library in your experiment. You can customize trimming, alignment, bigwig, and other parameters for each library.

    ```yaml
    sequencing_libraries:
      - name: "GSE123"
        protocol: "pe"
        trimmomatic:
          adaptive: true
        star: "--seedSearchStartLmax 30"
    ```

    `trimmomatic`
    : **Type**: `object/string` | **Default**: Standard flags
    : Custom trimming parameters. Set to `adaptive: true` to enable automated parameter derivation based on FastQC results.

    `star`
    : **Type**: `string` | **Default**: `""`
    : Extra CLI flags to pass directly to the STAR aligner for primary genome mapping.

    `bamCoverage`
    : **Type**: `string` | **Default**: `""`
    : Custom flags for `deeptools bamCoverage` (e.g., `--normalizeUsing CPM`).

=== "Global defaults"
    Shared settings for statistical comparisons.

    ```yaml
    defaults:
      strandedness: 2
      deseq2:
        test: "Wald"
        variable: "condition"
    ```

    `deseq2`
    : **Type**: `object` | **Required**: for DE
    : Configuration for group comparisons (Wald/LRT tests and reference levels).

    `strandedness`
    : **Type**: `integer` | **Default**: `0` (Unstranded)
    : Library preparation geometry. **`1`**: Forward, **`2`**: Reversely stranded.

---

## 3. Expert Parameters & Module Flags

Advanced settings for internal modules. These are hidden by default to prioritize scannability.

!!! abstract "Specialized Analysis Modules (Flags)"
    ```yaml
    disable_TE_analysis: true
    disable_tRNA_analysis: false
    ```
    Set these to `true` to skip specific parts of the pipeline and reduce computation time.

    - `disable_TE_analysis`: Skips STAR-TE and SalmonTE quantification.
    - `disable_salmonTE_analysis`: Skips secondary SalmonTE processing.
    - `disable_tRNA_analysis`: Skips specialized tRNA mapping.

???+ info "STAR-TE Internal Modes"
    ```yaml
    defaults:
      starTE_random:
        outFilterMultimapNmax: 5000
    ```
    Fine-grained control over multi-mapping read assignment (multi-hits) in TEs.

    - `starTE_random`: Settings for multimap assignment (e.g., `outFilterMultimapNmax: 5000`).
    - `starTE_multihit`: Settings for fractional counting mode.

???+ info "tRNA Expert Settings"
    ```yaml
    tRNA_quantification:
      method: "mim-tRNA-seq"
      mimseq_params:
        max_mismatches: 2
    ```
    Configuration for specialized tRNA sequencing protocols.

    - `method`: `standard` or `mim-tRNA-seq` (for clinical-grade tRNA kits).
    - `mimseq_params`: Sub-parameters like `max_mismatches` and `min_cov` for the mim-tRNA-seq logic.
