# Adaptive Trimming

The 3t-seq pipeline implements an **Adaptive Trimming** mechanism that automatically optimizes Trimmomatic parameters based on the specific quality profile of your sequencing data.

## Motivation

Standard fixed-threshold trimming can sometimes be too aggressive (losing valuable data) or too lenient (leaving artifacts that impact alignment). Adaptive trimming bridges this gap by analyzing [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) results before executing the trim step.

## How it Works

The adaptive trimming process consists of three main phases:

### 1. Quality Profiling

The pipeline first runs a "raw" FastQC pass on the input FASTQ files. The results are parsed to extract:

- **Adapter Content**: Detection of Illumina Universal, Nextera, or Small RNA adapters.
- **Per-base Quality**: Specifically targeting the \(10^{th}\) percentile (\(p_{10}\)) quality score. In simpler terms, the pipeline looks at the "worst 10%" of your sequence data to gauge how poor the quality gets, which helps it decide how aggressive the trimming needs to be.
- **Read Length**: The maximum observed read length in the library.

### 2. Parameter Derivation

Using the extracted statistics, the pipeline dynamically assembles a Trimmomatic command:

| Feature | Logic | Default / Derived |
| :--- | :--- | :--- |
| **Adapter Clipping** | Detects the most prevalent adapter type to ensure all synthetic sequences are removed. | `TruSeq3` or `NexteraPE` |
| **Quality Trimming** | Uses standard window trimming by default, but switches to a more aggressive, information-preserving method (`MAXINFO`) if the worst 10% of data has a quality score below 15. | Switches to `MAXINFO:40:0.5` if average \(p_{10} < 15\). |
| **Minimum Length** | Calculates a safe cutoff for downstream alignment. | \(\max(36, 0.35 \times \text{max_length})\) |
| **Base Trimming** | Ensures low-quality ends are removed. | `LEADING:3` and `TRAILING:3` |

### 3. Execution

The derived parameters are saved as a JSON metadata file and passed to Trimmomatic. This ensure that each sample in your experiment receives the optimal level of processing.

---

## User Overrides

While the system is designed to be fully automated, you can still provide fixed parameters in your `config.yaml`:

```yaml
sequencing_libraries:
  - name: MySample
    trimmomatic:
      adaptive: true
      extra_params: "MINLEN:50 CROP:100"
```

!!! note "Precedence"
    User-provided modules in `extra_params` always take precedence over derived ones. For example, if you provide `MINLEN:50`, the adaptive `MINLEN` calculation will be skipped for that sample.

## When to Use It

Adaptive trimming is enabled by default. It is particularly recommended for:

- Samples with variable quality across batches.
- Pilot studies where adapter content or quality decay is unknown.
- Large-scale integrative analyses of public datasets from different sources.
