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
- **Read Length**: The maximum read length, from the `Sequence length` line of FastQC's
  Basic Statistics (50 if it cannot be read).

### 2. Parameter Derivation

Using the extracted statistics, the pipeline dynamically assembles a Trimmomatic command:

| Feature | Logic | Default / Derived |
| :--- | :--- | :--- |
| **Adapter Clipping** | Takes the adapter with the highest Adapter Content peak across all mates. If it is Nextera and above 0.1 %, uses Nextera adapters; in every other case (Illumina Universal, Small RNA, or no adapter above 0.1 %) uses TruSeq3. | `ILLUMINACLIP:<file>:2:30:10` with `TruSeq3-PE.fa` / `TruSeq3-SE.fa`, or `NexteraPE-PE.fa` (also for single-end) |
| **Quality Trimming** | Averages the \(p_{10}\) of every base position over all mates (20 if FastQC has no quality data). Uses window trimming, or the information-preserving `MAXINFO` when the average is below 15. | `SLIDINGWINDOW:4:20` (the non-adaptive default is `4:15`); `MAXINFO:40:0.5` if average \(p_{10} < 15\). |
| **Minimum Length** | Calculates a safe cutoff for downstream alignment. | \(\max(36, 0.35 \times \text{max_length})\) |
| **Base Trimming** | Ensures low-quality ends are removed. | `LEADING:3` and `TRAILING:3` |

### 3. Execution

The derived parameters are saved per sample to
`<results_folder>/trim/_shared/<trim_hash>/<sample>.adaptive_params.json` and passed to
Trimmomatic. This ensures that each sample in your experiment receives its own level of
processing.

---

## User Overrides

While the system is designed to be fully automated, you can still provide fixed parameters in your `config.yaml`:

```yaml
comparisons:
  - name: MySample
    trimmomatic:
      adaptive: true
      extra_params: "MINLEN:50 CROP:100"
```

!!! note "Precedence"
    A step in `extra_params` replaces derived steps of the same type. For example, if you provide `MINLEN:50`, the adaptive `MINLEN` calculation is skipped for that sample. `ILLUMINACLIP` always runs first and `MINLEN` last; if you set `CROP:N` but no `MINLEN`, the derived `MINLEN` is lowered to N. See [how `extra_params` is merged](../modules/preprocessing.md#how-extra_params-is-merged).

## When to Use It

Adaptive trimming is off by default; enable it with `trimmomatic: {adaptive: true}`, or
the shorthand `trimmomatic: adaptive`. It is particularly recommended for:

- Samples with variable quality across batches.
- Pilot studies where adapter content or quality decay is unknown.
- Large-scale integrative analyses of public datasets from different sources.
