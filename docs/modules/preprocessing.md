# Preprocessing

The preprocessing module ensures that sequencing data is of high quality and free from technical artifacts before alignment.

## Workflow

```mermaid
graph LR
    A[Raw FASTQ] --> B[FastQC Raw]
    B --> C[Adaptive Trimming]
    C --> D[FastQC Trimmed]
    D --> E[Processed FASTQ]
```

## Quality Control (FastQC)

The pipeline automatically executes [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on every raw fastq file.

- Results are stored in `qc/fastqc-raw/`.
- Key metrics (Sequence Quality, Adapter Content) are used to drive the downstream trimming process.

## Trimming (Trimmomatic)

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is used to remove adapters and low-quality bases.

### Adaptive Trimming

By default, 3t-seq uses its [Adaptive Trimming](../adaptive-trimming.md) algorithm. This dynamically selects the best adapter file and quality thresholds based on the specific library's properties.

### Manual Configuration

You can provide fixed parameters for specific libraries in the `config.yaml`:

```yaml
sequencing_libraries:
  - name: Sample1
    trimmomatic:
      adaptive: false
      extra_params: "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
```

## Parameters & Defaults

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `trimmomatic.adaptive` | `false` | Enable/disable [Adaptive Trimming](../adaptive-trimming.md). |
| `trimmomatic.extra_params` | - | Fixed Trimmomatic modules (e.g., `CROP:100`). |

---

## Results

| Location | Description |
| :--- | :--- |
| `results/trim/` | Contains the trimmed FASTQ files. |
| `results/qc/fastqc-raw/` | FastQC reports for the raw input data. |
| `results/qc/fastqc-trimmed/` | FastQC reports for the data after trimming. |
| `results/trim/_shared/` | Deduplicated shared trimmed files (internal). |
