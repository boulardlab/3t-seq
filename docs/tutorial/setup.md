# Preparing Data & Samples

To run the pipeline on your own data, you need to organize your raw files and create a **Sample Sheet** CSV. This file maps your raw sequencing files to sample names and biological conditions.

---

## 1. Organizing Your Raw Files

The pipeline can work with FASTQ files anywhere on your system, but we recommend placing them in a `data/` or `reads/` directory within your project root.

- **Fastq files** should be gzipped (`.fastq.gz` or `.fq.gz`).
- **Paired-end** files typically use `_R1` / `_R2` or `_1` / `_2` suffixes.

---

## 2. Creating a Sample Sheet

Each sequencing library requires a corresponding CSV sample sheet. This file is the primary way the pipeline discovers your data.

### Column Reference

The following columns are grounded in the pipeline's internal logic:

| Column | Required | Description |
| :--- | :--- | :--- |
| `name` | **Yes** | A unique biological identifier for the sample (e.g., `WT_Rep1`). |
| `filename` | **SE** | The path to the raw sequence file for **Single-End** reads. |
| `filename_1` | **PE** | The path to the first mate for **Paired-End** reads. |
| `filename_2` | **PE** | The path to the second mate for **Paired-End** reads. |
| `condition` | No* | Metadata for **DESeq2** (e.g., `WT`, `KO`). You can add any number of additional metadata columns. |

!!! note
    *While additional columns like `condition` are not strictly required for alignment, they are essential for downstream Differential Expression (DE) analysis.

### Examples

#### Single-End (SE) Sample Sheet

```csv
name,filename,condition
Control_Rep1,reads/ctrl1.fq.gz,control
Control_Rep2,reads/ctrl2.fq.gz,control
Treated_Rep1,reads/treat1.fq.gz,treated
```

#### Paired-End (PE) Sample Sheet

```csv
name,filename_1,filename_2,treatment
Sample_A,data/sA_R1.fastq.gz,data/sA_R2.fastq.gz,Basal
Sample_B,data/sB_R1.fastq.gz,data/sB_R2.fastq.gz,Basal
Sample_C,data/sC_R1.fastq.gz,data/sC_R2.fastq.gz,TGFb
```

---

## 3. Library Preparation Reference

In your configuration, you will link these sample sheets to a library name and specify the protocol (`se` or `pe`).

```yaml
sequencing_libraries:
  - name: "MyExperiment"
    protocol: "pe" # Must match your sample sheet columns!
    sample_sheet: "my_samples.csv"
```

!!! warning
    If you specify `protocol: "pe"`, the pipeline expects `filename_1` and `filename_2` in the sample sheet. If you specify `protocol: "se"`, it expects `filename`.

---

## Next Steps

Once your data is prepared, learn how to use **Profiles** to encapsulate all your run settings:

- [**Advanced Profiles**](profiles.md): Logic of the `--profile` flag and encapsulation.
