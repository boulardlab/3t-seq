# Advanced Profiles

In **3t-seq**, we separate your **biological setup** (what files to process, the
`--configfile`) from your **computational strategy** (how to run them, the `--profile`).

---

## 1. What a profile holds

A Snakemake profile is a folder with a `config.yaml` whose keys are Snakemake command-line
options written in YAML (`jobs: 100` is `--jobs 100`). A profile typically sets:

- **Execution**: run jobs locally or submit them through Slurm (`executor: slurm`).
- **Concurrency**: how many cores and how many jobs at the same time.
- **Software deployment**: Conda environments and the Apptainer container
  (`software-deployment-method: [conda, apptainer]`), plus extra container options.
- **Resources**: default and per-rule resources, such as the Slurm account and partition.

Rules in 3t-seq set their own `threads`, `mem_mb` and `runtime`, scaled to the input size,
so a profile only needs the site-specific settings.

---

## 2. Profiles shipped with the repository

| Profile | Use it for |
| :--- | :--- |
| `profile/default` | Local runs on a workstation or a single node. |
| `profile/slurm-example` | Template for a Slurm cluster; copy it and fill in account and partition. |

### `profile/default`: local runs

```yaml
cores: 2
local-cores: 2
jobs: 1
restart-times: 0
printshellcmds: True
software-deployment-method:
  - conda
  - apptainer
apptainer-args: "--no-home"
```

It runs one job at a time on 2 cores. Raise `cores` (or pass `--cores N`) on a bigger
machine.

```bash
pixi run snakemake --profile profile/default --configfile my_config.yaml
```

### `profile/slurm-example`: Slurm clusters

The example uses the Slurm executor plugin (`snakemake-executor-plugin-slurm`, installed
with the Pixi environment). Copy it to `profile/slurm/`, which is git-ignored, and replace
the placeholders:

```bash
cp -r profile/slurm-example profile/slurm
```

The settings that matter:

```yaml
executor: slurm          # submit every job through Slurm
jobs: 100                # up to 100 jobs queued or running at once
software-deployment-method:
  - conda
  - apptainer
apptainer-args: "--no-home --bind /scratch"   # bind the file systems holding your data
latency-wait: 120        # wait for files on slow shared file systems
retries: 2               # resubmit failed jobs; heavy rules get more time and memory
default-resources:
  slurm_account: your_account
  slurm_partition: your_partition
```

The file also sets `printshellcmds`, `show-failed-logs`, conservative
`max-jobs-per-second` / `max-status-checks-per-second`, and shows commented
`set-resources` / `set-threads` examples, e.g. to send the STAR rules to a large-memory
partition.

```bash
pixi run snakemake --profile profile/slurm --configfile my_config.yaml
```

Run Snakemake itself on a login node or in a long-running, low-resource job (e.g. inside
`tmux`); it only submits and monitors jobs.

!!! note "Bind mounts"
    Jobs run inside the 3t-seq container. Apptainer mounts the working directory by
    default; any other location your FASTQ files, references or results live in (e.g.
    `/scratch`, a group share) must be added with `--bind` in `apptainer-args`.

Profiles written for 3t-seq 1.x with `slurm: True` no longer work in Snakemake 8; see
[Migrating from 1.x](../migration/index.md#profiles).

---

## 3. Your own profiles

Put any command-line options you would otherwise retype into a profile of your own, e.g.
`profile/my_cluster/config.yaml` based on `profile/slurm-example`:

```yaml
# profile/my_cluster/config.yaml
executor: slurm
jobs: 50
software-deployment-method:
  - conda
  - apptainer
apptainer-args: "--no-home --bind /scratch --bind /g/mylab"
latency-wait: 120
default-resources:
  slurm_account: "my_lab_account"
  slurm_partition: "htc"
```

```bash
pixi run snakemake --profile profile/my_cluster --configfile my_samples.yaml
```

Keep the analysis settings (comparisons, genome, DESeq2) in the `--configfile`; the
profile should only describe how and where jobs run.

---

## 4. Test profiles

`.tests/integration/profiles/` holds the profiles used by the Pixi tasks (`pixi run test`,
`dry-run`, `lint`) on the bundled test dataset. They are passed with `--workflow-profile`
and are not meant for your own data:

| Profile | Settings |
| :--- | :--- |
| `laptop` | 2 cores, 1 job, Conda + Apptainer, binds `$PIXI_PROJECT_ROOT` and `$HOME`. Default for all tasks. |
| `ci` | 4 cores, 2 jobs; used on GitHub Actions runners. |
| `hpc` | 2 cores, Conda + Apptainer; no Slurm settings. |

---

## Next Steps

Now that you understand profiles, learn how to execute the pipeline and generate reports:

- [**Running & Reporting**](running.md): execution patterns and the Snakemake report.
