# CI reference

This page describes each workflow in `.github/workflows/`.
For each workflow, it gives the trigger, the jobs, the failures and the fixes.

## Events and workflows

```mermaid
flowchart LR
    PRD["PR to dev"] --> L["lint.yml"]
    PRD --> T["test.yml"]
    PRD --> C["container-check.yml"]
    PRD -.->|Dockerfile changed| B3["build-3t-seq.yml"]
    PRM["PR to main"] --> L
    PRM --> T
    PRM --> C
    PRM -.->|Dockerfile changed| B3
    PRM -.->|release PR only| BS["build-salmonte.yml"]
    PD["push to dev"] -.->|path filter| B3
    PD -.->|path filter| BS
    PD -.->|path filter| D["docs.yml"]
    PM["push to main"] --> L
    PM --> T
    PM --> C
    PM --> RP["release-please.yml"]
    PM -.->|path filter| B3
    PM -.->|path filter| BS
    PM -.->|path filter| D
    TAG["push of tag v*"] --> B3
    TAG --> BS
    N["nightly 02:17 UTC"] --> T
```

A dashed arrow means that the workflow runs only when its condition is true.

| Workflow | PR to `dev` | PR to `main` | Push to `dev` | Push to `main` | Tag `v*` | Nightly | Manual |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `lint.yml` | Yes | Yes | No | Yes | No | No | No |
| `test.yml`, job `smoke` | Yes | Yes | No | Yes | No | No | No |
| `test.yml`, jobs `images` and `e2e` | No | Release PR only | No | No | No | Yes | Yes |
| `container-check.yml` | Yes | Yes | No | Yes | No | No | Yes |
| `build-3t-seq.yml` | Path filter | Path filter | Path filter | Path filter | Yes | No | Yes |
| `build-salmonte.yml` | No | Release PR only | Path filter | Path filter | Yes | No | Yes |
| `release-please.yml` | No | No | No | Yes | No | No | No |
| `docs.yml` | No | No | Path filter | Path filter | No | No | Yes |

`test.yml` also calls `container-check.yml` as a reusable workflow (`workflow_call`), in the job `images`.

## Required checks

The ruleset for `main` requires three checks. See [Ruleset for `main`](branches.md#ruleset-for-main).

| Check | Workflow |
| :--- | :--- |
| `check` | `container-check.yml` |
| `dockerfile-sync` | `container-check.yml` |
| `smoke` | `test.yml` |

The job `e2e` is not a required check. Do not merge a release PR when `e2e` fails.

When `test.yml` calls `container-check.yml`, GitHub shows the jobs as `images / check` and `images / dockerfile-sync`.
These are the same jobs. The ruleset uses the jobs from the direct `container-check.yml` run.

## `lint.yml` (Format and lint)

| | |
| :--- | :--- |
| Trigger | PR to `main` or `dev`; push to `main` |
| Job | `lint` |
| Steps | `pixi run -e dev format-check`, then `pixi run -e dev lint` |

**Failures:**

- *`format-check` prints a diff.* `snakefmt` wants to change the formatting of a file in `workflow/`.
  Run `pixi run -e dev format`. Commit the changed files.
- *`lint` reports warnings.* `snakemake --lint` found a problem in a rule.
  Run `pixi run -e dev lint` on your computer. Correct each warning.

## `test.yml` (Run tests)

| | |
| :--- | :--- |
| Trigger | PR to `main` or `dev`; push to `main`; nightly at 02:17 UTC; manual |
| Jobs | `smoke`; `images` (calls `container-check.yml`), then `e2e` |
| Profile | `ci` (`.tests/integration/profiles/ci/`): 4 cores, 2 jobs, for the GitHub runner |

### Job `smoke`

`smoke` is fast (about 30 seconds) and does not run the pipeline. It runs on PRs and on pushes to `main`.

1. It runs `pixi run -e dev lint <config> laptop` for the three configs.
2. It runs the unit tests: `pixi run -e dev unit-test`.
3. It runs `pixi run -e dev dry-run <config> ci -- --sdm conda` for the three configs.

The dry run uses conda only, because a dry run with Apptainer needs the pulled images.

**Failures:**

- *`lint` reports warnings.* Run `pixi run -e dev lint <config> laptop`. Correct each warning.
- *A unit test fails.* Run `pixi run -e dev unit-test` on your computer. Correct the helper or the test.
  See [Unit tests](local-checks.md#unit-tests).
- *The dry run fails.* Most often, the config is not valid or a rule input is missing.
  Run `pixi run -e dev dry-run <config> ci -- --sdm conda` on your computer. Read the error.
- *Pixi stops before Snakemake starts.* Read [Troubleshooting](troubleshooting.md).

### Job `e2e`

`e2e` runs the full pipeline on the `local-references` test data:

```bash
pixi run -e dev test local-references ci -- --singularity-args="--bind $ROOT --bind $HOME"
```

It runs on release PRs, every night and on a manual start. It starts only after `images` passes.
On a release PR, the new SalmonTE image can still be in the build.

The job also does these steps:

- It checks out the Git LFS test data.
- It runs on `ubuntu-24.04`, not `ubuntu-latest`. The PPA has no Ubuntu 26.04 build yet. On 26.04, apt installs Ubuntu's own Apptainer 1.4.5 instead. That build extracts the image to disk and fails with "Disk quota exceeded".
- It installs Apptainer from `ppa:apptainer/ppa`, not from conda-forge. See [Troubleshooting](troubleshooting.md#apptainer-job-hangs-in-a-container).
- It sets `kernel.apparmor_restrict_unprivileged_userns=0`. Rootless Apptainer needs this on Ubuntu 24.04.
- It keeps the pulled images in the GitHub Actions cache. The cache key is a hash of the pins.
- The job stops after 90 minutes.

To start `e2e` on a branch, run:

```bash
gh workflow run test.yml --ref <branch>
```

A manual start runs `images` and `e2e`, but not `smoke`.

### Read the job times

The step "Job times" runs `.github/scripts/snakemake_job_times.py` on the Snakemake log.
It writes a table to the job summary, also when the job fails.

1. Open the `test.yml` run: `gh run list --workflow test.yml`, then `gh run view <run-id> --web`.
2. Scroll down to the summary of the `e2e` job.
3. Read the table. Each row gives a rule, its number of jobs, the total minutes and the maximum minutes.
   The rule with the highest total is first. A rule marked "unfinished" did not end.
4. Read the "Snakemake wall time" line under the table.

To make the same table from a local run, run:

```bash
python3 .github/scripts/snakemake_job_times.py .tests/integration/.snakemake/log/<log-file>
```

**Failures:**

- *`images` fails.* Read [`container-check.yml`](#container-checkyml-container-check).
- *A rule fails.* The `ci` profile shows the log of the failed rule. The run also uploads the artifact `e2e-logs`.
  Download it with `gh run download <run-id> -n e2e-logs`.
- *The job stops at 90 minutes.* Read the job times. Find the rule that did not end.
  See [Troubleshooting](troubleshooting.md#apptainer-job-hangs-in-a-container).
- *To reproduce a failure*, run `pixi run -e dev test local-references laptop` on your computer.
  See [Local checks](local-checks.md).

## `container-check.yml` (Container check)

| | |
| :--- | :--- |
| Trigger | PR to `main` or `dev`; push to `main`; manual; `workflow_call` |
| Jobs | `dockerfile-sync` and `check`, in parallel |

### Job `dockerfile-sync`

1. It runs `pixi run containerize`.
2. It runs `git diff --exit-code -- Dockerfile workflow/Snakefile`.

**Failure:** *`git diff` shows a difference.*
A conda env changed, but the `Dockerfile` or the pin is old.
Run `pixi run containerize`. Commit `Dockerfile` and `workflow/Snakefile` together.

### Job `check`

1. It runs `container_pin.py check`. The pin in `workflow/Snakefile` must equal `env-<hash>` from the `Dockerfile`.
2. It runs `check-containers.sh`. All pinned images must have a fixed tag and must exist in `ghcr.io`.

On a PR, `check-containers.sh` waits for the build workflows of the PR head commit.
It checks the images every 60 seconds, for a maximum of 210 minutes.
It checks the build status only after 3 minutes.
After that time, it stops when a build fails or when no build runs.
It waits only on a PR. This includes the call from `test.yml` job `images` for a release PR.
On a push to `main`, a manual run or the nightly run, it checks one time and does not wait.

**Failures:**

| Message | Cause | Action |
| :--- | :--- | :--- |
| `pins ... but the Dockerfile needs env-...` | The pin and the `Dockerfile` do not agree. | Run `pixi run containerize`. Commit both files. |
| `mutable image refs` | A pin uses `:latest` or has no tag. | Pin `env-<hash>`, `X.Y.Z` or a digest. |
| `a build workflow failed for <sha>` | `build-3t-seq.yml` or `build-salmonte.yml` failed. | Correct the build. Run the failed jobs again. |
| `No build running for <sha>` and `missing images` | No workflow builds the pinned image. Example: a fork PR. | Build the image manually. See [Fork PRs](branches.md#fork-prs). |

## `build-3t-seq.yml` (Build 3t-seq)

| | |
| :--- | :--- |
| Trigger: PR | PR to `main` or `dev`, when `Dockerfile` or the workflow file changes |
| Trigger: push | Push to `main` or `dev`, with the same path filter; push of a `v*` tag |
| Trigger: manual | `workflow_dispatch` |
| Job | `docker` (skipped for fork PRs) |
| Platforms | `linux/amd64`, `linux/arm64` |

1. It reads the tag `env-<hash>` with `container_pin.py tag`.
2. If `env-<hash>` does not exist, it builds and pushes the image.
3. If `env-<hash>` exists, it adds the other tags to the existing image. It does not build.
4. On a `v*` tag, it adds the tag `X.Y.Z`, but only if `X.Y.Z` does not exist.

The [Containers](containers.md) page lists all tags.

**Failures:**

- *A conda env does not solve.* Most often, a package is missing for `linux/arm64`.
  Select a version that exists for both platforms. Run `pixi run containerize` again.
- *The build is slow.* The `arm64` build uses QEMU emulation. Wait for it to finish.

## `build-salmonte.yml` (Build SalmonTE)

| | |
| :--- | :--- |
| Trigger: PR | Release PR to `main` |
| Trigger: push | Push to `main` or `dev`, when `assets/docker/salmonTE/Dockerfile` or the workflow file changes; push of a `v*` tag |
| Trigger: manual | `workflow_dispatch`, with an optional `version` input |
| Job | `docker` (on PRs, only for `release-please--*` branches in this repository) |
| Platforms | `linux/amd64`, `linux/arm64` |

The workflow gets the version from these sources:

| Event | Version |
| :--- | :--- |
| Release PR | `.release-please-manifest.json` at the PR head |
| Manual | The `version` input |
| Tag `vX.Y.Z` | `X.Y.Z` |
| Push to a branch | None. It builds only branch tags. |

**Failures:**

| Message | Cause | Action |
| :--- | :--- | :--- |
| `invalid version` | The version is not `X.Y.Z`. | Give a version like `2.0.1`. |
| `vX.Y.Z is released and ... exists; refusing to overwrite` | A manual run or PR tried to replace a released image. | Do not replace released images. Release a new version. |
| `vX.Y.Z is already tagged; the release PR should propose a new version` | The release PR has the version of an existing git tag. | Correct `.release-please-manifest.json` or the `Release-As` footer. |

## `release-please.yml` (release-please)

| | |
| :--- | :--- |
| Trigger | Push to `main` |
| Job | `release-please` |
| Token | `RELEASE_PLEASE_TOKEN`; if the secret is empty, `GITHUB_TOKEN` |

- When `main` has new releasable commits, it opens or updates the release PR.
- When the release PR merges, it creates the git tag `vX.Y.Z` and a GitHub release.

Configuration: `release-please-config.json` and `.release-please-manifest.json`.
The `extra-files` are `pixi.toml` and `workflow/rules/salmonTE.smk`.

**Failures:**

- *`Bad credentials`.* The token expired, or its owner deleted it. See [Troubleshooting](troubleshooting.md#release-please-bad-credentials).
- *The release PR runs no checks.* The workflow used `GITHUB_TOKEN`.
  See [Troubleshooting](troubleshooting.md#release-pr-runs-no-checks).

The [Release procedure](release.md) gives the full steps.

## `docs.yml` (Docs)

| | |
| :--- | :--- |
| Trigger | Push to `main` or `dev` when `docs/`, `mkdocs.yml`, `pixi.toml` or the workflow file changes; manual |
| Jobs | `build` (`mkdocs build --clean --strict`), then `deploy` to GitHub Pages |

**Failure:** *`mkdocs build --strict` stops on a warning.*
Most often, a link is broken or a page is not in the `nav` of `mkdocs.yml`.
Run `pixi run -e docs docs-build` on your computer. Correct each warning.

!!! warning
    Pushes to `dev` and to `main` both deploy to the same GitHub Pages site.
    The site shows the last branch that CI deployed.

## Scripts

### `.github/scripts/container_pin.py`

Keeps the `containerized:` pin in `workflow/Snakefile` equal to the env hash in the `Dockerfile`.

| Command | Result |
| :--- | :--- |
| `python3 .github/scripts/container_pin.py tag` | Prints `env-<hash>` for the current `Dockerfile`. |
| `python3 .github/scripts/container_pin.py update` | Writes that tag into `workflow/Snakefile`. `pixi run containerize` runs this command. |
| `python3 .github/scripts/container_pin.py check` | Exits with 1 if the pin and the `Dockerfile` do not agree. |

### `.github/scripts/check-containers.sh`

Finds every `ghcr.io/boulardlab/*` image reference in `workflow/`. It skips comment lines.

| Option | Result |
| :--- | :--- |
| `--list` | Prints the references only. It does not access the registry. |
| `--wait-minutes N` | Waits up to `N` minutes for missing images. The default is 0 (check one time). |
| `--head-sha SHA` | Follows the build workflow runs for commit `SHA`. Needs `gh` and `GH_TOKEN`. |

The environment variable `GRACE_SECONDS` sets the time before the script looks for build runs. The default is 180.
