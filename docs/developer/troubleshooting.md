# Troubleshooting

Each entry gives the symptom, the cause and the action.
For errors in pipeline runs on user data, read [Usage: Troubleshooting](../usage/index.md#troubleshooting).

## Pixi: "received more arguments than expected"

**Symptom.** A pixi task stops before Snakemake starts:

```text
× task 'test' received more arguments than expected
help: use `--` to separate task arguments from extra passthrough arguments
```

**Cause.** The tasks `lint`, `dry-run`, `test` and `test-report` accept only two positional arguments, `config` and `profile`.
Pixi does not give other arguments to Snakemake.

**Action.** Write `--` before the Snakemake arguments:

```bash
pixi run -e dev test adaptive-trim laptop -- --cores 4
```

## Pixi: "could not find the task 'setup'"

**Symptom.** `pixi run test` stops with `could not find the task 'setup'`.

**Cause.** `test` depends on `setup`. `setup` exists only in the `dev` and `docs` environments.

**Action.** Run the task in `dev`: `pixi run -e dev test <config> <profile>`.

## Apptainer: "Could not write info to setgroups"

**Symptom.** Each container job fails. The log contains:

```text
Could not write info to setgroups: Permission denied
```

**Cause.** Ubuntu 24.04 and later use AppArmor to block unprivileged user namespaces.
Rootless Apptainer needs them.

**Action.** Allow unprivileged user namespaces:

```bash
sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0
```

The setting stays until the next reboot. The `e2e` job in `test.yml` does this step on the CI runner.

## Apptainer: job hangs in a container

**Symptom.** A container job starts but does not end, and its log stops changing. Trimmomatic is an example.
In CI, the `e2e` job stops at its 90-minute limit. The job times show the rule as "unfinished".

**Cause.** The Apptainer package from conda-forge uses `squashfuse_ll`.
That `squashfuse_ll` can hang when it reads files in a SIF image.
See [conda-forge/squashfuse-feedstock#21](https://github.com/conda-forge/squashfuse-feedstock/issues/21) and issue #39.
Older 3t-seq versions installed this Apptainer through pixi on `linux-64`.

**Action.**

1. Run `which apptainer`. If the path is in `.pixi/` or in a conda env, do not use that Apptainer.
2. Run `pixi install` to update the pixi environments. The current `pixi.toml` does not install Apptainer.
3. Install Apptainer 1.2.5 or later from your system packages or your cluster modules.
   See [Getting Started](../getting-started/index.md#apptainer).
4. Run `which apptainer` again. Make sure that it shows the system Apptainer.

## Rule: "TMPDIR: unbound variable"

**Symptom.** A FastQC job fails in the container. The log contains `TMPDIR: unbound variable`.

**Cause.** Snakemake runs shell commands with `set -u`. The container did not get `$TMPDIR`, so the variable is not set.

**Action.** Use `{resources.tmpdir}` in rule shell commands, not `$TMPDIR`. Snakemake writes the value into the command.
The FastQC rules use `{resources.tmpdir}` since PR #27. Snakemake always gives `tmpdir` a value.
The `laptop` test profile sets it to `/tmp` in `default-resources`.

## Container check waits or fails on missing images

**Symptom.** The `check` job prints `Waiting for: ghcr.io/boulardlab/...` for a long time.
Or it fails with `missing images in ghcr.io`.

**Cause.** One of these:

- A build workflow still runs for the PR head commit. The job waits for it, up to 210 minutes.
- A build workflow failed. The job prints `a build workflow failed for <sha>`.
- No workflow builds the image. The job prints `No build running for <sha>`. Examples:
    - The PR comes from a fork.
    - The pin changed, but the `Dockerfile` did not change, so `build-3t-seq.yml` did not start.
    - The `salmonte` pin has a version that is not built.

**Action.**

1. Open the `build-3t-seq.yml` or `build-salmonte.yml` run for the commit. Read its log.
2. If the build failed, correct the cause and push again.
3. If no build ran, start one with `gh workflow run build-3t-seq.yml --ref <branch>`.
4. When the image exists, run the failed `check` job again.

For fork PRs, read [Fork PRs](branches.md#fork-prs).

## `dockerfile-sync` shows a diff

**Symptom.** The job `dockerfile-sync` fails with:

```text
Dockerfile or its pin is stale. Run `pixi run containerize` and commit both files.
```

**Cause.** A conda env in `workflow/envs/` changed. The committed `Dockerfile` or the pin in `workflow/Snakefile` is old.
Another cause is a manual edit of one of these files.

**Action.**

1. Run `pixi run containerize`.
2. Commit `Dockerfile` and `workflow/Snakefile` together.
3. Push. Wait for `build-3t-seq.yml` to build the new `env-<hash>`.

## Snakemake: "Unable to find environment in container image"

**Symptom.** A job stops with `Unable to find environment in container image`.

**Cause.** The pinned image does not contain the conda env of the rule. The pin and the conda envs do not agree.

**Action.** Do the procedure in [Change a conda env](containers.md#procedure-change-a-conda-env).

## release-please: "Bad credentials"

**Symptom.** The `release-please` job on `main` fails with `Bad credentials`.

**Cause.** `RELEASE_PLEASE_TOKEN` expired, or its owner deleted it.

**Action.** Replace the token. See [Replace the token](release.md#replace-the-token).

## Release PR runs no checks

**Symptom.** The release PR or the git tag `vX.Y.Z` starts no workflows.
The required checks stay at "Expected — Waiting for status to be reported".

**Cause.** Release-please used `GITHUB_TOKEN`. GitHub does not start workflows for events from `GITHUB_TOKEN`.
`release-please.yml` uses `GITHUB_TOKEN` when the secret `RELEASE_PLEASE_TOKEN` is empty or missing.

**Action.**

1. Make sure that the secret exists: `gh secret list --repo boulardlab/3t-seq`.
2. If it does not exist, set it. See [Replace the token](release.md#replace-the-token).
3. Push a commit to the release PR branch, for example an edit of `CHANGELOG.md`.
   A push by a person starts the checks.
4. If the git tag exists but the image workflows did not run, start them on the tag:

    ```bash
    gh workflow run build-3t-seq.yml --ref vX.Y.Z
    gh workflow run build-salmonte.yml --ref vX.Y.Z
    ```
