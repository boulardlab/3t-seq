# Containers

The pipeline uses two images from the GitHub Container Registry (`ghcr.io`):

| Image | Contents | Pin | Pin location |
| :--- | :--- | :--- | :--- |
| `ghcr.io/boulardlab/3t-seq` | All conda envs of the pipeline | `:env-<hash>` | `containerized:` line in `workflow/Snakefile` |
| `ghcr.io/boulardlab/3t-seq/salmonte` | SalmonTE | `:X.Y.Z` (release version) | `container:` lines in `workflow/rules/salmonTE.smk` |

All pins must use a fixed tag. Do not use `:latest` and do not omit the tag.
The `check` job in `container-check.yml` enforces this rule.

## 3t-seq image

### How the pin works

1. `pixi run containerize` runs `snakemake --containerize`. The command writes the `Dockerfile`.
2. The `Dockerfile` contains the label `io.github.snakemake.conda_env_hash`. Its value is the env hash.
3. `container_pin.py update` writes the tag `env-<hash>` into the `containerized:` line of `workflow/Snakefile`.
4. `build-3t-seq.yml` builds the image from the `Dockerfile` and tags it `env-<hash>`.

### Why the pin must match

Snakemake finds each conda env inside the image by its content hash.
When you change a conda env, its hash changes.
If the pinned image does not contain the new env, Snakemake stops with an error similar to this one:

```text
Unable to find environment in container image
```

The text of the error can be different in other Snakemake versions.

The `env-<hash>` tag changes when any conda env changes.
Thus, each commit pins an image that contains exactly its conda envs.

### Procedure: change a conda env

1. Edit the file in `workflow/env/`.
2. Run `pixi run containerize`.
3. Make sure that `git diff` shows changes in `Dockerfile` and `workflow/Snakefile`.
4. Commit the env file, `Dockerfile` and `workflow/Snakefile` in the same commit.
5. Push the branch and open a PR from a branch in this repository.
6. Wait for `build-3t-seq.yml` to build `env-<hash>`. The `check` job waits for this build.

Do not edit the `Dockerfile` or the `containerized:` line by hand.
The `dockerfile-sync` job regenerates both files and fails on a difference.

Make sure that conda can solve each env for `linux/amd64` and `linux/arm64`.
`build-3t-seq.yml` builds both platforms.

## SalmonTE image

The SalmonTE image has its own `Dockerfile` in `assets/docker/salmonTE/`.
It does not use the conda envs, so its tag is the release version.

- Release-please updates the version in `workflow/rules/salmonTE.smk`.
  The lines have the marker `# x-release-please-version`.
- `build-salmonte.yml` builds `:X.Y.Z` from the release PR, before the PR merges.
- After the git tag `vX.Y.Z` exists, the workflow does not overwrite `:X.Y.Z`.

A change to `assets/docker/salmonTE/Dockerfile` does not change the pin.
PRs and the `e2e` job use the old `:X.Y.Z` image until the next release PR.
To test a new SalmonTE image earlier, do these steps:

1. Run `gh workflow run build-salmonte.yml --ref <branch>`. Do not give a `version` input.
2. Find the `<branch>-<short-sha>` tag in the build log. The workflow replaces `/` in the branch name with `-`.
3. Pin that tag in `workflow/rules/salmonTE.smk` on your branch, and run the tests.
4. Put back the `X.Y.Z` pin before the PR merges.

## When CI builds images

| Event | 3t-seq image | SalmonTE image |
| :--- | :--- | :--- |
| PR to `dev` or `main` that changes `Dockerfile` or the workflow file (branch in this repo) | Builds `env-<hash>` if it does not exist. Otherwise, retags the existing image. | No |
| Release PR to `main` | Only if `Dockerfile` changed | Builds `:X.Y.Z` |
| Push to `dev` or `main` | Builds or retags, only when `Dockerfile` or the workflow file changes | Builds, only when `assets/docker/salmonTE/Dockerfile` or the workflow file changes |
| Push of git tag `vX.Y.Z` | Adds tag `:X.Y.Z` to `env-<hash>` | Builds `:X.Y.Z` only if it is missing |
| Manual run (`workflow_dispatch`) | Builds or retags for the selected ref | Builds the selected ref, or `:<version>` |
| Fork PR | No | No |

When `env-<hash>` already exists, `build-3t-seq.yml` does not build again.
It points the other tags at the existing image.

## Tags

| Tag | Image | Set by | Can it change? |
| :--- | :--- | :--- | :--- |
| `env-<hash>` | 3t-seq | The first build of that env hash | No. CI never overwrites it. |
| `X.Y.Z` | 3t-seq | Push of git tag `vX.Y.Z` | No. CI skips it when it exists. |
| `X.Y.Z` | SalmonTE | Release PR or manual run | Yes, until git tag `vX.Y.Z` exists. Then no. |
| `<branch>`, `<branch>-<short-sha>` | Both | Push or manual run on a branch. PR builds do not set these tags. | Yes. It follows the branch. |
| `sha-<sha>` | Both | Each build or retag | No in practice. One commit gives one tag. |
| `latest` | Both | Push or manual run on `main` | Yes. |

Pipeline code must pin `env-<hash>` or `X.Y.Z`. The other tags are for people, not for pins.

### `latest`

- Only a workflow run on `main` moves `latest`.
- `build-3t-seq.yml` runs on a push to `main` only when the `Dockerfile` changes.
  Thus, `latest` can stay on an older image after a release.
- To move `latest` to the current `main`, run `gh workflow run build-3t-seq.yml --ref main`.
