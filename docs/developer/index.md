# Developer Guide

This guide is for developers and maintainers of 3t-seq.
Use it when you change the workflow, a conda env, a container, the CI or a release.

To run the pipeline on your own data, read the user documentation:
[Getting Started](../getting-started/index.md), [Usage](../usage/index.md) and the
[Tutorial](../tutorial/quickstart.md). To move a 1.x project to 2.0, read
[Migrating from 1.x](../migration/index.md).

## Pages in this guide

| Page | Contents |
| :--- | :--- |
| [Branches and PRs](branches.md) | Branch model, PR flow, fork PRs and the ruleset for `main` |
| [Local checks](local-checks.md) | The pixi tasks to run before you push, and the macOS setup |
| [CI reference](ci.md) | Each workflow in `.github/workflows`: trigger, jobs, failures and fixes |
| [Containers](containers.md) | How the pipeline pins its images, and when CI builds them |
| [Release procedure](release.md) | The steps from `dev` to a `vX.Y.Z` tag |
| [Troubleshooting](troubleshooting.md) | Known errors: symptom, cause and action |

## Change checklist

Use this table to find the checks for your change.

| You change | Do this before you push |
| :--- | :--- |
| A rule or a script in `workflow/` | Run `pixi run -e dev format`, then `pixi run -e dev format-check` and `pixi run -e dev lint`. |
| A conda env in `workflow/envs/` | Run `pixi run containerize`. Commit `Dockerfile` and `workflow/Snakefile` together. See [Containers](containers.md). |
| The SalmonTE image in `assets/docker/salmonTE/` | Read [SalmonTE image](containers.md#salmonte-image). The new image ships with the next release. |
| The config schema | Read [Config schema](#config-schema). |
| A file in `docs/` or `mkdocs.yml` | Run `pixi run -e docs docs-build`. |
| A workflow in `.github/workflows/` | Read the [CI reference](ci.md) for the workflow. Open the PR from a branch in this repository. |

Also update the user documentation in `docs/` when your change affects users.

## Config schema

The file `workflow/schemas/config.schema.yaml` is the single source of truth for config keys and default values.

- `validate_config()` in `workflow/rules/common.smk` checks the config against the schema.
  It reports all errors together, then applies the schema defaults.
- The schema sets `additionalProperties: false` on most objects.
  Thus, a key that is not in the schema causes an error.
- When you add a config key, add it to the schema with a `default`.
  Then document it in [Configuration](../configuration/index.md).
- When you remove or rename a key, add a row to the [migration guide](../migration/index.md).
- `AUTO_FETCH_GENOMES` in `common.smk` lists the genomes that the pipeline can download (`mm10`, `mm39`).

!!! note
    `validate_config()` comes from PR #35. Before #35 merges, the Snakefile calls
    `snakemake.utils.validate()` and `apply_schema_defaults()`.

## Glossary

This guide uses these terms with one meaning only.

| Term | Meaning |
| :--- | :--- |
| PR | A GitHub pull request. |
| Base branch | The branch that a PR merges into (`dev` or `main`). |
| Feature branch | A short-lived branch for one change, for example `fix/20-lint`. |
| Fork PR | A PR from a branch in a fork, not in `boulardlab/3t-seq`. |
| Maintainer | A person with write or admin access to `boulardlab/3t-seq`. |
| Workflow | A YAML file in `.github/workflows/`. |
| Job | One job in a workflow, for example `dockerfile-sync`. A job reports one check on a PR. |
| Gate | A job that must pass before a PR can merge into `main`. |
| Ruleset | The GitHub rules for a branch. The ruleset for `main` is `main: release gate`. |
| Conda env | A YAML file in `workflow/envs/` that a rule uses. |
| Env hash | The hash of all conda envs. Snakemake writes it into the `Dockerfile` as the label `io.github.snakemake.conda_env_hash`. |
| Image | A container image in `ghcr.io/boulardlab/`. |
| Tag | An image tag (`:env-<hash>`, `:2.0.0`) or a git tag (`v2.0.0`). The text always says which. |
| Pin | The image reference with a fixed tag in `workflow/`. |
| Release PR | The PR that release-please opens on `main` to make a release. |
| PAT | A GitHub personal access token. |
