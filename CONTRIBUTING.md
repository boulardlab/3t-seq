# Contributing to 3t-seq

Thank you for your help. This file gives the short version.
The full rules are in the [Developer Guide](docs/developer/index.md).

1. Open an issue that describes the problem or the feature.
2. Create a branch from `dev`, for example `fix/42-short-name`.
3. Write commit messages as [Conventional Commits](https://www.conventionalcommits.org/), for example `fix: ...` or `feat: ...`.
4. Run the local checks. See [Local checks](docs/developer/local-checks.md).
5. If you changed a conda env, run `pixi run containerize`. Commit `Dockerfile` and `workflow/Snakefile` together.
6. Update the documentation in `docs/` if your change affects users.
7. Open a PR to `dev`, not to `main`.

PRs from forks have limits. CI does not build images for forks.
A maintainer must approve the CI runs of first-time contributors.
Read [Fork PRs](docs/developer/branches.md#fork-prs).

Pages in the Developer Guide:

- [Branches and PRs](docs/developer/branches.md)
- [CI reference](docs/developer/ci.md)
- [Containers](docs/developer/containers.md)
- [Release procedure](docs/developer/release.md)
- [Troubleshooting](docs/developer/troubleshooting.md)
