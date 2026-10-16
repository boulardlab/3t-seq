# Release procedure

This page gives the steps to release version `X.Y.Z` from `dev`.
[Release-please](https://github.com/googleapis/release-please) makes the release PR, the git tag and the GitHub release.

## Overview

| Step | Who | Result |
| :--- | :--- | :--- |
| 1. Prepare `dev` | Maintainer | `dev` passes all checks. |
| 2. PR `dev` → `main` | Maintainer + second approver | One squash commit on `main` |
| 3. Release PR | release-please | PR with the new version and changelog |
| 4. Merge the release PR | Maintainer + second approver | Git tag `vX.Y.Z`, GitHub release, image tags `:X.Y.Z` |
| 5. After the release | Maintainer | `dev` contains `main` again |

## Secret `RELEASE_PLEASE_TOKEN`

Release-please needs a PAT. With `GITHUB_TOKEN`, GitHub does not start other workflows for the release PR and the git tag.
Thus, the images and the required checks do not run.

The token must have these properties:

| Property | Value |
| :--- | :--- |
| Type | Fine-grained PAT |
| Resource owner | `boulardlab` |
| Repository access | Only `boulardlab/3t-seq` |
| Permission "Contents" | Read and write |
| Permission "Pull requests" | Read and write |
| Permission "Issues" | Read and write (release-please adds labels) |
| Expiration | A date. The token stops working on that date. |

The release PR has the token owner as author. Thus, the token owner cannot approve the release PR.

### Replace the token

Do these steps before the token expires, or when release-please reports `Bad credentials`:

1. Create a new fine-grained PAT with the properties in the table.
2. If the organization requires it, ask an organization owner to approve the token.
3. Run `gh secret set RELEASE_PLEASE_TOKEN --repo boulardlab/3t-seq`. Paste the token.
4. Run `gh run rerun <run-id> --failed` for the failed `release-please` run.
5. Write the new expiration date in your calendar.

## Step 1: Prepare `dev`

1. Merge all PRs for the release into `dev`.
2. Make sure that the last CI runs on `dev` passed.
3. If `main` has commits that `dev` does not have, merge `main` into `dev`.
4. Make sure that the `RELEASE_PLEASE_TOKEN` is still valid.
5. Run the dry runs for all configs:

    ```bash
    pixi run -e dev dry-run remote-references laptop
    pixi run -e dev dry-run local-references laptop
    pixi run -e dev dry-run adaptive-trim laptop
    ```

6. Run `pixi run -e dev test adaptive-trim laptop`.

## Step 2: PR from `dev` to `main`

1. Open the PR: `gh pr create --base main --head dev`.
2. Give the PR a Conventional Commits title, for example `feat: 3t-seq 2.0`.
   The title becomes the subject of the squash commit.
3. Wait for the checks `check`, `dockerfile-sync` and `run-integration-tests` to pass.
4. If GitHub says that the branch is out of date, click **Update branch**. Wait for the checks again.
5. Ask a second maintainer to approve the PR.
6. Select **Squash and merge**.
7. To force a version, end the commit message with this line:

    ```text
    Release-As: X.Y.Z
    ```

8. Confirm the merge.

Release-please selects the version from the commit types.
`feat` gives a minor version, `fix` gives a patch version.
A `!` after the type, or a `BREAKING CHANGE:` footer, gives a major version.
Use `Release-As` when these rules do not give the version that you want.

## Step 3: Release PR

The push to `main` starts `release-please.yml`. Release-please opens a PR from a `release-please--*` branch.
The release PR changes these files:

- `CHANGELOG.md`: a new section for `X.Y.Z`.
- `.release-please-manifest.json`: the new version.
- `pixi.toml`: the `version` line.
- `workflow/rules/salmonTE.smk`: the SalmonTE pins, `salmonte:X.Y.Z`.

The release PR starts these workflows:

- `build-salmonte.yml` builds `ghcr.io/boulardlab/3t-seq/salmonte:X.Y.Z`.
- `container-check.yml` waits for that image.
- `test.yml` runs the integration test with the new SalmonTE pin.

Then do these steps:

1. Read the new `CHANGELOG.md` section in the release PR.
2. A squash merge gives a changelog with one entry. Edit the section to add the real release notes.
   Commit the edit to the release PR branch.
3. Wait for the three required checks to pass.

## Step 4: Merge the release PR

1. Ask a maintainer who is not the token owner to approve the release PR.
   Alternatively, an admin can bypass the rule.
2. Squash-merge the release PR.
3. Wait for `release-please.yml` on `main`. It creates the git tag `vX.Y.Z` and the GitHub release.
4. Copy the release notes into the GitHub release body if they are missing.

The git tag `vX.Y.Z` starts the image workflows:

- `build-3t-seq.yml` adds the tag `:X.Y.Z` to the existing `env-<hash>` image.
- `build-salmonte.yml` does nothing, because `salmonte:X.Y.Z` exists.

## Step 5: After the release

1. Make sure that the images exist:

    ```bash
    docker manifest inspect ghcr.io/boulardlab/3t-seq:X.Y.Z
    docker manifest inspect ghcr.io/boulardlab/3t-seq/salmonte:X.Y.Z
    ```

2. Merge `main` back into `dev`:

    ```bash
    git fetch origin
    git switch dev
    git merge --ff-only origin/dev
    git merge origin/main
    git push origin dev
    ```

3. If `latest` must show the release, run `gh workflow run build-3t-seq.yml --ref main`.
   See [`latest`](containers.md#latest).
4. Close the milestone and the release tracking issue.

!!! note "Release 2.0.0"
    Issue #24 tracks the extra steps for 2.0.0. Example: close the old release PR #11.
