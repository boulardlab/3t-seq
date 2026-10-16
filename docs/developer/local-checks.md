# Local checks

Run the same checks as CI on your computer before you push.
All checks use [pixi](https://pixi.sh) tasks from `pixi.toml`.

## Environments

| Environment | Use it for |
| :--- | :--- |
| `default` | `containerize`, `build-conda-envs` |
| `dev` | `format`, `format-check`, `lint`, `dry-run`, `setup`, `test`, `test-report` |
| `docs` | `docs-build`, `docs-serve` |

Select the environment with `-e`, for example `pixi run -e dev lint`.

!!! warning
    `test` and `test-report` depend on the `setup` task. `setup` exists only in `dev` and `docs`.
    Thus, `pixi run test` without `-e dev` fails with "could not find the task 'setup'".

## Tasks

| Task | Command | CI equivalent |
| :--- | :--- | :--- |
| Format | `pixi run -e dev format` | None. This task changes files. |
| Format check | `pixi run -e dev format-check` | `lint.yml` |
| Lint | `pixi run -e dev lint <config> <profile>` | `lint.yml` |
| Dry run | `pixi run -e dev dry-run <config> <profile>` | None |
| Test data | `pixi run -e dev setup` | `test.yml` (checkout with LFS) |
| Integration test | `pixi run -e dev test <config> <profile>` | `test.yml`, job `run-integration-tests` |
| Containerize | `pixi run containerize` | `container-check.yml`, job `dockerfile-sync` |
| Pin check | `python3 .github/scripts/container_pin.py check` | `container-check.yml`, job `check` |
| Image list | `bash .github/scripts/check-containers.sh --list` | `container-check.yml`, job `check` |
| Docs | `pixi run -e docs docs-build` | `docs.yml` |
| Pre-commit hooks | `pixi run -e dev hooks` | None |

The tasks `lint`, `dry-run`, `test` and `test-report` accept two positional arguments:

| Argument | Values | Default |
| :--- | :--- | :--- |
| `config` | `remote-references`, `local-references`, `adaptive-trim` | `remote-references` (`lint`, `dry-run`), `local-references` (`test`, `test-report`) |
| `profile` | `laptop`, `hpc` | `laptop` |

The files are in `.tests/integration/configs/<config>.yaml` and `.tests/integration/profiles/<profile>/`.

### Extra Snakemake arguments

To give more arguments to Snakemake, write `--` after the task arguments:

```bash
pixi run -e dev test adaptive-trim laptop -- --cores 4
```

Without `--`, pixi stops with "received more arguments than expected".

## Recommended sequence

Before you push a workflow change, run these commands:

1. Run `pixi run -e dev format`.
2. Run `pixi run -e dev format-check`.
3. Run `pixi run -e dev lint`.
4. Run `pixi run -e dev dry-run` for each config that your change affects.
5. If you changed a conda env, run `pixi run containerize`.
6. Run `pixi run -e dev test adaptive-trim laptop`. CI runs this config.

The integration test downloads references and runs containers. It needs Apptainer or Singularity.

## Integration test on Linux

The `default` environment installs Apptainer from conda-forge on `linux-64`.

On Ubuntu 24.04 and later, AppArmor blocks unprivileged user namespaces.
Rootless Apptainer needs them. To allow them until the next reboot, run:

```bash
sudo sysctl -w kernel.apparmor_restrict_unprivileged_userns=0
```

CI does the same step in `test.yml`.

## Integration test on macOS

Apptainer does not run on macOS. Use a Linux virtual machine (VM) with [Lima](https://lima-vm.io/).
Snakemake calls `singularity` on the host, and a wrapper script sends the call to the VM.

### Set up the VM

1. Install Lima: `brew install lima`.
2. Write a Lima template for an Ubuntu VM. Set `arch: x86_64` to run `linux/amd64` images.
3. In the template, mount your home directory with `writable: true`.
   Snakemake and the containers must see the same paths on the host and in the VM.
4. In the template, add a `provision` script that installs Apptainer or SingularityCE.
5. Create and start the VM, for example `limactl start --name=singularity-amd64 <template>.yaml`.
6. Make sure that the project directory is inside the mounted directory.

!!! tip
    An `x86_64` VM on Apple Silicon uses emulation, so the integration test is slow.
    CI runs the test only on `linux/amd64`.

### Add the wrapper

1. Create a file named `singularity` in a directory on your `PATH`, for example `~/.local/bin`.
2. Put this text in the file. Change the instance name and the binary path if yours are different.

    ```bash
    #!/bin/bash
    exec limactl shell singularity-amd64 -- /usr/bin/apptainer "$@"
    ```

3. Make the file executable: `chmod +x ~/.local/bin/singularity`.
4. Run `singularity --version`. The command must print the version from the VM.

Then run the integration test as usual:

```bash
pixi run -e dev test adaptive-trim laptop
```

The `laptop` profile binds `$PIXI_PROJECT_ROOT` and `$HOME` into the containers.
For more macOS notes, read [Apple Silicon & macOS notes](../getting-started/index.md#apple-silicon-macos-notes).
