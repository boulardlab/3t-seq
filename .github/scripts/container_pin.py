#!/usr/bin/env python3
"""Keep the Snakefile `containerized:` pin in sync with the generated Dockerfile.

The 3t-seq image is tagged ``env-<hash>``, where ``<hash>`` is the
``io.github.snakemake.conda_env_hash`` label that ``snakemake --containerize``
writes into the Dockerfile. Snakemake looks environments up inside the image by
the same content hash, so a commit that changes a conda env needs a new image,
and this tag makes the pin change with it.

Usage:
    container_pin.py tag      print the image tag for the current Dockerfile
    container_pin.py update   rewrite the Snakefile pin to that tag
    container_pin.py check    exit 1 if the Snakefile pin does not match
"""

import re
import sys
from pathlib import Path

DOCKERFILE = Path("Dockerfile")
SNAKEFILE = Path("workflow/Snakefile")
IMAGE = "ghcr.io/boulardlab/3t-seq"

HASH_LABEL_RE = re.compile(r'^LABEL io\.github\.snakemake\.conda_env_hash="([0-9a-f]+)"$', re.M)
PIN_RE = re.compile(r'^(containerized: "docker://' + re.escape(IMAGE) + r':)([^"]+)(")', re.M)


def image_tag(dockerfile: Path = DOCKERFILE) -> str:
    """Return the ``env-<hash>`` tag derived from the Dockerfile label.

    Raises:
        ValueError: if the Dockerfile has no conda_env_hash label.
    """
    match = HASH_LABEL_RE.search(dockerfile.read_text())
    if match is None:
        raise ValueError(f"no conda_env_hash label in {dockerfile}; run `pixi run containerize`")
    return f"env-{match.group(1)}"


def current_pin(snakefile: Path = SNAKEFILE) -> str:
    """Return the tag currently pinned by the Snakefile ``containerized:`` line.

    Raises:
        ValueError: if no ``containerized:`` line for the 3t-seq image is found.
    """
    match = PIN_RE.search(snakefile.read_text())
    if match is None:
        raise ValueError(f'no `containerized: "docker://{IMAGE}:<tag>"` line in {snakefile}')
    return match.group(2)


def update_pin(tag: str, snakefile: Path = SNAKEFILE) -> bool:
    """Rewrite the Snakefile pin to ``tag``; return True if the file changed."""
    current_pin(snakefile)  # fail loudly if the line is missing
    text = snakefile.read_text()
    new_text = PIN_RE.sub(lambda m: f"{m.group(1)}{tag}{m.group(3)}", text, count=1)
    if new_text == text:
        return False
    snakefile.write_text(new_text)
    return True


def main(argv: list[str]) -> int:
    """Dispatch the ``tag`` / ``update`` / ``check`` subcommands."""
    if len(argv) != 2 or argv[1] not in {"tag", "update", "check"}:
        print(__doc__, file=sys.stderr)
        return 2
    try:
        tag = image_tag()
        if argv[1] == "tag":
            print(tag)
        elif argv[1] == "update":
            changed = update_pin(tag)
            print(f"{SNAKEFILE}: {IMAGE}:{tag}" + ("" if changed else " (unchanged)"))
        else:
            pinned = current_pin()
            if pinned != tag:
                print(
                    f"::error::{SNAKEFILE} pins {IMAGE}:{pinned} but the Dockerfile needs {tag}. "
                    "Run `pixi run containerize` and commit both files.",
                    file=sys.stderr,
                )
                return 1
            print(f"{SNAKEFILE} pin matches Dockerfile: {tag}")
    except (OSError, ValueError) as err:
        print(f"::error::{err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
