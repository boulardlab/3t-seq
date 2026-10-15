#!/usr/bin/env bash
# Verify that every ghcr.io/boulardlab/* image pinned under workflow/ uses an
# immutable tag (not :latest, not untagged) and exists in the registry.
#
# Usage: check-containers.sh [--list] [--wait-minutes N] [--head-sha SHA]
#   --list          only print the extracted refs, no registry access
#   --wait-minutes  poll until all refs exist (release-please PRs, where the
#                   build workflows run concurrently); 0 = check once
#   --head-sha      with --wait-minutes, fail fast if a build workflow run for
#                   this commit concluded with failure (needs `gh` + GH_TOKEN)
set -euo pipefail

workflow_dir="workflow"
build_workflows=("build-3t-seq.yml" "build-salmonte.yml")
list_only=0
wait_minutes=0
head_sha=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --list) list_only=1; shift ;;
    --wait-minutes) wait_minutes="$2"; shift 2 ;;
    --head-sha) head_sha="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

# Covers the Snakefile `containerized:` line and every `container:` directive
# (the URI sits on the line after the directive). Commented-out lines skipped.
mapfile -t refs < <(
  grep -rhE 'ghcr\.io/boulardlab/' "$workflow_dir" \
    | grep -vE '^[[:space:]]*#' \
    | grep -oE 'ghcr\.io/boulardlab/[A-Za-z0-9._/:@-]+' \
    | sort -u
)

if [[ ${#refs[@]} -eq 0 ]]; then
  echo "::error::no ghcr.io/boulardlab/* image refs found under $workflow_dir/"
  exit 1
fi

echo "Pinned images:"
printf '  %s\n' "${refs[@]}"
[[ $list_only -eq 1 ]] && exit 0

# Immutability: digest pins are fine; tags must exist and must not be latest.
mutable=()
for ref in "${refs[@]}"; do
  [[ "$ref" == *@sha256:* ]] && continue
  name="${ref##*/}"
  if [[ "$name" != *:* ]]; then
    mutable+=("$ref (no tag)")
  elif [[ "${name##*:}" == "latest" ]]; then
    mutable+=("$ref (:latest)")
  fi
done
if [[ ${#mutable[@]} -gt 0 ]]; then
  echo "::error::mutable image refs (pin a version tag or digest):"
  printf '  %s\n' "${mutable[@]}"
  exit 1
fi

missing_refs() {
  local ref
  for ref in "${refs[@]}"; do
    docker manifest inspect "$ref" > /dev/null 2>&1 || echo "$ref"
  done
}

build_failed() {
  [[ -z "$head_sha" ]] && return 1
  local wf
  for wf in "${build_workflows[@]}"; do
    if gh run list --workflow "$wf" --commit "$head_sha" \
      --json conclusion --jq '.[].conclusion' | grep -qx failure; then
      echo "::error::$wf failed for $head_sha"
      return 0
    fi
  done
  return 1
}

deadline=$(( $(date +%s) + wait_minutes * 60 ))
while true; do
  mapfile -t missing < <(missing_refs)
  if [[ ${#missing[@]} -eq 0 ]]; then
    echo "All pinned images exist."
    exit 0
  fi
  if [[ $(date +%s) -ge $deadline ]] || build_failed; then
    echo "::error::missing images in ghcr.io:"
    printf '  %s\n' "${missing[@]}"
    exit 1
  fi
  echo "Waiting for: ${missing[*]}"
  sleep 60
done
