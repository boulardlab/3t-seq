#!/usr/bin/env bash
# Verify that every ghcr.io/boulardlab/* image pinned under workflow/ uses an
# immutable tag (not :latest, not untagged) and exists in the registry.
#
# Usage: check-containers.sh [--list] [--wait-minutes N] [--head-sha SHA]
#   --list          only print the extracted refs, no registry access
#   --wait-minutes  poll for up to N minutes while a build workflow run for
#                   --head-sha is still in progress (PRs, where the builds run
#                   concurrently with this check); 0 = check once
#   --head-sha      commit whose build-* workflow runs to follow (needs `gh`
#                   and GH_TOKEN); without it, waits until the deadline
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

# Status of the build-* runs for head_sha: "failed", "running" or "idle".
build_status() {
  local wf states
  states=$(for wf in "${build_workflows[@]}"; do
    gh run list --workflow "$wf" --commit "$head_sha" \
      --json status,conclusion --jq '.[] | "\(.status) \(.conclusion)"'
  done)
  if [[ -z "$states" ]]; then
    echo idle
  elif grep -q ' failure$' <<< "$states"; then
    echo failed
  elif grep -qv '^completed ' <<< "$states"; then
    echo running
  else
    echo idle
  fi
}

start=$(date +%s)
deadline=$(( start + wait_minutes * 60 ))
# Build runs triggered by the same event may not be registered yet.
grace=$(( start + ${GRACE_SECONDS:-180} ))
while true; do
  mapfile -t missing < <(missing_refs)
  if [[ ${#missing[@]} -eq 0 ]]; then
    echo "All pinned images exist."
    exit 0
  fi
  now=$(date +%s)
  give_up=0
  if [[ $now -ge $deadline ]]; then
    give_up=1
  elif [[ -n "$head_sha" && $now -ge $grace ]]; then
    case "$(build_status)" in
      failed) echo "::error::a build workflow failed for $head_sha"; give_up=1 ;;
      idle) echo "No build running for $head_sha."; give_up=1 ;;
    esac
  fi
  if [[ $give_up -eq 1 ]]; then
    echo "::error::missing images in ghcr.io:"
    printf '  %s\n' "${missing[@]}"
    echo "Build them with the build-3t-seq / build-salmonte workflows (see their workflow_dispatch)."
    exit 1
  fi
  echo "Waiting for: ${missing[*]}"
  sleep 60
done
