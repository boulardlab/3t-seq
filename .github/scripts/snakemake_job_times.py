#!/usr/bin/env python3
"""Summarise job wall times per rule from a Snakemake log.

Usage: snakemake_job_times.py <snakemake.log>

Prints a Markdown table (rule, jobs, total and max minutes), sorted by total
time, followed by the overall wall time. Jobs that never finished are listed
as running.
"""

import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

TIMESTAMP = re.compile(r"^\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]$")
RULE = re.compile(r"^(?:local)?(?:checkpoint|rule) (\S+):$")
JOBID = re.compile(r"^\s+jobid: (\d+)$")
FINISHED = re.compile(r"^Finished job (\d+)\.$")


def parse(log: Path) -> tuple[dict[str, list[float]], list[str], float]:
    """Return per-rule durations (min), unfinished rules and total wall time (min)."""
    now: datetime | None = None
    first: datetime | None = None
    rule: str | None = None
    started: dict[str, tuple[str, datetime]] = {}
    durations: dict[str, list[float]] = defaultdict(list)
    for line in log.read_text(errors="replace").splitlines():
        if m := TIMESTAMP.match(line):
            now = datetime.strptime(m.group(1), "%a %b %d %H:%M:%S %Y")
            first = first or now
        elif m := RULE.match(line):
            rule = m.group(1)
        elif (m := JOBID.match(line)) and rule and now:
            started[m.group(1)] = (rule, now)
            rule = None
        elif (m := FINISHED.match(line)) and now and m.group(1) in started:
            name, t0 = started.pop(m.group(1))
            durations[name].append((now - t0).total_seconds() / 60)
    total = (now - first).total_seconds() / 60 if now and first else 0.0
    return durations, [r for r, _ in started.values()], total


def main() -> None:
    durations, running, total = parse(Path(sys.argv[1]))
    print("| rule | jobs | total min | max min |")
    print("|---|---:|---:|---:|")
    for name, d in sorted(durations.items(), key=lambda kv: -sum(kv[1])):
        print(f"| {name} | {len(d)} | {sum(d):.1f} | {max(d):.1f} |")
    for name in running:
        print(f"| {name} (unfinished) | 1 | | |")
    print(f"\nSnakemake wall time: {total:.1f} min")


if __name__ == "__main__":
    main()
