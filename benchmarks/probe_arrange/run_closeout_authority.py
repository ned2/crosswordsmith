#!/usr/bin/env python3
"""Run the final unseeded 500M operation authority battery serially."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "benchmarks/probe_arrange/run.pl"
BUDGET = 500_000_000

CASES = (
    ("guard-09x09-17w", "fixtures/ladder_09x09_17w.pl", 9, 17, 11, "guard", True),
    ("guard-15x15-40w", "fixtures/ladder_15x15_40w.pl", 15, 40, 11, "guard", True),
    ("guard-21x21-82w", "fixtures/ladder_21x21_82w.pl", 21, 82, 11, "guard", True),
    (
        "cliff-09x09-18w-seed12",
        "benchmarks/probe_arrange/fixtures/cliff_09x09_18w_seed12.pl",
        9,
        18,
        12,
        "cliff",
        False,
    ),
    (
        "cliff-15x15-44w-seed11",
        "benchmarks/probe_arrange/fixtures/cliff_15x15_44w_seed11.pl",
        15,
        44,
        11,
        "cliff",
        False,
    ),
    (
        "cliff-21x21-88w-seed11",
        "benchmarks/probe_arrange/fixtures/cliff_21x21_88w_seed11.pl",
        21,
        88,
        11,
        "cliff",
        False,
    ),
)


def main() -> int:
    for operation, fixture, grid, count, fixture_seed, role, must_place in CASES:
        print(f"closeout-authority heartbeat operation={operation} phase=start", file=sys.stderr, flush=True)
        command = [
            "swipl",
            "-q",
            str(RUNNER),
            "--",
            "authority-operation",
            fixture,
            str(grid),
            str(count),
            str(fixture_seed),
            "none",
            str(BUDGET),
            operation,
            "0",
            f"closeout_{role}",
            "120",
        ]
        completed = subprocess.run(
            command,
            cwd=ROOT,
            text=True,
            capture_output=True,
            timeout=150,
            check=True,
        )
        sys.stderr.write(completed.stderr)
        row = json.loads(completed.stdout)
        if must_place and row["outcome"] != "placed":
            raise RuntimeError(f"robust guard did not place: {operation}: {row['outcome']}")
        print(json.dumps(row, sort_keys=True), flush=True)
        print(f"closeout-authority heartbeat operation={operation} phase=done", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
