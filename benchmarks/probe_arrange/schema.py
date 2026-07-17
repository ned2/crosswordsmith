#!/usr/bin/env python3
"""Validate arrange-probe JSONL and guard against mixed-rig aggregation."""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict

REQUIRED = {
    "rig", "limit_kind", "operation_id", "attempt_index", "fixture",
    "fixture_seed", "search_seed", "corner", "arm", "cutoff", "termination",
    "outcome",
    "success_inferences", "censored", "max_depth", "places", "unplaces",
    "wipeouts", "reward", "layout_signature", "swi_version", "commit",
}
ENUMS = {
    "rig": {"authority", "instrumented"},
    "limit_kind": {"inferences", "nodes", "decisions", "none"},
    "termination": {"ok", "budget", "exhausted", "interrupted"},
    "outcome": {"placed", "not_proven", "infeasible", "interrupted"},
}
COUNTERS = ("max_depth", "places", "unplaces", "wipeouts", "nodes", "decisions")


def validate(row: dict) -> None:
    missing = REQUIRED - row.keys()
    if missing:
        raise ValueError(f"missing fields: {sorted(missing)}")
    for key, allowed in ENUMS.items():
        if row[key] not in allowed:
            raise ValueError(f"{key}={row[key]!r}, expected one of {sorted(allowed)}")
    if not isinstance(row["operation_id"], str) or not row["operation_id"]:
        raise ValueError("operation_id must be a non-empty string")
    if not isinstance(row["attempt_index"], int) or row["attempt_index"] < 0:
        raise ValueError("attempt_index must be a non-negative integer")
    if not isinstance(row["censored"], bool):
        raise ValueError("censored must be boolean")
    cutoff = row["cutoff"]
    if row["limit_kind"] == "none":
        if cutoff is not None:
            raise ValueError("limit_kind=none requires cutoff=null")
    elif isinstance(cutoff, bool) or not isinstance(cutoff, int) or cutoff < 0:
        raise ValueError("limited rows require a non-negative integer cutoff")
    if row["outcome"] == "placed":
        if row["reward"] is None or row["layout_signature"] is None:
            raise ValueError("placed rows require reward and layout_signature")
    elif row["reward"] is not None or row["layout_signature"] is not None:
        raise ValueError("non-placed rows require null reward/layout_signature")
    if row["rig"] == "authority":
        if row["limit_kind"] != "inferences":
            raise ValueError("authority rows require limit_kind=inferences")
        if any(row.get(field) is not None for field in COUNTERS):
            raise ValueError("authority mechanism counters must be null")
        if row["outcome"] == "placed" and row["success_inferences"] is None:
            raise ValueError("placed authority rows require success_inferences")
    else:
        if row["success_inferences"] is not None:
            raise ValueError("instrumented success_inferences must be null")
        if row["limit_kind"] == "inferences":
            raise ValueError("instrumented rows cannot use inference limits")
    expected_termination = {
        "ok": ("placed", False),
        "budget": ("not_proven", True),
        "exhausted": ("infeasible", False),
        "interrupted": ("interrupted", True),
    }
    expected_outcome, expected_censored = expected_termination[row["termination"]]
    if (row["outcome"], row["censored"]) != (expected_outcome, expected_censored):
        raise ValueError(
            f"termination={row['termination']} requires outcome={expected_outcome} "
            f"and censored={expected_censored}"
        )
    for field in COUNTERS:
        value = row.get(field)
        if value is not None and (not isinstance(value, int) or value < 0):
            raise ValueError(f"{field} must be null or a non-negative integer")


def guard_groups(rows: list[dict], fields: list[str]) -> dict[tuple, list[dict]]:
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        key = tuple(row.get(field) for field in fields)
        groups[key].append(row)
    for key, group in groups.items():
        rigs = {row["rig"] for row in group}
        if len(rigs) != 1:
            raise ValueError(f"mixed-rig aggregation rejected for group {key}: {sorted(rigs)}")
    return groups


def load_jsonl(stream) -> list[dict]:
    rows = []
    for line_no, line in enumerate(stream, 1):
        if not line.strip():
            continue
        row = json.loads(line)
        try:
            validate(row)
        except ValueError as exc:
            raise ValueError(f"line {line_no}: {exc}") from exc
        rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs="-", default=["-"])
    parser.add_argument("--group-by", default="fixture,corner,arm")
    args = parser.parse_args()
    rows = []
    for name in args.file:
        if name == "-":
            rows.extend(load_jsonl(sys.stdin))
        else:
            with open(name, encoding="utf-8") as stream:
                rows.extend(load_jsonl(stream))
    guard_groups(rows, [part for part in args.group_by.split(",") if part])
    print(f"schema: OK ({len(rows)} rows; mixed-rig guard passed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
