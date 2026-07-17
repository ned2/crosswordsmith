#!/usr/bin/env python3
"""Validate and summarize the fixed P-R0 authority pilot."""
from __future__ import annotations

import argparse
import json
import pathlib
import statistics
import sys
from collections import defaultdict

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmarks.probe_arrange.schema import guard_groups, load_jsonl  # noqa: E402

BASE = "1bccf47917fea074a404bf467e51865ce988b8b8"
CORNERS = ("topleft_across", "topright")
CLIFFS = (
    "cliff_09x09_18w_seed12.pl",
    "cliff_15x15_44w_seed11.pl",
    "cliff_21x21_84w_seed12.pl",
    "cliff_21x21_88w_seed11.pl",
)
CONTROLS = (
    "ladder_09x09_08w.pl",
    "ladder_15x15_12w.pl",
    "ladder_21x21_25w.pl",
)
CURVE_CUTOFFS = (1_000_000, 10_000_000, 50_000_000, 100_000_000, 250_000_000, 500_000_000)


def check_matrix(rows: list[dict]) -> None:
    if len(rows) != 131:
        raise SystemExit(f"expected 131 rows, found {len(rows)}")
    operation_ids = [row["operation_id"] for row in rows]
    if len(set(operation_ids)) != len(operation_ids):
        raise SystemExit("duplicate operation_id")
    for row in rows:
        if row["rig"] != "authority" or row["limit_kind"] != "inferences":
            raise SystemExit(f"non-authority inference row: {row['operation_id']}")
        if row["commit"] != BASE or row["cutoff"] != 500_000_000:
            raise SystemExit(f"provenance/cutoff mismatch: {row['operation_id']}")
        if row["attempt_index"] != 0 or row["seed_index"] not in range(16):
            raise SystemExit(f"attempt/partition mismatch: {row['operation_id']}")
    cliff_rows = [row for row in rows if row["subject_kind"] == "cliff"]
    controls = [row for row in rows if row["subject_kind"] == "control"]
    expected = {(fixture, corner, index) for fixture in CLIFFS for corner in CORNERS for index in range(16)}
    actual = {(row["fixture"], row["corner"], row["seed_index"]) for row in cliff_rows}
    if actual != expected:
        raise SystemExit(f"cliff matrix mismatch: missing={expected-actual}, extra={actual-expected}")
    if len(controls) != 3 or {row["fixture"] for row in controls} != set(CONTROLS):
        raise SystemExit("control matrix mismatch")
    guard_groups(rows, ["fixture", "corner", "arm"])


def group_summary(rows: list[dict]) -> dict:
    success = [row for row in rows if row["outcome"] == "placed"]
    inferences = [row["success_inferences"] for row in success]
    return {
        "rows": len(rows),
        "placed": len(success),
        "censored": sum(row["censored"] for row in rows),
        "exhausted": sum(row["termination"] == "exhausted" for row in rows),
        "interrupted": sum(row["termination"] == "interrupted" for row in rows),
        "success_inferences_min": min(inferences) if inferences else None,
        "success_inferences_median": statistics.median(inferences) if inferences else None,
        "success_inferences_max": max(inferences) if inferences else None,
        "distinct_success_layouts": len({row["layout_signature"] for row in success}),
        "reward_min": min((row["reward"] for row in success), default=None),
        "reward_median": statistics.median([row["reward"] for row in success]) if success else None,
        "reward_max": max((row["reward"] for row in success), default=None),
    }


def success_curves(rows: list[dict]) -> dict:
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["fixture"], row["corner"])].append(row)
    curves = {}
    for (fixture, corner), group in sorted(grouped.items()):
        curves[f"{fixture}:{corner}"] = {
            str(cutoff): sum(
                row["outcome"] == "placed" and row["success_inferences"] <= cutoff
                for row in group
            )
            for cutoff in CURVE_CUTOFFS
        }
    return curves


def complementarity(rows: list[dict]) -> dict:
    indexed = {(row["fixture"], row["corner"], row["seed_index"]): row for row in rows}
    result = {}
    for fixture in CLIFFS:
        counts = {"both": 0, "topleft_only": 0, "topright_only": 0, "neither": 0}
        for seed_index in range(16):
            left = indexed[(fixture, "topleft_across", seed_index)]["outcome"] == "placed"
            right = indexed[(fixture, "topright", seed_index)]["outcome"] == "placed"
            key = "both" if left and right else "topleft_only" if left else "topright_only" if right else "neither"
            counts[key] += 1
        left_total = counts["both"] + counts["topleft_only"]
        right_total = counts["both"] + counts["topright_only"]
        union = 16 - counts["neither"]
        counts.update(
            topleft_total=left_total,
            topright_total=right_total,
            either_corner=union,
            union_gain_over_best=union - max(left_total, right_total),
        )
        result[fixture] = counts
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=pathlib.Path)
    args = parser.parse_args()
    with args.results.open(encoding="utf-8") as stream:
        rows = load_jsonl(stream)
    check_matrix(rows)
    cliff_rows = [row for row in rows if row["subject_kind"] == "cliff"]
    controls = [row for row in rows if row["subject_kind"] == "control"]
    groups = defaultdict(list)
    for row in cliff_rows:
        groups[(row["fixture"], row["corner"])].append(row)
    summary = {
        "rows": len(rows),
        "cliff_rows": len(cliff_rows),
        "control_rows": len(controls),
        "runner_child_cpu_seconds": sum(row["runner_child_cpu_seconds"] for row in rows),
        "authority_search_cpu_seconds": sum(row["cpu_seconds"] for row in rows),
        "early_stop": False,
        "controls": {row["fixture"]: group_summary([row]) for row in controls},
        "fixture_corners": {
            f"{fixture}:{corner}": group_summary(group)
            for (fixture, corner), group in sorted(groups.items())
        },
        "success_curves": success_curves(cliff_rows),
        "complementarity": complementarity(cliff_rows),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
