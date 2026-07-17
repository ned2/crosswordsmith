#!/usr/bin/env python3
import json
import pathlib
import subprocess
import unittest

from benchmarks.probe_arrange.schema import guard_groups, validate

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUN = ROOT / "benchmarks" / "probe_arrange" / "run.pl"


def row(rig="authority", limit_kind=None, cutoff=None):
    if limit_kind is None:
        limit_kind = "inferences" if rig == "authority" else "none"
    if cutoff is None and limit_kind != "none":
        cutoff = 500_000_000
    return {
        "rig": rig, "limit_kind": limit_kind,
        "operation_id": "op", "attempt_index": 0, "fixture": "f.pl",
        "fixture_seed": 11, "search_seed": None, "corner": None, "arm": "control",
        "cutoff": cutoff, "termination": "ok", "outcome": "placed",
        "success_inferences": 10 if rig == "authority" else None,
        "censored": False, "max_depth": None if rig == "authority" else 2,
        "places": None if rig == "authority" else 2,
        "unplaces": None if rig == "authority" else 0,
        "wipeouts": None if rig == "authority" else 0,
        "nodes": None if rig == "authority" else 2,
        "decisions": None if rig == "authority" else 2,
        "reward": 12, "layout_signature": "sig", "swi_version": "10.1.10",
        "commit": "a" * 40,
    }


def terminated(termination):
    value = row()
    outcome, censored = {
        "ok": ("placed", False),
        "budget": ("not_proven", True),
        "exhausted": ("infeasible", False),
        "interrupted": ("interrupted", True),
    }[termination]
    value.update(termination=termination, outcome=outcome, censored=censored)
    if outcome != "placed":
        value.update(reward=None, layout_signature=None, success_inferences=None)
    return value


def emitted(*args):
    result = subprocess.run(
        ["swipl", "-q", str(RUN), "--", *map(str, args)],
        cwd=ROOT, check=True, text=True, capture_output=True,
    )
    return json.loads(result.stdout)


class SchemaTests(unittest.TestCase):
    def test_all_limit_kinds(self):
        for value in (
            row("authority", "inferences", 500_000_000),
            row("instrumented", "nodes", 100),
            row("instrumented", "decisions", 200),
            row("instrumented", "none", None),
        ):
            validate(value)

    def test_limit_kind_cutoff_mismatches_rejected(self):
        cases = [
            (row("instrumented", "none", None), "none requires cutoff=null"),
            (row("instrumented", "nodes", 10), "non-negative numeric cutoff"),
            (row("instrumented", "decisions", 10), "non-negative numeric cutoff"),
            (row("instrumented", "nodes", -1), "non-negative numeric cutoff"),
            (row("authority"), "authority rows require limit_kind=inferences"),
            (row("instrumented", "inferences", 100), "cannot use inference limits"),
        ]
        cases[0][0]["cutoff"] = 0
        cases[1][0]["cutoff"] = None
        cases[2][0]["cutoff"] = "10"
        cases[4][0]["limit_kind"] = "nodes"
        for value, message in cases:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, message):
                validate(value)

    def test_placed_budget_exhausted_and_interrupted(self):
        for termination in ("ok", "budget", "exhausted", "interrupted"):
            validate(terminated(termination))

    def test_termination_outcome_mismatch_rejected(self):
        value = terminated("budget")
        value["outcome"] = "infeasible"
        with self.assertRaisesRegex(ValueError, "termination=budget requires"):
            validate(value)

    def test_authority_mechanism_fields_remain_null(self):
        value = row()
        value["places"] = 0
        with self.assertRaisesRegex(ValueError, "must be null"):
            validate(value)

    def test_instrumented_success_inferences_remain_null(self):
        value = row("instrumented")
        value["success_inferences"] = 1
        with self.assertRaisesRegex(ValueError, "must be null"):
            validate(value)

    def test_actual_emitted_authority_rows(self):
        placed = emitted("authority-operation", "fixtures/bundled_17_clues.pl",
                          17, 6, 0, "none", 500_000_000,
                          "schema-authority-placed", 0, "control", 120)
        budget = emitted("authority-operation", "fixtures/bundled_17_clues.pl",
                          17, 6, 0, "none", 1000,
                          "schema-authority-budget", 0, "control", 120)
        corner = emitted("authority-corner", "fixtures/bundled_17_clues.pl",
                         17, 6, 0, "none", 500_000_000, "topright",
                         "schema-authority-corner", 0, "control", 120)
        validate(placed)
        validate(budget)
        validate(corner)
        self.assertEqual((placed["cutoff"], placed["termination"], placed["outcome"]),
                         (500_000_000, "ok", "placed"))
        self.assertEqual((budget["cutoff"], budget["termination"], budget["outcome"]),
                         (1000, "budget", "not_proven"))
        self.assertEqual((corner["cutoff"], corner["termination"], corner["outcome"]),
                         (500_000_000, "ok", "placed"))

    def test_actual_emitted_instrumented_rows(self):
        unbounded = emitted("instrumented", "fixtures/bundled_17_clues.pl",
                            17, 6, 0, "none", "topright", "lean", "none", 0,
                            "schema-inst-none", 0, "control", 120)
        limited = emitted("instrumented", "fixtures/bundled_17_clues.pl",
                          17, 6, 0, "none", "topright", "lean", "decisions", 1,
                          "schema-inst-limit", 0, "control", 120)
        node_limited = emitted("instrumented", "fixtures/bundled_17_clues.pl",
                               17, 6, 0, "none", "topright", "lean", "nodes", 1,
                               "schema-inst-nodes", 0, "control", 120)
        validate(unbounded)
        validate(limited)
        validate(node_limited)
        self.assertEqual((unbounded["cutoff"], unbounded["termination"]),
                         (None, "ok"))
        self.assertEqual((limited["cutoff"], limited["termination"], limited["outcome"]),
                         (1, "budget", "not_proven"))
        self.assertEqual((node_limited["cutoff"], node_limited["termination"],
                          node_limited["outcome"]),
                         (1, "budget", "not_proven"))

    def test_actual_interruption_preserves_cutoff(self):
        value = emitted("authority-operation", "fixtures/ladder_09x09_16w.pl",
                        9, 16, 11, 7, 500_000_000,
                        "schema-interrupted", 0, "control", 0.000001)
        validate(value)
        self.assertEqual((value["cutoff"], value["termination"], value["outcome"],
                          value["censored"]),
                         (500_000_000, "interrupted", "interrupted", True))

    def test_mixed_rig_rejected(self):
        with self.assertRaisesRegex(ValueError, "mixed-rig"):
            guard_groups([row(), row("instrumented")], ["fixture", "corner", "arm"])


if __name__ == "__main__":
    unittest.main()
