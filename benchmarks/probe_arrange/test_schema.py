#!/usr/bin/env python3
import unittest

from benchmarks.probe_arrange.schema import guard_groups, validate


def row(rig="authority"):
    base = {
        "rig": rig, "limit_kind": "inferences" if rig == "authority" else "none",
        "operation_id": "op", "attempt_index": 0, "fixture": "f.pl",
        "fixture_seed": 11, "search_seed": None, "corner": None, "arm": "control",
        "cutoff": "ok", "outcome": "placed",
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
    return base


class SchemaTests(unittest.TestCase):
    def test_valid_authority_and_instrumented(self):
        validate(row())
        validate(row("instrumented"))

    def test_null_not_false_zero(self):
        bad = row()
        bad["places"] = 0
        with self.assertRaisesRegex(ValueError, "must be null"):
            validate(bad)

    def test_mixed_rig_rejected(self):
        with self.assertRaisesRegex(ValueError, "mixed-rig"):
            guard_groups([row(), row("instrumented")], ["fixture", "corner", "arm"])


if __name__ == "__main__":
    unittest.main()
