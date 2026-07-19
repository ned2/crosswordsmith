#!/usr/bin/env python3

import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import check_ac_fill_12
import gen_grids
import score_fill


class ScoreFillTests(unittest.TestCase):
    def test_report_agreement_checks_all_sidecar_fields(self):
        stats = {"n": 2, "mean": 49.2, "min": 30, "below50": 1, "junk": 0}
        report = {
            "n": 2,
            "mean": 49.3,
            "min": 30,
            "belowThreshold": 1,
            "threshold": 50,
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.json"
            path.write_text(json.dumps(report))
            self.assertTrue(score_fill.report_agrees(path, stats))
            report["threshold"] = 30
            path.write_text(json.dumps(report))
            self.assertFalse(score_fill.report_agrees(path, stats))

    def test_quality_floor_rejects_regressions(self):
        case = check_ac_fill_12.CASES[0]
        passing = {
            "n": 54,
            "mean": 45.0,
            "min": 30,
            "below50": 22,
            "junk": 0,
        }
        self.assertEqual(check_ac_fill_12.quality_errors(case, passing), [])

        regressed = dict(passing, mean=44.9, min=29, junk=1)
        errors = check_ac_fill_12.quality_errors(case, regressed)
        self.assertEqual(len(errors), 3)

    def test_gate_uses_default_deterministic_path(self):
        case = check_ac_fill_12.CASES[0]
        command = check_ac_fill_12.fill_command(
            Path("crosswordsmith"),
            Path("stw.txt"),
            case,
            Path("fill.json"),
            Path("report.json"),
        )
        self.assertNotIn("--budget", command)
        self.assertNotIn("--seed", command)
        self.assertNotIn("--shuffle", command)

    def test_layout_validation_rejects_grid_word_disagreement(self):
        layout = json.loads(
            (check_ac_fill_12.ROOT / "tests/golden/fill_3.json").read_text()
        )
        mask = ["...", "...", "..."]
        check_ac_fill_12.validate_layout(layout, mask)

        layout["grid"][0][0]["letter"] = "X"
        with self.assertRaisesRegex(ValueError, "answer disagrees"):
            check_ac_fill_12.validate_layout(layout, mask)

    def test_layout_validation_enforces_canonical_numbers_and_fields(self):
        source = check_ac_fill_12.ROOT / "tests/golden/fill_3.json"
        mask = ["...", "...", "..."]

        layout = json.loads(source.read_text())
        layout["words"][1]["number"] = 99
        with self.assertRaisesRegex(ValueError, "row-major clue number"):
            check_ac_fill_12.validate_layout(layout, mask)

        layout = json.loads(source.read_text())
        del layout["grid"][1][1]["number"]
        with self.assertRaisesRegex(ValueError, "missing canonical fields"):
            check_ac_fill_12.validate_layout(layout, mask)

    def test_scorer_rejects_non_ascii_words(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "scores.txt"
            path.write_text("PLAIN\nSCORED;50\n")
            self.assertEqual(score_fill.load_scores(path)["PLAIN"], 1)

            path.write_text("CAFE;50\nCAFÉ;60\n")
            with self.assertRaisesRegex(ValueError, "supports ASCII words only"):
                score_fill.load_scores(path)


class GridGeneratorTests(unittest.TestCase):
    def test_amer11_comes_from_shipped_grid(self):
        shipped = json.loads((gen_grids.ROOT / "grids/amer11.json").read_text())
        self.assertEqual(gen_grids.grid_documents()["amer11"], shipped)


if __name__ == "__main__":
    unittest.main()
