#!/usr/bin/env python3
"""Permanent external-dictionary gate for design-spec AC-FILL-12."""

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

import score_fill


ROOT = Path(__file__).resolve().parents[2]
KNOWN_STW_SHA256 = (
    "9aa4563eb82f9af45f3a9c4fbc7fe3c8c3e27da6b76cc48f5e09c499bd17e4ee"
)
CASES = (
    {
        "label": "blocked_13a@30",
        "min_score": 30,
        "expected_n": 54,
        "mean_floor": 45.0,
        "min_floor": 30,
    },
    {
        "label": "blocked_13a@1",
        "min_score": 1,
        "expected_n": 54,
        "mean_floor": 38.7,
        "min_floor": 1,
    },
)


class GateFailure(Exception):
    pass


def fingerprint(path):
    digest = hashlib.sha256()
    lines = 0
    with path.open("rb") as stream:
        for line in stream:
            digest.update(line)
            lines += 1
    return digest.hexdigest(), lines


def fill_command(crosswordsmith, dictionary, case, layout, report):
    return [
        str(crosswordsmith),
        "fill",
        "--grid", str(ROOT / "grids/blocked_13a.json"),
        "--dict", str(dictionary),
        "--min-score", str(case["min_score"]),
        "--report-json", str(report),
        "--out", str(layout),
    ]


def quality_errors(case, observed):
    errors = []
    if observed["n"] != case["expected_n"]:
        errors.append(
            f"n={observed['n']} (expected {case['expected_n']})"
        )
    if observed["mean"] < case["mean_floor"]:
        errors.append(
            f"mean={observed['mean']} (floor {case['mean_floor']:.1f})"
        )
    if observed["min"] < case["min_floor"]:
        errors.append(
            f"min={observed['min']} (floor {case['min_floor']})"
        )
    if observed["junk"]:
        errors.append(f"junk={observed['junk']} (expected 0)")
    return errors


def validate_layout(doc, mask):
    """Validate a canonical complete fill against the reference mask."""
    size = len(mask)
    if not isinstance(doc, dict) or doc.get("gridLength") != size:
        raise ValueError(f"gridLength must be {size}")

    grid = doc.get("grid")
    if (not isinstance(grid, list) or len(grid) != size
            or any(not isinstance(row, list) or len(row) != size for row in grid)):
        raise ValueError(f"grid must be {size}x{size}")

    for row, mask_row in zip(grid, mask):
        for cell, mark in zip(row, mask_row):
            if mark == "#":
                if cell is not None:
                    raise ValueError("blocked mask cell is not null")
                continue
            if not isinstance(cell, dict):
                raise ValueError("open mask cell is not a canonical cell object")
            required = {"letter", "number", "across", "down"}
            if not required.issubset(cell):
                raise ValueError("open cell is missing canonical fields")
            letter = cell.get("letter")
            if (not isinstance(letter, str) or len(letter) != 1
                    or not ("A" <= letter <= "Z")):
                raise ValueError("open cell does not contain one A-Z letter")
            for field in ("number", "across", "down"):
                value = cell.get(field)
                if value is not None and type(value) is not int:
                    raise ValueError(f"cell {field} must be an integer or null")

    expected = set()
    for cells in score_fill.slots_from_mask(mask):
        direction = "across" if cells[0][0] == cells[1][0] else "down"
        expected.add((direction, tuple(cells)))
    start_numbers = {
        cell: number
        for number, cell in enumerate(
            sorted({cells[0] for _, cells in expected}), start=1
        )
    }
    expected_numbers = {
        slot: start_numbers[slot[1][0]]
        for slot in expected
    }

    words = doc.get("words")
    if not isinstance(words, list):
        raise ValueError("words must be an array")

    seen = set()
    answers = set()
    references = {}
    for word in words:
        if not isinstance(word, dict):
            raise ValueError("word entry must be an object")
        answer = word.get("answer")
        direction = word.get("direction")
        cells = word.get("cells")
        number = word.get("number")
        if not isinstance(answer, str) or direction not in ("across", "down"):
            raise ValueError("word needs a string answer and valid direction")
        if type(number) is not int or not isinstance(word.get("meta"), dict):
            raise ValueError("word needs an integer number and object meta")
        if not isinstance(cells, list) or not cells:
            raise ValueError("word cells must be a non-empty array")

        parsed_cells = []
        for pair in cells:
            if (not isinstance(pair, list) or len(pair) != 2
                    or any(type(value) is not int for value in pair)):
                raise ValueError("word cell must be an integer [row, col] pair")
            row, col = pair
            if not (0 <= row < size and 0 <= col < size):
                raise ValueError("word cell lies outside the grid")
            parsed_cells.append((row, col))

        slot = (direction, tuple(parsed_cells))
        if slot not in expected:
            raise ValueError("word cells are not a maximal slot in the mask")
        if slot in seen:
            raise ValueError("slot appears more than once")
        if number != expected_numbers[slot]:
            raise ValueError("word has the wrong row-major clue number")
        seen.add(slot)

        folded = score_fill.fold(answer)
        letters = "".join(grid[row][col]["letter"] for row, col in parsed_cells)
        if folded != letters:
            raise ValueError("word answer disagrees with grid letters")
        if folded in answers:
            raise ValueError("answer appears more than once")
        answers.add(folded)

        for cell in parsed_cells:
            refs = references.setdefault(cell, {})
            if direction in refs:
                raise ValueError("cell has duplicate direction references")
            refs[direction] = number

    if seen != expected:
        raise ValueError(f"layout covers {len(seen)} of {len(expected)} slots")

    for row in range(size):
        for col in range(size):
            if mask[row][col] == "#":
                continue
            cell = grid[row][col]
            refs = references.get((row, col), {})
            if cell.get("across") != refs.get("across"):
                raise ValueError("cell across reference disagrees with words")
            if cell.get("down") != refs.get("down"):
                raise ValueError("cell down reference disagrees with words")
            expected_number = start_numbers.get((row, col))
            if cell.get("number") != expected_number:
                raise ValueError("cell number disagrees with starting words")


def run_case(crosswordsmith, dictionary, scores, case, work):
    layout = work / f"{case['label']}.json"
    report = work / f"{case['label']}.report.json"
    command = fill_command(crosswordsmith, dictionary, case, layout, report)

    print(f"RUN  {case['label']} (default budget, no user seed)", flush=True)
    completed = subprocess.run(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise GateFailure(
            f"{case['label']}: fill exited {completed.returncode}"
            + (f": {detail}" if detail else "")
        )
    if not layout.is_file() or not report.is_file():
        raise GateFailure(f"{case['label']}: fill did not write layout and report")

    try:
        layout_doc = json.loads(layout.read_text())
        mask_doc = json.loads((ROOT / "grids/blocked_13a.json").read_text())
        validate_layout(layout_doc, mask_doc["mask"])
        words = score_fill.cs_words(layout)
        observed = score_fill.stats(words, scores)
        report_ok = score_fill.report_agrees(report, observed)
        json.loads(report.read_text())
    except (
        OSError, ValueError, KeyError, TypeError, AttributeError, IndexError
    ) as exc:
        raise GateFailure(f"{case['label']}: malformed fill/report: {exc}") from exc

    errors = quality_errors(case, observed)
    if not report_ok:
        errors.append("engine sidecar disagrees with independent post-hoc scoring")
    if errors:
        raise GateFailure(f"{case['label']}: " + "; ".join(errors))

    print(
        f"PASS {case['label']}: n={observed['n']} mean={observed['mean']:.1f} "
        f"min={observed['min']} below50={observed['below50']}",
        flush=True,
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dict",
        dest="dictionary",
        default=os.environ.get("STW"),
        help="STW-class word;score file (default: $STW)",
    )
    parser.add_argument(
        "--crosswordsmith",
        type=Path,
        default=ROOT / "crosswordsmith",
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args(argv)
    if not args.dictionary:
        parser.error(
            "set STW or pass --dict with a Spread the Wordlist word;score file "
            "(CC BY-NC-SA; not bundled)"
        )
    args.dictionary = Path(args.dictionary).expanduser().resolve()
    return args


def main(argv=None):
    args = parse_args(argv)
    if not args.dictionary.is_file():
        raise SystemExit(f"dictionary is not a readable file: {args.dictionary}")
    if not args.crosswordsmith.is_file():
        raise SystemExit(f"crosswordsmith executable not found: {args.crosswordsmith}")

    digest, lines = fingerprint(args.dictionary)
    print(f"dictionary: {args.dictionary}")
    print(f"dictionary lines: {lines}")
    print(f"dictionary sha256: {digest}")
    if digest != KNOWN_STW_SHA256:
        raise SystemExit(
            "dictionary does not match the recorded AC-FILL-12 reference "
            f"snapshot (expected sha256 {KNOWN_STW_SHA256})"
        )

    scores = score_fill.load_scores(args.dictionary)
    if not scores:
        raise SystemExit("dictionary contains no scored entries")

    try:
        with tempfile.TemporaryDirectory(prefix="crosswordsmith-ac-fill-12-") as tmp:
            work = Path(tmp)
            for case in CASES:
                run_case(
                    args.crosswordsmith.resolve(), args.dictionary, scores, case, work
                )
    except GateFailure as exc:
        print(f"FAIL AC-FILL-12: {exc}")
        return 1

    print("PASS AC-FILL-12")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
