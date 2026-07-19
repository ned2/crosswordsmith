"""Semantic diff units and CLI contracts (spec §6.5, D12)."""

import copy
import json
import subprocess
import sys
from pathlib import Path

import pytest

from xword.board import Board, Cell, derive_words
from xword.convert import convert_text
from xword.diff import diff_boards, diff_texts, format_diff_text
from xword.formats import parse_board
from xword.util import dump_json


def board_from_mask(mask: list[str]) -> Board:
    grid = [
        [None if ch == "#" else Cell(letter=ch if ch != "." else None) for ch in row]
        for row in mask
    ]
    return Board(height=len(grid), width=len(grid[0]), grid=grid, words=derive_words(grid))


def run(*args: str, stdin: str | None = None):
    return subprocess.run(
        [sys.executable, "-m", "xword.cli", *args],
        input=stdin,
        capture_output=True,
        text=True,
    )


def run_bytes(*args: str, stdin: bytes):
    return subprocess.run(
        [sys.executable, "-m", "xword.cli", *args],
        input=stdin,
        capture_output=True,
    )


def paths(report):
    return [difference["path"] for difference in report["differences"]]


def test_native_and_ipuz_are_semantically_equal(fixtures: Path):
    native = (fixtures / "bundled_17.native.json").read_text()
    ipuz = (fixtures / "bundled_17.ipuz.json").read_text()
    output, different = diff_texts(native, ipuz)
    assert not different
    assert output == "equal: no semantic differences\n"


def test_payload_metadata_answers_and_empty_clues_are_ignored(fixtures: Path):
    left = parse_board((fixtures / "bundled_17.native.json").read_text())
    right = copy.deepcopy(left)
    right.diagnostics = {"producer": {"note": "changed"}}
    for word in right.words:
        word.answer = None
        word.meta = {"private": "changed"}
    left.words[0].clue = None
    right.words[0].clue = ""
    assert diff_boards(left, right) == {"differences": [], "equal": True}


@pytest.mark.parametrize(
    ("mutation", "expected_path"),
    [
        (lambda b: setattr(b.grid[0][0], "letter", "Z"), "grid[0][0].letter"),
        (lambda b: setattr(b.grid[0][0], "number", 99), "grid[0][0].number"),
        (lambda b: setattr(b.grid[0][0], "across_num", 99), "grid[0][0].across"),
        (lambda b: setattr(b.grid[0][0], "down_num", 99), "grid[0][0].down"),
        (lambda b: setattr(b.grid[0][0], "circle", True), "grid[0][0].circle"),
        (lambda b: setattr(b.grid[0][0], "bar_right", True), "grid[0][0].barRight"),
        (lambda b: setattr(b.grid[0][0], "bar_below", True), "grid[0][0].barBelow"),
        (lambda b: setattr(b.grid[0][0], "prefilled", True), "grid[0][0].prefilled"),
        (lambda b: setattr(b.grid[0][0], "style", {"color": "red"}), "grid[0][0].style"),
        (lambda b: setattr(b.words[0], "clue", "Changed"), "words[across:0,0-0,1].clue"),
        (lambda b: setattr(b.words[0], "enumeration", "(1,1)"), "words[across:0,0-0,1].enumeration"),
        (lambda b: setattr(b.words[0], "number", 99), "words[across:0,0-0,1].number"),
        (lambda b: b.meta.__setitem__("title", "Changed"), "meta.title"),
        (lambda b: b.meta.__setitem__("author", "Changed"), "meta.author"),
    ],
)
def test_structural_fields_are_compared(mutation, expected_path):
    left = board_from_mask(["AB", "CD"])
    right = copy.deepcopy(left)
    mutation(right)
    assert expected_path in paths(diff_boards(left, right))


def test_block_topology_and_word_geometry_are_localised():
    left = board_from_mask(["AB", "CD"])
    right = board_from_mask(["A#", "CD"])
    report = diff_boards(left, right)
    assert "grid[0][1]" in paths(report)
    assert any(path.startswith("words[") for path in paths(report))


def test_missing_enumeration_defaults_to_cell_count():
    left = board_from_mask(["AB", "CD"])
    right = copy.deepcopy(left)
    right.words[0].enumeration = "(2)"
    assert diff_boards(left, right)["equal"]


def test_duplicate_ipuz_circle_and_bar_style_is_normalised():
    left = board_from_mask(["AB", "CD"])
    right = copy.deepcopy(left)
    left.grid[0][0].circle = right.grid[0][0].circle = True
    left.grid[0][0].bar_right = right.grid[0][0].bar_right = True
    right.grid[0][0].style = {"shapebg": "circle", "barred": "R"}
    assert diff_boards(left, right)["equal"]


def test_ipuz_left_and_right_bar_spellings_are_semantically_equal():
    base = {
        "version": "http://ipuz.org/v2",
        "kind": ["http://ipuz.org/crossword#1"],
        "dimensions": {"width": 2, "height": 2},
        "puzzle": [[0, 0], [0, 0]],
    }
    left_bar = copy.deepcopy(base)
    left_bar["puzzle"][0][1] = {"cell": 0, "style": {"barred": "L"}}
    right_bar = copy.deepcopy(base)
    right_bar["puzzle"][0][0] = {"cell": 0, "style": {"barred": "R"}}
    left = parse_board(json.dumps(left_bar), "ipuz")
    right = parse_board(json.dumps(right_bar), "ipuz")
    assert left.grid[0][0].bar_right
    assert diff_boards(left, right)["equal"]


def test_rebus_token_difference_is_reported_as_one_cell():
    left = board_from_mask(["AB", "CD"])
    right = copy.deepcopy(left)
    right.grid[0][0].letter = "AL"
    report = diff_boards(left, right)
    assert report["differences"] == [
        {"left": "A", "path": "grid[0][0].letter", "right": "AL"}
    ]


def test_exolve_shade_parse_loss_is_reported(fixtures: Path):
    source = (fixtures / "sample_shaded.ipuz.json").read_text()
    converted, _warnings = convert_text(source, "exolve", "ipuz")
    output, different = diff_texts(source, converted, "ipuz", "exolve")
    assert different
    assert ".style:" in output


def test_difference_order_and_json_are_stable():
    left = board_from_mask(["AB", "CD"])
    right = copy.deepcopy(left)
    right.height = 3
    right.grid[0][1].letter = "Z"
    right.words[0].clue = "Changed"
    right.meta["title"] = "Title"
    report = diff_boards(left, right)
    assert paths(report) == [
        "dimensions.height",
        "grid[0][1].letter",
        "words[across:0,0-0,1].clue",
        "meta.title",
    ]
    assert dump_json(report) == dump_json(json.loads(dump_json(report)))
    assert format_diff_text(report).startswith("different: 4 semantic differences\n")


def test_cli_exit_zero_for_equal_and_one_after_report(fixtures: Path, tmp_path: Path):
    native = fixtures / "bundled_17.native.json"
    ipuz = fixtures / "bundled_17.ipuz.json"
    equal = run("diff", str(native), str(ipuz))
    assert equal.returncode == 0
    assert equal.stdout == "equal: no semantic differences\n"

    changed = json.loads(ipuz.read_text())
    changed["title"] = "Changed"
    changed_path = tmp_path / "changed.ipuz"
    changed_path.write_text(json.dumps(changed))
    out = tmp_path / "diff.json"
    different = run("diff", "--json", str(native), str(changed_path), "--out", str(out))
    assert different.returncode == 1
    assert different.stdout == ""
    assert json.loads(out.read_text())["differences"][-1]["path"] == "meta.title"


def test_cli_one_stdin_operand_and_independent_format_overrides(fixtures: Path):
    native = fixtures / "bundled_17.native.json"
    ipuz = fixtures / "bundled_17.ipuz.json"
    result = run(
        "diff",
        "--from-left",
        "native",
        "--from-right",
        "ipuz",
        "-",
        str(ipuz),
        stdin=native.read_text(),
    )
    assert result.returncode == 0


def test_cli_rejects_two_stdin_operands_with_exit_two():
    result = run("diff", "-", "-", stdin="anything")
    assert result.returncode == 2
    assert result.stdout == ""
    assert "only one operand" in result.stderr


def test_cli_parse_error_exits_two_and_leaves_out_absent(tmp_path: Path):
    good = tmp_path / "good.json"
    bad = tmp_path / "bad.json"
    out = tmp_path / "diff.json"
    good.write_text('{"grid": []}')
    bad.write_text("garbage")
    result = run("diff", str(good), str(bad), "--out", str(out))
    assert result.returncode == 2
    assert result.stdout == ""
    assert not out.exists()


def test_cli_malformed_native_schema_and_invalid_utf8_exit_two(tmp_path: Path, fixtures: Path):
    malformed = tmp_path / "malformed.json"
    malformed.write_text('{"gridLength": 1, "grid": [[null]], "words": 7}')
    native = fixtures / "bundled_17.native.json"
    schema_result = run("diff", "--from-left", "native", str(malformed), str(native))
    assert schema_result.returncode == 2
    assert "words must be an array" in schema_result.stderr

    invalid_utf8 = tmp_path / "invalid.json"
    invalid_utf8.write_bytes(b"\xff")
    encoding_result = run("diff", str(invalid_utf8), str(native))
    assert encoding_result.returncode == 2
    assert "cannot read" in encoding_result.stderr

    stdin_result = run_bytes("diff", "-", str(native), stdin=b"\xff")
    assert stdin_result.returncode == 2
    assert b"cannot read stdin" in stdin_result.stderr


def test_cli_output_error_exits_two_not_difference(tmp_path: Path, fixtures: Path):
    native = fixtures / "bundled_17.native.json"
    result = run(
        "diff",
        str(native),
        str(native),
        "--out",
        str(tmp_path / "missing" / "diff.txt"),
    )
    assert result.returncode == 2
    assert result.stdout == ""
    assert "cannot write" in result.stderr


def test_cli_missing_file_exits_two(tmp_path: Path, fixtures: Path):
    result = run(
        "diff",
        str(tmp_path / "missing.json"),
        str(fixtures / "bundled_17.native.json"),
    )
    assert result.returncode == 2
    assert "cannot read" in result.stderr
