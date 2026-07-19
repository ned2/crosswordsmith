"""Descriptive stats units and CLI contracts (spec §6.4, D11)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from xword import XwordError
from xword.board import Board, Cell, derive_words
from xword.formats import parse_board
from xword.stats import board_stats, format_stats_text, stats_text
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


def test_small_grid_exact_metrics():
    report = board_stats(board_from_mask(["AB#", "C#D", "EFG"]))
    assert report == {
        "blocks": 2,
        "cells": {
            "checked": 3,
            "filled": 7,
            "light": 7,
            "orphan": 0,
            "unchecked": 4,
            "unfilled": 0,
        },
        "connectivity": {"components": 1, "connected": True},
        "dimensions": {"height": 3, "width": 3},
        "fillState": "solved",
        "symmetry": {"deficit": 1, "rot180": False},
        "words": {
            "across": 2,
            "down": 2,
            "enumerations": [
                {"count": 2, "enumeration": "(2)"},
                {"count": 2, "enumeration": "(3)"},
            ],
            "lengths": [{"count": 2, "length": 2}, {"count": 2, "length": 3}],
            "total": 4,
        },
    }


@pytest.mark.parametrize(
    ("mask", "state"),
    [(["..", ".."], "blank"), (["A.", ".."], "partial"), (["AB", "CD"], "solved")],
)
def test_fill_states(mask, state):
    assert board_stats(board_from_mask(mask))["fillState"] == state


def test_orphans_are_separate_from_unchecked():
    report = board_stats(board_from_mask(["AB#", "###", "##C"]))
    assert report["cells"]["unchecked"] == 2
    assert report["cells"]["orphan"] == 1
    assert sum(report["cells"][key] for key in ("checked", "unchecked", "orphan")) == 3


def test_rebus_is_one_filled_cell_and_word_length_uses_cells():
    board = board_from_mask(["AB", "CD"])
    board.grid[0][0].letter = "AL"
    report = board_stats(board)
    assert report["cells"]["filled"] == 4
    assert report["words"]["lengths"] == [{"count": 4, "length": 2}]


def test_connectivity_and_symmetry_deficit():
    report = board_stats(board_from_mask(["A##", "###", "##B"]))
    assert report["connectivity"] == {"components": 2, "connected": False}
    assert report["symmetry"] == {"deficit": 0, "rot180": True}

    asymmetric = board_stats(board_from_mask(["A##", "###", "###"]))
    assert asymmetric["symmetry"] == {"deficit": 1, "rot180": False}


def test_all_blocks_are_blank_connected_and_symmetric():
    report = board_stats(board_from_mask(["##", "##"]))
    assert report["fillState"] == "blank"
    assert report["connectivity"] == {"components": 0, "connected": True}
    assert report["symmetry"] == {"deficit": 0, "rot180": True}


def test_bars_fail_instead_of_reporting_block_derived_entries():
    board = board_from_mask(["AB", "CD"])
    board.grid[0][0].bar_right = True
    with pytest.raises(XwordError, match="bar-aware"):
        board_stats(board)

    styled = board_from_mask(["AB", "CD"])
    styled.grid[0][1].style = {"barred": "L"}
    with pytest.raises(XwordError, match="bar-aware"):
        board_stats(styled)


def test_native_and_ipuz_have_identical_stats(fixtures: Path):
    native = parse_board((fixtures / "bundled_17.native.json").read_text())
    ipuz = parse_board((fixtures / "bundled_17.ipuz.json").read_text())
    exolve = parse_board((fixtures / "bundled_17.exolve").read_text())
    assert board_stats(native) == board_stats(ipuz) == board_stats(exolve)


def test_human_and_json_output_are_deterministic():
    board = board_from_mask(["AB", "CD"])
    report = board_stats(board)
    human = format_stats_text(report)
    assert human == format_stats_text(report)
    assert human.endswith("\n")
    encoded = dump_json(report)
    assert encoded == dump_json(json.loads(encoded))


def test_stats_cli_file_stdin_json_and_out(fixtures: Path, tmp_path: Path):
    native = fixtures / "bundled_17.native.json"
    positional = run("stats", str(native))
    piped = run("stats", "--json", stdin=native.read_text())
    out = tmp_path / "stats.json"
    written = run("stats", "--json", str(native), "--out", str(out))

    assert positional.returncode == 0
    assert positional.stdout.startswith("Dimensions: 17x17\n")
    assert piped.returncode == 0
    assert json.loads(piped.stdout)["words"]["total"] == 6
    assert written.returncode == 0
    assert written.stdout == ""
    assert out.read_text() == piped.stdout


def test_stats_cli_failure_leaves_out_absent(tmp_path: Path):
    out = tmp_path / "stats.json"
    result = run("stats", "--json", "--out", str(out), stdin="garbage")
    assert result.returncode == 1
    assert result.stdout == ""
    assert not out.exists()


def test_stats_text_format_override(fixtures: Path):
    text = (fixtures / "bundled_17.native.json").read_text()
    with pytest.raises(XwordError, match="ipuz"):
        stats_text(text, "ipuz")
