"""Board model units: derived numbering + enumeration normalisation (spec §5.1)."""

import pytest

from xword.board import (
    Cell,
    derive_words,
    enumeration_cell_count,
    enumeration_from_answer,
    extract_clue_enum,
    split_clue_enum,
)


def grid_from_mask(mask: list[str]):
    """'#' = block, anything else = a light cell holding that letter."""
    return [
        [None if ch == "#" else Cell(letter=ch if ch != "." else None) for ch in row]
        for row in mask
    ]


class TestDerivedNumbering:
    def test_row_major_scan_across_before_down(self):
        grid = grid_from_mask(
            [
                "AB#",
                "C#D",
                "EFG",
            ]
        )
        words = derive_words(grid)
        got = [(w.number, w.direction, w.cells) for w in words]
        assert got == [
            (1, "across", [(0, 0), (0, 1)]),
            (1, "down", [(0, 0), (1, 0), (2, 0)]),
            (2, "down", [(1, 2), (2, 2)]),
            (3, "across", [(2, 0), (2, 1), (2, 2)]),
        ]

    def test_single_cell_runs_are_not_words(self):
        # middle row is isolated single lights: no across or down runs of >= 2
        grid = grid_from_mask(
            [
                "A#B",
                "###",
                "C#D",
            ]
        )
        assert derive_words(grid) == []
        assert all(cell is None or cell.number is None for row in grid for cell in row)

    def test_cell_annotations(self):
        grid = grid_from_mask(["AB", "CD"])
        derive_words(grid)
        assert grid[0][0].number == 1
        assert grid[0][1].number == 2
        assert grid[1][0].number == 3
        assert grid[1][1].number is None
        assert (grid[1][1].across_num, grid[1][1].down_num) == (3, 2)


class TestEnumerationFromAnswer:
    @pytest.mark.parametrize(
        ("answer", "enum"),
        [
            ("FLOW", "(4)"),
            ("OMEGA POINT", "(5,5)"),
            ("X-RAY", "(1-3)"),
            ("A B-C", "(1,1-1)"),
            (" PADDED ", "(6)"),  # leading/trailing separators dropped
            ("A  B", "(1,1)"),  # adjacent separators collapse
            (" - ", "(0)"),  # all-separator answer collapses to 0 (engine rule)
        ],
    )
    def test_matches_engine_rule(self, answer, enum):
        assert enumeration_from_answer(answer) == enum


class TestExolveEnumParseBack:
    def test_clue_internal_parens_kept(self):
        assert split_clue_enum("Send (a letter) to HQ (4)") == ("Send (a letter) to HQ", "(4)")

    def test_no_enum(self):
        assert split_clue_enum("Ring (2 wds)") == ("Ring (2 wds)", None)

    def test_separator_variants(self):
        assert split_clue_enum("Multi (3,4,5)") == ("Multi", "(3,4,5)")
        assert split_clue_enum("Hyphen (4-5)") == ("Hyphen", "(4-5)")
        assert split_clue_enum("Apostrophe (1'4)") == ("Apostrophe", "(1'4)")
        assert split_clue_enum("Nested (a (b) c) (3)") == ("Nested (a (b) c)", "(3)")

    def test_cell_count(self):
        assert enumeration_cell_count("(5,5)") == 10
        assert enumeration_cell_count("(1-3)") == 4

    def test_cross_check_accepts_matching_enum(self):
        assert extract_clue_enum("Transcending entropy (5,5)", 10) == (
            "Transcending entropy",
            "(5,5)",
        )

    def test_cross_check_rejects_mismatch_and_derives(self):
        # trailing (5) is clue text on a 4-cell word: fold back, derive from grid
        assert extract_clue_enum("Beethoven's Fifth (5)", 4) == (
            "Beethoven's Fifth (5)",
            "(4)",
        )

    def test_no_enum_derives_from_cells(self):
        assert extract_clue_enum("Ring (2 wds)", 4) == ("Ring (2 wds)", "(4)")
