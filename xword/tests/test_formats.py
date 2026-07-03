"""Parser round-trips (spec §11): identity where lossless, Board-equal otherwise."""

import json

import pytest

from xword import XwordError
from xword.board import ACROSS, DOWN
from xword.formats import parse_board, serialize_board


class TestNative:
    def test_round_trip_identity(self, fixtures):
        src = (fixtures / "bundled_17.native.json").read_text()
        out = serialize_board(parse_board(src, "native"), "native")
        assert json.loads(out) == json.loads(src)

    def test_serialize_idempotent(self, fixtures):
        src = (fixtures / "bundled_17.native.json").read_text()
        once = serialize_board(parse_board(src, "native"), "native")
        twice = serialize_board(parse_board(once, "native"), "native")
        assert once == twice

    def test_board_shape(self, fixtures):
        board = parse_board((fixtures / "bundled_17.native.json").read_text())
        assert (board.height, board.width) == (17, 17)
        assert board.meta == {}  # native carries no title/author
        omega = next(w for w in board.words if w.number == 1)
        assert (omega.direction, omega.answer) == (ACROSS, "OMEGA POINT")
        assert omega.enumeration == "(5,5)"  # derived from the answer's space
        assert omega.clue == "Transcending entropy"
        assert "link" in omega.meta  # metadata passthrough kept

    def test_crossing_cell_carries_both_words(self, fixtures):
        board = parse_board((fixtures / "bundled_17.native.json").read_text())
        cell = board.grid[0][3]  # OMEGA POINT x GNOSTIC GOSPELS
        assert (cell.across_num, cell.down_num, cell.number) == (1, 2, 2)

    def test_words_grid_disagreement_is_an_error(self, fixtures):
        obj = json.loads((fixtures / "bundled_17.native.json").read_text())
        obj["words"][0]["cells"] = [[0, 1], [0, 2]]  # not a derived word
        with pytest.raises(XwordError, match="does not match"):
            parse_board(json.dumps(obj), "native")


class TestIpuz:
    def test_round_trip_identity_vs_engine(self, fixtures):
        # our serializer reproduces the engine's ipuz export structurally
        src = (fixtures / "bundled_17.ipuz.json").read_text()
        out = serialize_board(parse_board(src, "ipuz"), "ipuz")
        assert json.loads(out) == json.loads(src)

    def test_engine_int_empty_marker_is_read_not_hardcoded(self, fixtures):
        # the engine declares empty: 0 (int); real-world files use "0" (str)
        board = parse_board((fixtures / "bundled_17.ipuz.json").read_text())
        assert board.grid[1][3].letter == "N"  # non-start white cell, int 0 in puzzle
        assert board.grid[1][0] is None  # "#" block

    def test_numbering_rederived_not_trusted(self, fixtures):
        board = parse_board((fixtures / "bundled_17.ipuz.json").read_text())
        got = {(w.number, w.direction) for w in board.words}
        assert got == {(1, ACROSS), (4, ACROSS), (5, ACROSS), (6, ACROSS), (2, DOWN), (3, DOWN)}

    def test_sample_board_round_trip(self, fixtures):
        src = (fixtures / "sample.ipuz.json").read_text()
        b1 = parse_board(src, "ipuz")
        b2 = parse_board(serialize_board(b1, "ipuz"), "ipuz")
        assert b1 == b2

    def test_sample_encodings(self, fixtures):
        board = parse_board((fixtures / "sample.ipuz.json").read_text())
        assert board.meta == {"title": "Probe Sample", "author": "S4 Spike"}
        assert board.grid[1][0].letter == "PH"  # rebus cell is multi-char
        assert board.grid[0][1].circle  # styled cell object {"cell", "style"}
        assert board.grid[2][2] is None  # str "#" block
        down_words = {w.number: w for w in board.words if w.direction == DOWN}
        assert down_words[1].clue == "First down"  # [number, "text"] clue shape
        assert down_words[1].enumeration == "(3)"  # derived: shape carries no enum
        across = {w.number: w for w in board.words if w.direction == ACROSS}
        assert across[1].enumeration == "(3)"  # object shape: taken from the field

    def test_rebus_blocks_native_serialization(self, fixtures):
        board = parse_board((fixtures / "sample.ipuz.json").read_text())
        # title and the circled cell block first (both structural too, D7);
        # strip them so the assertion pins the rebus check specifically
        board.meta.clear()
        for row in board.grid:
            for cell in row:
                if cell is not None:
                    cell.circle = False
                    cell.style = None
        with pytest.raises(XwordError, match="rebus"):
            serialize_board(board, "native")


class TestExolve:
    def test_board_round_trip(self, fixtures):
        src = (fixtures / "bundled_17.exolve").read_text()
        b1 = parse_board(src, "exolve")
        b2 = parse_board(serialize_board(b1, "exolve"), "exolve")
        assert b1 == b2

    def test_engine_export_parses(self, fixtures):
        board = parse_board((fixtures / "bundled_17.exolve").read_text())
        assert (board.height, board.width) == (17, 17)
        assert board.meta == {"title": "Untitled"}  # in the source, so kept
        omega = next(w for w in board.words if w.number == 1)
        assert omega.clue == "Transcending entropy"
        assert omega.enumeration == "(5,5)"  # parsed back from the clue line
        assert len(omega.cells) == 10

    def test_decorated_sample(self, fixtures):
        board = parse_board((fixtures / "sample_decorated.exolve").read_text())
        assert board.meta == {"title": "Decorated Sample", "author": "S5 Spike"}
        r0 = board.grid[0]
        assert r0[0].bar_right and r0[0].letter == "R"  # R|
        assert r0[1].bar_below  # E_
        assert board.grid[2][1].circle  # 0@
        assert board.grid[2][1].letter is None  # unfilled light = 0
        assert board.grid[1][1] is None  # . block
        assert board.grid[0][0].number == 1  # numbering from geometry, not stored

    def test_decorated_round_trip(self, fixtures):
        src = (fixtures / "sample_decorated.exolve").read_text()
        b1 = parse_board(src, "exolve")
        b2 = parse_board(serialize_board(b1, "exolve"), "exolve")
        assert b1 == b2

    def test_enum_collision_folds_into_clue(self, fixtures):
        # "Beethoven's Fifth (5)" on a 4-cell word: (5) is clue text, enum derived
        board = parse_board((fixtures / "sample_decorated.exolve").read_text())
        d1 = next(w for w in board.words if w.number == 1 and w.direction == DOWN)
        assert d1.clue == "Beethoven's Fifth (5)"
        assert d1.enumeration == "(4)"

    def test_unknown_clue_number_is_an_error(self, fixtures):
        src = (fixtures / "sample_decorated.exolve").read_text()
        with pytest.raises(XwordError, match="no derived word"):
            parse_board(src.replace("2 Adjoin (4)", "9 Adjoin (4)"), "exolve")
