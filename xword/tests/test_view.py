"""`view` golden snapshots via the S3 capture recipe (byte-stable, env-invariant)."""

from xword.formats import parse_board
from xword.render.geometry import board_geometry
from xword.render.terminal import print_view, render_grid
from xword.util import capture_console


def capture_view(text: str, *, blank: bool = False) -> str:
    console = capture_console()
    print_view(console, board_geometry(parse_board(text), blank=blank))
    return console.export_text(styles=False)


def test_solved_golden(fixtures, golden):
    src = (fixtures / "bundled_17.native.json").read_text()
    assert capture_view(src) == (golden / "view_bundled_17_solved.txt").read_text()


def test_blank_golden(fixtures, golden):
    src = (fixtures / "bundled_17.native.json").read_text()
    assert capture_view(src, blank=True) == (golden / "view_bundled_17_blank.txt").read_text()


def test_capture_is_stable_across_runs(fixtures):
    src = (fixtures / "bundled_17.native.json").read_text()
    assert capture_view(src) == capture_view(src)


def test_same_view_for_all_formats_modulo_title(fixtures):
    # ipuz/exolve exports of the same layout render the identical grid + clues;
    # they differ only in the engine-invented "Untitled" title line up front
    def from_grid(text: str) -> str:
        return text[text.index("┌"):]

    native = capture_view((fixtures / "bundled_17.native.json").read_text())
    ipuz = capture_view((fixtures / "bundled_17.ipuz.json").read_text())
    exolve = capture_view((fixtures / "bundled_17.exolve").read_text())
    assert from_grid(native) == from_grid(ipuz) == from_grid(exolve)


def test_grid_never_wraps():
    # a grid wider than the console must overflow, not re-wrap (S2/§6.1)
    from xword.board import Board, Cell, derive_words

    grid = [[Cell(letter="A") for _ in range(30)] for _ in range(2)]
    words = derive_words(grid)
    board = Board(height=2, width=30, grid=grid, words=words)
    console = capture_console(width=40)
    console.print(render_grid(board_geometry(board)), soft_wrap=True)
    lines = console.export_text(styles=False).splitlines()
    # 2 rows -> 5 canvas lines; wrapping would multiply that
    assert len(lines) == 5
    assert all(len(line) == 30 * 4 + 1 for line in lines)
