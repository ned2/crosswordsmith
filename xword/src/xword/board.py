"""The `Board` model — the union superset every command consumes (spec §5.1).

Numbering is *derived* by row-major scan, never trusted from the source
(ipuz numbers are advisory; Exolve stores none). Enumeration is *normalised*:
derived from the answer's spaces/hyphens (native), taken from the field
(ipuz), or parsed back from the clue text with a cell-count cross-check
(Exolve). "Invent nothing": a word with no clue carries no clue.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from . import XwordError

ACROSS = "across"
DOWN = "down"

Coord = tuple[int, int]  # (row, col), 0-indexed
Meta = dict[str, Any]


@dataclass
class Cell:
    """A light (non-block) cell. A block is represented as ``None`` in the grid."""

    letter: str | None = None  # multi-char = rebus
    number: int | None = None  # derived; non-null iff a word starts here
    across_num: int | None = None  # derived
    down_num: int | None = None  # derived
    circle: bool = False
    bar_right: bool = False
    bar_below: bool = False
    prefilled: bool = False
    style: Meta | None = None  # ipuz style passthrough (verbatim)


Grid = list[list[Cell | None]]  # row-major; None = block


@dataclass
class Word:
    number: int
    direction: str  # ACROSS | DOWN
    cells: list[Coord]
    answer: str | None = None  # display form; spaces/hyphens preserved
    enumeration: str | None = None  # normalised, e.g. "(5,5)" / "(4-5)"
    clue: str | None = None  # None = source carried no clue
    meta: Meta = field(default_factory=dict)  # native passthrough (clue, link, …)


@dataclass
class Board:
    height: int
    width: int
    grid: Grid
    words: list[Word] = field(default_factory=list)
    meta: Meta = field(default_factory=dict)  # puzzle-level: title?, author?


def derive_words(grid: Grid) -> list[Word]:
    """Standard row-major numbering over a block pattern.

    Mutates the cells' `number`/`across_num`/`down_num` in place and returns
    skeleton words (number, direction, cells) in scan order, across before
    down at a shared start cell. A word is a run of >= 2 light cells.
    """
    height = len(grid)
    width = len(grid[0]) if height else 0
    for row in grid:
        for cell in row:
            if cell is not None:
                cell.number = cell.across_num = cell.down_num = None

    def light(r: int, c: int) -> bool:
        return 0 <= r < height and 0 <= c < width and grid[r][c] is not None

    words: list[Word] = []
    n = 0
    for r in range(height):
        for c in range(width):
            start = grid[r][c]
            if start is None:
                continue
            starts_across = not light(r, c - 1) and light(r, c + 1)
            starts_down = not light(r - 1, c) and light(r + 1, c)
            if not (starts_across or starts_down):
                continue
            n += 1
            start.number = n
            if starts_across:
                cells: list[Coord] = []
                cc = c
                while light(r, cc):
                    cells.append((r, cc))
                    cell = grid[r][cc]
                    assert cell is not None
                    cell.across_num = n
                    cc += 1
                words.append(Word(number=n, direction=ACROSS, cells=cells))
            if starts_down:
                cells = []
                rr = r
                while light(rr, c):
                    cells.append((rr, c))
                    cell = grid[rr][c]
                    assert cell is not None
                    cell.down_num = n
                    rr += 1
                words.append(Word(number=n, direction=DOWN, cells=cells))
    return words


def grid_answer(board: Board, word: Word) -> str | None:
    """The word's letters read off the grid, or None if any cell is unfilled."""
    letters: list[str] = []
    for r, c in word.cells:
        cell = board.grid[r][c]
        if cell is None or cell.letter is None:
            return None
        letters.append(cell.letter)
    return "".join(letters)


# --- Enumeration -------------------------------------------------------------
# Derivation from the answer's display form mirrors the engine's export
# (space -> ",", hyphen -> "-"; zero-length runs and their separators are
# dropped; an all-separator answer collapses to "(0)").

_SEPARATORS = {" ": ",", "-": "-"}


def enumeration_from_answer(answer: str) -> str:
    segments: list[str] = []
    run = 0
    pending_sep: str | None = None
    for ch in answer:
        if ch in _SEPARATORS:
            if run > 0:
                if segments:
                    segments.append(pending_sep or "")
                segments.append(str(run))
                run = 0
            pending_sep = _SEPARATORS[ch]
        else:
            run += 1
    if run > 0:
        if segments:
            segments.append(pending_sep or "")
        segments.append(str(run))
    body = "".join(segments) if segments else "0"
    return f"({body})"


def enumeration_from_cells(cells: list[Coord]) -> str:
    return f"({len(cells)})"


# The normative Exolve parse-back rule (spec §5.1 / S5): the enumeration is
# the LAST parenthesised group anchored at end-of-string whose interior is
# only the enum alphabet — digits and the separators , - . ' and space.
_ENUM_RE = re.compile(r"\s*\(([0-9]+(?:[,\-. ']*[0-9]+)*)\)\s*$")


def split_clue_enum(text: str) -> tuple[str, str | None]:
    """Split ``"clue text (5,5)"`` into (clue, enumeration) — no cross-check."""
    text = text.rstrip()
    m = _ENUM_RE.search(text)
    if not m:
        return (text.strip(), None)
    return (text[: m.start()].rstrip(), f"({m.group(1)})")


def enumeration_cell_count(enumeration: str) -> int:
    return sum(int(n) for n in re.findall(r"[0-9]+", enumeration))


def extract_clue_enum(text: str, cell_count: int) -> tuple[str, str]:
    """Split a clue with the grid as disambiguator (spec §5.1).

    Accept a trailing parenthesised enum only if its integer runs sum to the
    word's cell count; otherwise the parens are clue text and the enumeration
    is derived from the cell count.
    """
    clue, enum = split_clue_enum(text)
    if enum is not None and enumeration_cell_count(enum) == cell_count:
        return (clue, enum)
    return (text.strip(), f"({cell_count})")


def validate_grid(grid: Grid) -> None:
    if not grid or not grid[0]:
        raise XwordError("empty grid")
    width = len(grid[0])
    for i, row in enumerate(grid):
        if len(row) != width:
            raise XwordError(f"ragged grid: row {i} has {len(row)} cells, expected {width}")
