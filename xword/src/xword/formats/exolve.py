"""Exolve <-> Board (S5 findings).

Grid rows are TOKENISED (main char + trailing decorator suffixes), never
column-indexed — `0@BR|S` is 4 cells. An unfilled light cell is `0`. Exolve
stores no numbering (it re-derives from geometry, as we do). The enumeration
is embedded in the clue text; parse-back uses the normative end-anchored
regex with the cell-count cross-check (spec §5.1).
"""

from __future__ import annotations

import re

from .. import XwordError
from ..board import (
    ACROSS,
    DOWN,
    Board,
    Cell,
    Grid,
    enumeration_from_cells,
    derive_words,
    extract_clue_enum,
    validate_grid,
)

_DECORATORS = "|_@!*~"
_CLUE_LINE_RE = re.compile(r"^(\d+)\s+(.*)$")


def parse(text: str) -> Board:
    width: int | None = None
    height: int | None = None
    meta: dict[str, str] = {}
    sections: dict[str, list[str]] = {"grid": [], "across": [], "down": []}
    section: str | None = None

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("exolve-"):
            directive, _, value = line.partition(":")
            directive = directive[len("exolve-"):].strip()
            value = value.strip()
            section = None
            if directive == "width":
                width = _int_directive("width", value)
            elif directive == "height":
                height = _int_directive("height", value)
            elif directive == "title":
                meta["title"] = value
            elif directive == "setter":
                meta["author"] = value
            elif directive in sections:
                section = directive
            # other directives (begin/end/id/option/…) carry nothing Board models
            continue
        if section is not None:
            sections[section].append(line)

    if width is None or height is None:
        raise XwordError("exolve: missing exolve-width/exolve-height")
    rows = sections["grid"]
    if len(rows) != height:
        raise XwordError(f"exolve: expected {height} grid rows, got {len(rows)}")

    grid: Grid = []
    for r, row in enumerate(rows):
        tokens = _tokenize_row(row, r)
        if len(tokens) != width:
            raise XwordError(f"exolve: grid row {r} has {len(tokens)} cells, expected {width}")
        cells: list[Cell | None] = []
        for main, decorators in tokens:
            if main == ".":
                cells.append(None)
                continue
            letter = main if main.isalpha() else None  # '0'/'?' = unfilled light
            cells.append(
                Cell(
                    letter=letter,
                    circle="@" in decorators,
                    bar_right="|" in decorators,
                    bar_below="_" in decorators,
                    prefilled="!" in decorators,
                )
            )
        grid.append(cells)
    validate_grid(grid)

    words = derive_words(grid)
    board = Board(height=height, width=width, grid=grid, words=words, meta=dict(meta))
    _attach_clues(board, ACROSS, sections["across"])
    _attach_clues(board, DOWN, sections["down"])
    return board


def _int_directive(name: str, value: str) -> int:
    try:
        return int(value)
    except ValueError:
        raise XwordError(f"exolve: exolve-{name} must be an integer, got {value!r}") from None


def _tokenize_row(row: str, r: int) -> list[tuple[str, str]]:
    tokens: list[tuple[str, str]] = []
    i = 0
    while i < len(row):
        main = row[i]
        i += 1
        if main.isspace():
            continue
        decorators = ""
        while i < len(row):
            if row[i] in _DECORATORS:
                decorators += row[i]
                i += 1
            elif row[i] == "{":  # {k} SVG-decoration suffix
                end = row.find("}", i)
                if end == -1:
                    raise XwordError(f"exolve: unterminated {{...}} decorator in grid row {r}")
                i = end + 1
            else:
                break
        tokens.append((main, decorators))
    return tokens


def _attach_clues(board: Board, direction: str, lines: list[str]) -> None:
    by_number = {w.number: w for w in board.words if w.direction == direction}
    for line in lines:
        m = _CLUE_LINE_RE.match(line)
        if not m:
            raise XwordError(f"exolve: bad {direction} clue line (no leading number): {line!r}")
        number = int(m.group(1))
        word = by_number.get(number)
        if word is None:
            raise XwordError(
                f"exolve: clue {number} {direction} matches no derived word "
                "(numbering is re-derived from the grid)"
            )
        clue, enum = extract_clue_enum(m.group(2), len(word.cells))
        word.clue = clue
        word.enumeration = enum


def serialize(board: Board) -> str:
    lines = ["exolve-begin"]
    lines.append(f"  exolve-width: {board.width}")
    lines.append(f"  exolve-height: {board.height}")
    if "title" in board.meta:
        lines.append(f"  exolve-title: {board.meta['title']}")
    if "author" in board.meta:
        lines.append(f"  exolve-setter: {board.meta['author']}")
    lines.append("  exolve-grid:")
    for r, row in enumerate(board.grid):
        chars: list[str] = []
        for c, cell in enumerate(row):
            if cell is None:
                chars.append(".")
                continue
            if cell.letter is not None and len(cell.letter) > 1:
                raise XwordError(
                    f"exolve output does not support rebus cells yet "
                    f"({cell.letter!r} at [{r},{c}])"
                )
            chars.append(cell.letter or "0")
            if cell.circle:
                chars.append("@")
            if cell.prefilled:
                chars.append("!")
            if cell.bar_right:
                chars.append("|")
            if cell.bar_below:
                chars.append("_")
        lines.append("    " + "".join(chars))
    for header, direction in (("across", ACROSS), ("down", DOWN)):
        words = [w for w in board.words if w.direction == direction]
        if not any(w.clue is not None for w in words):
            continue
        lines.append(f"  exolve-{header}:")
        for word in sorted(words, key=lambda w: w.number):
            enum = word.enumeration or enumeration_from_cells(word.cells)
            parts = [str(word.number)] + ([word.clue] if word.clue else []) + [enum]
            lines.append("    " + " ".join(parts))
    lines.append("exolve-end")
    return "\n".join(lines) + "\n"
