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
    Meta,
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
    rebus = False

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
            elif directive == "option" and "rebus-cells" in value.split():
                rebus = True  # grid rows are space-separated multi-char tokens
            elif directive in sections:
                section = directive
            # other directives (begin/end/id/colour/…) carry nothing Board models
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
        tokens = _tokenize_rebus_row(row, r) if rebus else _tokenize_row(row, r)
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


def _tokenize_rebus_row(row: str, r: int) -> list[tuple[str, str]]:
    """Split a rebus-mode grid row: cells are whitespace-separated, each a
    multi-char main (letters / `0` / `.`) followed by any decorator suffixes."""
    tokens: list[tuple[str, str]] = []
    for chunk in row.split():
        i = 0
        while i < len(chunk) and chunk[i] not in _DECORATORS and chunk[i] != "{":
            i += 1
        main, rest = chunk[:i], chunk[i:]
        if not main:
            raise XwordError(f"exolve: rebus cell token {chunk!r} has no content in grid row {r}")
        decorators = ""
        j = 0
        while j < len(rest):
            if rest[j] in _DECORATORS:
                decorators += rest[j]
                j += 1
            elif rest[j] == "{":  # {k} SVG-decoration suffix
                end = rest.find("}", j)
                if end == -1:
                    raise XwordError(f"exolve: unterminated {{...}} decorator in grid row {r}")
                j = end + 1
            else:
                raise XwordError(f"exolve: bad rebus cell token {chunk!r} in grid row {r}")
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


# The canonical ipuz StyleSpec background-shade key. Exolve's colour model is
# one background colour + a cell list (`exolve-colour`), strictly coarser than
# the full StyleSpec, so only this background-fill key maps; text/border colour
# (`colortext`, `colorborder`), `imagebg`, `mark`, boolean `highlight`, … have
# no Exolve home and stay fail-strict.
_SHADE_KEY = "color"


def _style_expressible(style: Meta) -> bool:
    """True when every ipuz style key has an Exolve home.

    `shapebg: circle` and `barred` ride on the Cell's circle/bar decorator
    flags; the background shade `color` emits a top-level `exolve-colour`
    directive. Any other key genuinely has no Exolve home — fail-strict.
    """
    for key, value in style.items():
        if key == "shapebg" and value == "circle":
            continue
        if key == "barred":
            continue
        if key == _SHADE_KEY and isinstance(value, str) and value:
            continue
        return False
    return True


def _shade_colour(style: Meta | None) -> str | None:
    """The Exolve background colour for a cell, or None when it has no shade."""
    if not style:
        return None
    colour = style.get(_SHADE_KEY)
    return colour if isinstance(colour, str) and colour else None


# Note: no exolve-id is emitted (invent-nothing). Exolve treats the directive
# as optional and derives an id from a signature of the grid and clues when
# absent, so consumers lose nothing but cross-edit state continuity (spec §6.2).
# The title is also emitted only when present — this serializer never invents
# one; the Q5 default (`exolve-title: Untitled` for title-less boards, an Exet
# crash workaround) is injected at the `convert` boundary with a warning, so
# library round-trips through serialize() stay faithful (convert.py, spec §14).
def serialize(board: Board) -> str:
    lines = ["exolve-begin"]
    lines.append(f"  exolve-width: {board.width}")
    lines.append(f"  exolve-height: {board.height}")
    # A multi-char (rebus) cell flips the whole grid into Exolve's rebus-cells
    # mode: every row becomes space-separated tokens (§4.1). With no rebus cell
    # the compact, non-spaced grid is emitted byte-for-byte as before (D6).
    has_rebus = any(
        cell is not None and cell.letter is not None and len(cell.letter) > 1
        for row in board.grid
        for cell in row
    )
    if has_rebus:
        lines.append("  exolve-option: rebus-cells")
    if "title" in board.meta:
        lines.append(f"  exolve-title: {board.meta['title']}")
    if "author" in board.meta:
        lines.append(f"  exolve-setter: {board.meta['author']}")
    shades: list[tuple[str, int, int]] = []  # (colour, row, col), 0-indexed
    lines.append("  exolve-grid:")
    for r, row in enumerate(board.grid):
        tokens: list[str] = []
        for c, cell in enumerate(row):
            if cell is None:
                tokens.append(".")
                continue
            if cell.style and not _style_expressible(cell.style):
                raise XwordError(
                    f"exolve cannot hold ipuz cell styling ({cell.style!r} at "
                    f"[{r},{c}]) — target ipuz instead"
                )
            colour = _shade_colour(cell.style)
            if colour is not None:
                shades.append((colour, r, c))
            token = cell.letter or "0"  # multi-char = rebus, emitted verbatim
            if cell.circle:
                token += "@"
            if cell.prefilled:
                token += "!"
            if cell.bar_right:
                token += "|"
            if cell.bar_below:
                token += "_"
            tokens.append(token)
        lines.append("    " + (" " if has_rebus else "").join(tokens))
    # Background shades are top-level exolve-colour directives (not grid
    # decorators): one line per colour, cells grouped, deterministically sorted
    # by (colour, row, col) with 1-indexed rNcM cell-specs. Emitted only when a
    # shaded cell exists, so the no-shade path stays byte-identical (D6).
    by_colour: dict[str, list[tuple[int, int]]] = {}
    for colour, r, c in shades:
        by_colour.setdefault(colour, []).append((r, c))
    for colour in sorted(by_colour):
        specs = " ".join(f"r{r + 1}c{c + 1}" for r, c in sorted(by_colour[colour]))
        lines.append(f"  exolve-colour: {colour} {specs}")
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
