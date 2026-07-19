"""ipuz v2 <-> Board (S4 findings).

`empty`/`block` are *declared* top-level fields with defaults "0"/"#" — read,
never hardcode (the engine emits int 0 for empty; real-world files use "0").
A styled/labelled cell is an object {"cell": …, "style": …}. Clue items come
in three shapes: {number, clue, enumeration}, [number, "text"], bare string.
Source numbers are advisory; numbering is re-derived from the block pattern.
"""

from __future__ import annotations

import json
from typing import Any

from .. import XwordError
from ..board import (
    ACROSS,
    DOWN,
    Board,
    Cell,
    Grid,
    enumeration_from_cells,
    derive_words,
    validate_grid,
)
from ..util import dump_json

VERSION = "http://ipuz.org/v2"
KIND = "http://ipuz.org/crossword#1"


def _matches(value: Any, marker: Any) -> bool:
    return value == marker or str(value) == str(marker)


def parse(text: str) -> Board:
    try:
        obj = json.loads(text)
    except json.JSONDecodeError as e:
        raise XwordError(f"ipuz: invalid JSON: {e}") from e
    if not isinstance(obj, dict):
        raise XwordError("ipuz: top level must be a JSON object")

    dims = obj.get("dimensions")
    if not isinstance(dims, dict) or "width" not in dims or "height" not in dims:
        raise XwordError("ipuz: missing dimensions.width/height")
    width, height = dims["width"], dims["height"]
    block = obj.get("block", "#")
    empty = obj.get("empty", "0")

    puzzle = obj.get("puzzle")
    if not isinstance(puzzle, list) or len(puzzle) != height:
        raise XwordError("ipuz: puzzle must be `height` rows")
    solution = obj.get("solution")
    if solution is not None and (not isinstance(solution, list) or len(solution) != height):
        raise XwordError("ipuz: solution must be `height` rows")

    grid: Grid = []
    for r, row in enumerate(puzzle):
        if not isinstance(row, list) or len(row) != width:
            raise XwordError(f"ipuz: puzzle row {r} must have `width` cells")
        cells: list[Cell | None] = []
        for c, pc in enumerate(row):
            style: dict[str, Any] | None = None
            label = pc
            value: Any = None
            if isinstance(pc, dict):
                label = pc.get("cell")
                value = pc.get("value")
                raw_style = pc.get("style")
                if isinstance(raw_style, dict):
                    style = raw_style
            # ipuz null is a void cell, not strictly a block; the Board does
            # not model voids, so both render/serialize as blocks here.
            if label is None or _matches(label, block):
                cells.append(None)
                continue
            circle = bool(style) and style.get("shapebg") == "circle"
            barred = style.get("barred") if style else None
            barred = barred.upper() if isinstance(barred, str) else ""
            # A `value` on the puzzle cell is a given/prefilled letter (§4.2);
            # the matching solution letter reconciles below (they agree).
            prefilled = isinstance(value, str) and not (
                _matches(value, block) or _matches(value, empty)
            )
            cells.append(
                Cell(
                    letter=value if prefilled else None,
                    circle=circle,
                    bar_right="R" in barred,
                    bar_below="B" in barred,
                    prefilled=prefilled,
                    style=style,
                )
            )
        grid.append(cells)
    validate_grid(grid)

    # Canonicalise left/top bar spellings onto the neighbouring cell's
    # right/below edge. Board stores each interior edge exactly once; ipuz
    # permits naming that edge from either adjacent cell.
    for r, row in enumerate(grid):
        for c, cell in enumerate(row):
            if cell is None or not cell.style:
                continue
            barred = cell.style.get("barred")
            if not isinstance(barred, str):
                continue
            barred = barred.upper()
            if "L" in barred and c > 0 and grid[r][c - 1] is not None:
                grid[r][c - 1].bar_right = True
            if "T" in barred and r > 0 and grid[r - 1][c] is not None:
                grid[r - 1][c].bar_below = True

    if solution is not None:
        for r, row in enumerate(solution):
            if not isinstance(row, list) or len(row) != width:
                raise XwordError(f"ipuz: solution row {r} must have `width` cells")
            for c, sc in enumerate(row):
                cell = grid[r][c]
                if cell is None:
                    continue
                value = sc.get("value") if isinstance(sc, dict) else sc
                if value is None or _matches(value, block) or _matches(value, empty):
                    continue
                if not isinstance(value, str):
                    raise XwordError(f"ipuz: bad solution cell at [{r},{c}]: {sc!r}")
                cell.letter = value  # multi-char = rebus

    words = derive_words(grid)
    board = Board(height=height, width=width, grid=grid, words=words)

    clues = obj.get("clues", {})
    if not isinstance(clues, dict):
        raise XwordError("ipuz: clues must be an object")
    for key, direction in (("Across", ACROSS), ("Down", DOWN)):
        items = clues.get(key, [])
        if not isinstance(items, list):
            raise XwordError(f"ipuz: clues.{key} must be a list")
        _attach_clues(board, direction, items)

    for key in ("title", "author"):
        if isinstance(obj.get(key), str):
            board.meta[key] = obj[key]
    return board


def _attach_clues(board: Board, direction: str, items: list[Any]) -> None:
    words = [w for w in board.words if w.direction == direction]
    by_number = {w.number: w for w in words}
    for i, item in enumerate(items):
        number: int | None = None
        clue: str = ""
        enum: str | None = None
        if isinstance(item, dict):
            n = item.get("number")
            number = n if isinstance(n, int) else None
            raw_clue = item.get("clue")
            clue = raw_clue if isinstance(raw_clue, str) else ""
            raw_enum = item.get("enumeration")
            enum = f"({raw_enum.strip('()')})" if isinstance(raw_enum, str) else None
        elif isinstance(item, list) and len(item) == 2 and isinstance(item[0], int):
            number, clue = item[0], str(item[1])
        elif isinstance(item, str):
            clue = item  # bare string: positional — the i-th word of this direction
        else:
            raise XwordError(f"ipuz: unrecognised clues.{direction} item: {item!r}")

        if number is not None:
            word = by_number.get(number)
            if word is None:
                raise XwordError(
                    f"ipuz: clue {number} {direction} matches no derived word "
                    "(numbering is re-derived from the block pattern)"
                )
        elif i < len(words):
            word = words[i]
        else:
            raise XwordError(f"ipuz: more {direction} clues than {direction} words")
        word.clue = clue
        word.enumeration = enum if enum is not None else enumeration_from_cells(word.cells)


def _synth_style(cell: Cell) -> dict[str, Any] | None:
    """An ipuz StyleSpec for decorations carried as Cell flags (Exolve sources)."""
    style: dict[str, Any] = {}
    if cell.circle:
        style["shapebg"] = "circle"
    barred = ("R" if cell.bar_right else "") + ("B" if cell.bar_below else "")
    if barred:
        style["barred"] = barred
    return style or None


def serialize(board: Board) -> str:
    puzzle: list[list[Any]] = []
    solution: list[list[Any]] = []
    any_letter = False
    for r, row in enumerate(board.grid):
        prow: list[Any] = []
        srow: list[Any] = []
        for c, cell in enumerate(row):
            if cell is None:
                prow.append("#")
                srow.append("#")
                continue
            label: Any = cell.number if cell.number is not None else 0
            style = cell.style if cell.style is not None else _synth_style(cell)
            # A prefilled (given) cell shows its letter in the puzzle grid: ipuz
            # carries it as `value` on the puzzle cell object, alongside the
            # solution letter (dump_json sorts keys, so the object stays canonical).
            value = cell.letter if cell.prefilled else None
            if style or value is not None:
                pc: dict[str, Any] = {"cell": label}
                if value is not None:
                    pc["value"] = value
                if style:
                    pc["style"] = style
                prow.append(pc)
            else:
                prow.append(label)
            if cell.letter is None:
                srow.append(0)
            else:
                any_letter = True
                srow.append(cell.letter)
        puzzle.append(prow)
        solution.append(srow)

    obj: dict[str, Any] = {
        "version": VERSION,
        "kind": [KIND],
        "dimensions": {"width": board.width, "height": board.height},
        "empty": 0,
        "puzzle": puzzle,
    }
    if any_letter:
        obj["solution"] = solution
    for key in ("title", "author"):
        if key in board.meta:
            obj[key] = board.meta[key]

    if any(w.clue is not None for w in board.words):
        across: list[dict[str, Any]] = []
        down: list[dict[str, Any]] = []
        for word in sorted(board.words, key=lambda w: w.number):
            enum = word.enumeration or enumeration_from_cells(word.cells)
            item = {"number": word.number, "clue": word.clue or "", "enumeration": enum}
            (across if word.direction == ACROSS else down).append(item)
        obj["clues"] = {"Across": across, "Down": down}

    return dump_json(obj)
