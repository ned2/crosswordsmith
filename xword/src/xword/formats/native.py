"""native canonical JSON <-> Board (docs/json-output-spec.md §6).

Numbering/word geometry are re-derived from the block pattern; the source's
`words[]` entries are attached to the derived words by (direction, cells).
"""

from __future__ import annotations

import json
from typing import Any

from .. import XwordError
from ..board import (
    Board,
    Cell,
    Grid,
    derive_words,
    enumeration_from_answer,
    grid_answer,
    validate_grid,
)
from ..util import dump_json


def parse(text: str) -> Board:
    try:
        obj = json.loads(text)
    except json.JSONDecodeError as e:
        raise XwordError(f"native: invalid JSON: {e}") from e
    if not isinstance(obj, dict):
        raise XwordError("native: top level must be a JSON object")
    for key in ("gridLength", "grid", "words"):
        if key not in obj:
            raise XwordError(f"native: missing required key {key!r}")
    n = obj["gridLength"]
    rows = obj["grid"]
    if not isinstance(n, int) or not isinstance(rows, list) or len(rows) != n:
        raise XwordError("native: grid must be gridLength rows")

    grid: Grid = []
    for r, row in enumerate(rows):
        if not isinstance(row, list) or len(row) != n:
            raise XwordError(f"native: grid row {r} must have gridLength cells")
        cells: list[Cell | None] = []
        for c, entry in enumerate(row):
            if entry is None:
                cells.append(None)
            elif isinstance(entry, dict) and isinstance(entry.get("letter"), str):
                cells.append(Cell(letter=entry["letter"]))
            else:
                raise XwordError(f"native: bad cell at [{r},{c}]: {entry!r}")
        grid.append(cells)
    validate_grid(grid)

    words = derive_words(grid)
    by_start = {(w.direction, w.cells[0]): w for w in words}
    for entry in obj["words"]:
        if not isinstance(entry, dict):
            raise XwordError(f"native: bad words[] entry: {entry!r}")
        direction = entry.get("direction")
        src_cells = [tuple(cell) for cell in entry.get("cells", [])]
        if not src_cells or direction not in ("across", "down"):
            raise XwordError(f"native: bad words[] entry: {entry!r}")
        word = by_start.get((direction, src_cells[0]))
        if word is None or word.cells != src_cells:
            raise XwordError(
                f"native: {direction} word at {src_cells[0]} does not match the "
                "grid's block pattern (words[] and grid disagree)"
            )
        answer = entry.get("answer")
        if not isinstance(answer, str):
            raise XwordError(f"native: word at {src_cells[0]} has no answer string")
        meta = entry.get("meta", {})
        if not isinstance(meta, dict):
            raise XwordError(f"native: word at {src_cells[0]} meta must be an object")
        word.answer = answer
        word.enumeration = enumeration_from_answer(answer)
        word.meta = meta
        clue = meta.get("clue")
        if isinstance(clue, str):
            word.clue = clue

    return Board(height=n, width=n, grid=grid, words=words)


def serialize(board: Board) -> str:
    if board.height != board.width:
        raise XwordError(
            f"native cannot hold a rectangular grid ({board.width}x{board.height}); "
            "native is square-only — target ipuz or exolve instead"
        )
    grid_out: list[list[Any]] = []
    for r, row in enumerate(board.grid):
        row_out: list[Any] = []
        for c, cell in enumerate(row):
            if cell is None:
                row_out.append(None)
                continue
            if cell.letter is None:
                raise XwordError(
                    f"native cannot hold an unfilled cell (at [{r},{c}]); "
                    "native carries solved layouts only"
                )
            if len(cell.letter) > 1:
                raise XwordError(
                    f"native cannot hold a rebus cell ({cell.letter!r} at [{r},{c}])"
                )
            row_out.append(
                {
                    "across": cell.across_num,
                    "down": cell.down_num,
                    "letter": cell.letter,
                    "number": cell.number,
                }
            )
        grid_out.append(row_out)

    words_out: list[dict[str, Any]] = []
    for word in sorted(board.words, key=lambda w: (w.number, w.direction)):
        answer = word.answer if word.answer is not None else grid_answer(board, word)
        if answer is None:
            raise XwordError(
                f"native cannot hold word {word.number} {word.direction}: unfilled cells"
            )
        meta = dict(word.meta)
        if word.clue is not None and "clue" not in meta:
            meta["clue"] = word.clue
        words_out.append(
            {
                "answer": answer,
                "cells": [list(cell) for cell in word.cells],
                "direction": word.direction,
                "meta": meta,
                "number": word.number,
            }
        )

    return dump_json({"grid": grid_out, "gridLength": board.height, "words": words_out})
