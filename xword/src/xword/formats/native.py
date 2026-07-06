"""native canonical JSON <-> Board (docs/json-output-spec.md §6).

Numbering/word geometry are re-derived from the block pattern; the source's
`words[]` entries are attached to the derived words by (direction, cells).

The serializer is the structural-fidelity authority for the native target
(D7, spec §10): native is the narrowest box, so anything it cannot represent
— rectangular grids, circled/barred/prefilled/styled cells, rebus, unfilled
cells — raises an XwordError naming the property and a capable target. The
optional top-level title/author anchors (json-output-spec §6.5) are the
exception: native now carries them (parsed into / emitted from board.meta),
so title/author no longer fail-strict. Native's enumeration lives in the
answer's spaces/hyphens,
so an answer-less word (ipuz/Exolve source) is reconstructed from the grid
letters split per its enumeration; on enum/grid mismatch the grid wins
(same principle as the Exolve enum disambiguator, spec §5.1).
"""

from __future__ import annotations

import json
import re
from typing import Any

from .. import XwordError
from ..board import (
    Board,
    Cell,
    Grid,
    Word,
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

    diagnostics = obj.get("diagnostics")
    if diagnostics is not None and not isinstance(diagnostics, dict):
        raise XwordError("native: diagnostics must be an object")
    board = Board(height=n, width=n, grid=grid, words=words, diagnostics=diagnostics)
    # Optional top-level anchors (json-output-spec §6.5): read title/author into
    # board.meta, mirroring ipuz/exolve — the puzzle-level meta store all three
    # formats round-trip through. A non-string (or absent) value is simply not
    # carried; nothing is invented at parse time.
    for key in ("title", "author"):
        if isinstance(obj.get(key), str):
            board.meta[key] = obj[key]
    return board


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
            if cell.circle:
                raise XwordError(
                    f"native cannot hold a circled cell (at [{r},{c}]) — "
                    "target ipuz or exolve instead"
                )
            if cell.bar_right or cell.bar_below:
                raise XwordError(
                    f"native cannot hold a barred cell (at [{r},{c}]) — "
                    "target ipuz or exolve instead"
                )
            if cell.prefilled:
                raise XwordError(
                    f"native cannot hold a prefilled cell (at [{r},{c}]) — "
                    "target exolve instead"
                )
            if cell.style:
                raise XwordError(
                    f"native cannot hold cell styling ({cell.style!r} at [{r},{c}]) — "
                    "target ipuz instead"
                )
            if cell.letter is None:
                raise XwordError(
                    f"native cannot hold an unfilled cell (at [{r},{c}]); native "
                    "carries solved layouts only — target ipuz or exolve instead"
                )
            if len(cell.letter) > 1:
                raise XwordError(
                    f"native cannot hold a rebus cell ({cell.letter!r} at [{r},{c}]) — "
                    "target ipuz instead"
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
        answer = word.answer
        if answer is None:
            letters = grid_answer(board, word)
            if letters is None:  # unreachable after the grid scan; kept as a guard
                raise XwordError(
                    f"native cannot hold word {word.number} {word.direction}: unfilled cells"
                )
            answer = _answer_from_enumeration(letters, word)
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

    payload: dict[str, Any] = {
        "grid": grid_out,
        "gridLength": board.height,
        "words": words_out,
    }
    # Optional top-level anchors (json-output-spec §6.5): emit title/author when
    # the board carries them. native is now a title/author carrier (matching
    # ipuz/exolve), so native → native is anchor-lossless and a source that
    # gained a title on an exolve/ipuz hop carries it back. (dump_json sorts
    # keys, so placement here is cosmetic.)
    for key in ("title", "author"):
        if key in board.meta:
            payload[key] = board.meta[key]
    if board.diagnostics is not None:
        # native -> native keeps diagnostics: the target can hold it, so D7's
        # drop-when-homeless rule never triggers (spec §10)
        payload["diagnostics"] = board.diagnostics
    return dump_json(payload)


# native encodes the enumeration in the answer's display form: "," (word
# break) as a space, "-" as a hyphen. Other enum separators (. ' and space)
# have no native encoding and are a structural failure.
_ENUM_SEP_TO_CHAR = {",": " ", "-": "-"}


def _answer_from_enumeration(letters: str, word: Word) -> str:
    """Rebuild the display answer from bare grid letters + the enumeration.

    On enum/grid mismatch the grid is the authority (spec §5.1): the plain
    letters are returned and the enumeration re-derives as `(n)`.
    """
    if word.enumeration is None:
        return letters
    parts = re.findall(r"[0-9]+|[^0-9]", word.enumeration.strip().strip("()"))
    if sum(int(p) for p in parts if p.isdigit()) != len(letters):
        return letters
    out: list[str] = []
    i = 0
    for p in parts:
        if p.isdigit():
            out.append(letters[i : i + int(p)])
            i += int(p)
        else:
            sep = _ENUM_SEP_TO_CHAR.get(p)
            if sep is None:
                raise XwordError(
                    f"native cannot hold the enumeration {word.enumeration!r} of word "
                    f"{word.number} {word.direction}: its answers encode only word "
                    "breaks (,) and hyphens (-) — target ipuz or exolve instead"
                )
            out.append(sep)
    return "".join(out)
