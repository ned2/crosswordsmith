"""`convert` — parse → D7 fidelity → serialize (spec §6.2, §10, D7).

Structural fidelity lives in the serializers: each raises an XwordError
naming the unsupported property and suggesting a capable target. This module
adds D7's metadata side — native-only annotations (per-word `link`/arbitrary
`meta` keys, top-level `diagnostics`) that ipuz/Exolve have no home for are
dropped by those serializers, and each drop is reported as a warning line the
CLI prints to stderr (unless `-q`). The native target keeps all of it, so
native → native is payload-lossless and warns about nothing.

This boundary also owns the ONE sanctioned addition (Q5, spec §14): a
title-less board targeting Exolve gains the engine's default title, because
real Exolve consumers require one (Exet's Save crashes on a null title —
docs/exet-verification.md). The serializers themselves stay invent-nothing.
"""

from __future__ import annotations

from .board import Board
from .detect import EXOLVE, NATIVE
from .formats import parse_board, serialize_board

# The engine's invented default (export.pl emits `exolve-title: Untitled`);
# byte-matching it keeps `convert` output line-identical to `crosswordsmith
# export` on the title line (spec §11 engine cross-check).
EXOLVE_DEFAULT_TITLE = "Untitled"


def metadata_drop_warnings(board: Board, to_fmt: str) -> list[str]:
    """The D7 drop-and-warn lines for serializing `board` to `to_fmt`."""
    if to_fmt == NATIVE:
        return []  # native holds word meta verbatim and diagnostics (§10)
    warnings: list[str] = []
    counts: dict[str, int] = {}
    for word in board.words:
        for key in word.meta:
            if key == "clue":
                continue  # has a home: the target's clue text
            counts[key] = counts.get(key, 0) + 1
    for key in sorted(counts):
        n = counts[key]
        warnings.append(
            f"{to_fmt} has no home for per-word meta key {key!r}; "
            f"dropped from {n} word{'s' if n != 1 else ''}"
        )
    if board.diagnostics is not None:
        warnings.append(f"{to_fmt} has no home for top-level 'diagnostics'; dropped")
    return warnings


def exolve_title_default(board: Board, to_fmt: str) -> list[str]:
    """Q5 (spec §14, resolved): inject the engine's default title on a
    title-less board targeting Exolve, and say so with a warning line.

    Exolve-only: Exet's Save crashes on a null title (exet-verification.md),
    the same evidence that made the engine's `export` inject its default.
    ipuz stays invent-nothing — `title` is optional there and nothing breaks.
    Mutates `board.meta` (parse_board built it for this conversion alone).
    """
    if to_fmt == EXOLVE and "title" not in board.meta:
        board.meta["title"] = EXOLVE_DEFAULT_TITLE
        return [
            "exolve requires a title (Exet's Save crashes without one); "
            f"emitted default {EXOLVE_DEFAULT_TITLE!r}"
        ]
    return []


def convert_text(text: str, to_fmt: str, from_fmt: str | None = None) -> tuple[str, list[str]]:
    """Convert `text` to `to_fmt`, returning (output, warning lines)."""
    board = parse_board(text, from_fmt)
    warnings = metadata_drop_warnings(board, to_fmt) + exolve_title_default(board, to_fmt)
    return serialize_board(board, to_fmt), warnings
