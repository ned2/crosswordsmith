"""`convert` — parse → D7 fidelity → serialize (spec §6.2, §10, D7).

Structural fidelity lives in the serializers: each raises an XwordError
naming the unsupported property and suggesting a capable target. This module
adds D7's metadata side — native-only annotations (per-word `link`/arbitrary
`meta` keys, top-level `diagnostics`) that ipuz/Exolve have no home for are
dropped by those serializers, and each drop is reported as a warning line the
CLI prints to stderr (unless `-q`). The native target keeps all of it, so
native → native is payload-lossless and warns about nothing.
"""

from __future__ import annotations

from .board import Board
from .detect import NATIVE
from .formats import parse_board, serialize_board


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


def convert_text(text: str, to_fmt: str, from_fmt: str | None = None) -> tuple[str, list[str]]:
    """Convert `text` to `to_fmt`, returning (output, drop warnings)."""
    board = parse_board(text, from_fmt)
    return serialize_board(board, to_fmt), metadata_drop_warnings(board, to_fmt)
