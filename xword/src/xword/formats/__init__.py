"""Data plane: `parse(fmt) -> Board` and `Board -> serialize(fmt)` (spec §5.2)."""

from __future__ import annotations

from .. import XwordError
from ..board import Board
from ..detect import EXOLVE, FORMATS, IPUZ, NATIVE, detect
from . import exolve, ipuz, native

_PARSERS = {NATIVE: native.parse, IPUZ: ipuz.parse, EXOLVE: exolve.parse}
_SERIALIZERS = {NATIVE: native.serialize, IPUZ: ipuz.serialize, EXOLVE: exolve.serialize}


def parse_board(text: str, fmt: str | None = None) -> Board:
    """Parse `text` into a Board, detecting the format unless `fmt` is given."""
    if fmt is None:
        fmt = detect(text)
    if fmt not in _PARSERS:
        raise XwordError(f"unknown format {fmt!r}; expected one of {', '.join(FORMATS)}")
    return _PARSERS[fmt](text)


def serialize_board(board: Board, fmt: str) -> str:
    if fmt not in _SERIALIZERS:
        raise XwordError(f"unknown format {fmt!r}; expected one of {', '.join(FORMATS)}")
    return _SERIALIZERS[fmt](board)
