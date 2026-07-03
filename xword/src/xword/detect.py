"""Format detection ladder (spec §7). `--from` always overrides."""

from __future__ import annotations

import json

from . import XwordError

NATIVE = "native"
IPUZ = "ipuz"
EXOLVE = "exolve"

FORMATS = (NATIVE, IPUZ, EXOLVE)


def _is_ipuz_url(value: object) -> bool:
    return isinstance(value, str) and "ipuz.org" in value


def detect(text: str) -> str:
    """Return the format name for `text`, or raise a clear XwordError."""
    try:
        obj = json.loads(text)
    except (json.JSONDecodeError, UnicodeDecodeError):
        obj = None
    if isinstance(obj, dict):
        # ipuz first (S4): the signal is a version/kind VALUE containing the
        # ipuz.org host, not the mere presence of a `version` key.
        kind = obj.get("kind")
        if _is_ipuz_url(obj.get("version")) or (
            isinstance(kind, list) and any(_is_ipuz_url(k) for k in kind)
        ):
            return IPUZ
        if {"grid", "words", "gridLength"} <= obj.keys():
            return NATIVE
    if any(
        stripped.startswith(("exolve-begin", "exolve-width"))
        for stripped in (line.strip() for line in text.splitlines())
    ):
        return EXOLVE
    raise XwordError(
        "unrecognised input: not ipuz (no ipuz.org version/kind), not native "
        "(no grid+words+gridLength), not Exolve (no exolve-begin/exolve-width "
        "marker); use --from to name the format explicitly"
    )
