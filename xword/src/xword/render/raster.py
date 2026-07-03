"""Raster backend: master SVG -> PNG/PDF via cairosvg (spec §6.3, §12; S6).

cairosvg is the optional ``xword[raster]`` extra; it is imported lazily so the
base install (and every non-raster command) never touches it. Its output is
byte-deterministic on a machine (S6) but host-font dependent across machines,
so raster artefacts are tested by dimensions/magic, not byte-golden (§11).

The two failure modes are reported distinctly (§8, S6):

* the extra is not installed        -> ``ImportError`` -> "install xword[raster]"
* the extra is installed but the
  system cairo library is missing    -> ``OSError`` (cairocffi failing to
  ``dlopen`` ``libcairo.so.2``)      -> "apt install libcairo2"
"""

from __future__ import annotations

from typing import Any

from .. import XwordError

_MISSING_EXTRA = (
    "PNG/PDF output needs the optional raster extra. Install it with "
    "`pip install 'xword[raster]'` (or run via `uv run --extra raster xword …`)."
)
_MISSING_CAIRO = (
    "the raster extra is installed but the system cairo library is missing "
    "(cairocffi could not load libcairo.so.2). Install it with "
    "`apt install libcairo2` (or your distro's equivalent)."
)


def _cairosvg() -> Any:
    """Import cairosvg, mapping its two failure modes to clear domain errors."""
    try:
        import cairosvg  # type: ignore[import-untyped]
    except ImportError as e:
        raise XwordError(_MISSING_EXTRA) from e
    except OSError as e:  # cairocffi dlopen of libcairo.so.2 failed
        raise XwordError(_MISSING_CAIRO) from e
    return cairosvg


def to_png(svg: str) -> bytes:
    return _cairosvg().svg2png(bytestring=svg.encode("utf-8"))


def to_pdf(svg: str) -> bytes:
    return _cairosvg().svg2pdf(bytestring=svg.encode("utf-8"))
