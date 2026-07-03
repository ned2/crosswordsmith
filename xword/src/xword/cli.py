"""The `xword` CLI (spec §6, §8; cyclopts idioms per S1).

cyclopts renders parse/validation errors to stderr with exit 1 by itself;
only domain errors (XwordError) need the thin `main()` wrapper. The whole
output is built as one string and written once, so `--out` is written only
on success and never leaves a partial file.
"""

from __future__ import annotations

import sys
from enum import Enum
from pathlib import Path
from typing import Annotated, Optional

from cyclopts import App, Parameter
from cyclopts.types import ExistingFile
from rich.console import Console

from . import XwordError, __version__
from .convert import convert_text
from .formats import parse_board
from .render.geometry import board_geometry
from .render.terminal import print_view
from .util import capture_console, make_view_console

app = App(
    name="xword",
    help="Terminal viewer & format multitool for crossword layouts.",
    version=__version__,
)
_stderr = Console(stderr=True)


class DataFormat(Enum):
    native = "native"
    ipuz = "ipuz"
    exolve = "exolve"


class RenderFormat(Enum):
    svg = "svg"
    html = "html"
    png = "png"
    pdf = "pdf"


@app.default
def _bare() -> None:
    # A bare `xword` prints usage and exits non-zero (§8) — cyclopts' own
    # default prints help to stdout and exits 0, so override it.
    app.help_print(console=_stderr)
    raise SystemExit(2)


def _read_input(file: Optional[Path]) -> str:
    if file is not None:
        return file.read_text(encoding="utf-8")
    if sys.stdin.isatty():
        raise XwordError("no input: pass a FILE argument or pipe a layout to stdin")
    return sys.stdin.read()


@app.command
def view(
    file: Annotated[
        Optional[ExistingFile], Parameter(help="Layout file; omit to read stdin.")
    ] = None,
    *,
    from_: Annotated[
        Optional[DataFormat], Parameter(help="Input format; overrides detection.")
    ] = None,
    blank: Annotated[
        bool, Parameter(help="Render the empty reader's view (numbers, no letters).")
    ] = False,
    out: Annotated[
        Optional[Path], Parameter(name=("--out", "-o"), help="Write to FILE instead of stdout.")
    ] = None,
) -> None:
    """See a layout in the terminal: grid + Across/Down clue lists."""
    text = _read_input(file)
    board = parse_board(text, from_.value if from_ else None)
    geom = board_geometry(board, blank=blank)
    if out is not None:
        console = capture_console()
        print_view(console, geom)
        out.write_text(console.export_text(styles=False), encoding="utf-8")
    else:
        print_view(make_view_console(), geom)


@app.command
def convert(
    file: Annotated[
        Optional[ExistingFile], Parameter(help="Layout file; omit to read stdin.")
    ] = None,
    *,
    to: Annotated[DataFormat, Parameter(name="--to", help="Target format (required).")],
    from_: Annotated[
        Optional[DataFormat], Parameter(help="Input format; overrides detection.")
    ] = None,
    out: Annotated[
        Optional[Path], Parameter(name=("--out", "-o"), help="Write to FILE instead of stdout.")
    ] = None,
    quiet: Annotated[
        bool, Parameter(name=("--quiet", "-q"), help="Silence metadata-drop warnings.")
    ] = False,
) -> None:
    """Convert between data formats (any -> any)."""
    text = _read_input(file)
    result, warnings = convert_text(text, to.value, from_.value if from_ else None)
    if not quiet:
        for warning in warnings:
            print(f"xword: warning: {warning}", file=sys.stderr)
    if out is not None:
        out.write_text(result, encoding="utf-8")
    else:
        sys.stdout.write(result)


@app.command
def render(
    file: Annotated[
        Optional[ExistingFile], Parameter(help="Layout file; omit to read stdin.")
    ] = None,
    *,
    to: Annotated[RenderFormat, Parameter(name="--to", help="Target format (required).")],
    blank: Annotated[
        bool, Parameter(help="Render the empty reader's view (numbers, no letters).")
    ] = False,
    out: Annotated[
        Optional[Path], Parameter(name=("--out", "-o"), help="Write to FILE instead of stdout.")
    ] = None,
) -> None:
    """Render a layout to a visual artefact (svg/html/png/pdf). [not yet implemented]"""
    raise XwordError("'render' is not implemented yet (Phase 3; spec §13)")


def main() -> None:
    try:
        app()
    except XwordError as e:
        print(f"xword: error: {e}", file=sys.stderr)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
