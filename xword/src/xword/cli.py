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
from .diff import diff_texts
from .formats import parse_board
from .render import raster
from .render.geometry import board_geometry
from .render.html import html_document
from .render.svg import master_svg
from .render.terminal import print_view
from .stats import stats_text
from .util import capture_console, make_view_console

app = App(
    name="xword",
    help="Viewer, format, rendering, and analysis multitool for crossword layouts.",
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
        try:
            return file.read_text(encoding="utf-8")
        except (OSError, UnicodeError) as e:
            detail = e.strerror if isinstance(e, OSError) else str(e)
            raise XwordError(f"cannot read {str(file)!r}: {detail or e}") from e
    if sys.stdin.isatty():
        raise XwordError("no input: pass a FILE argument or pipe a layout to stdin")
    return _read_stdin()


def _read_stdin(*, prefix: str = "") -> str:
    try:
        return sys.stdin.read()
    except (OSError, UnicodeError) as e:
        detail = e.strerror if isinstance(e, OSError) else str(e)
        raise XwordError(f"{prefix}cannot read stdin: {detail or e}") from e


def _write_text(text: str, out: Optional[Path]) -> None:
    """Write a text artefact to --out (built fully first, so only on success) or stdout."""
    try:
        if out is not None:
            out.write_text(text, encoding="utf-8")
        else:
            sys.stdout.write(text)
    except (OSError, UnicodeError) as e:
        detail = e.strerror if isinstance(e, OSError) else str(e)
        target = str(out) if out is not None else "stdout"
        raise XwordError(f"cannot write {target!r}: {detail or e}") from e


def _read_diff_inputs(left: str, right: str) -> tuple[str, str]:
    if left == "-" and right == "-":
        raise XwordError("diff: only one operand may be '-' (stdin)")

    def read(operand: str) -> str:
        if operand == "-":
            if sys.stdin.isatty():
                raise XwordError("diff: stdin operand '-' has no piped input")
            return _read_stdin(prefix="diff: ")
        try:
            return Path(operand).read_text(encoding="utf-8")
        except (OSError, UnicodeError) as e:
            detail = e.strerror if isinstance(e, OSError) else str(e)
            raise XwordError(f"diff: cannot read {operand!r}: {detail or e}") from e

    return read(left), read(right)


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
        bool, Parameter(name=("--quiet", "-q"), help="Silence conversion warnings.")
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
    """Render a layout to a visual artefact (svg/html/png/pdf)."""
    text = _read_input(file)
    board = parse_board(text, from_.value if from_ else None)
    geom = board_geometry(board, blank=blank)

    if to is RenderFormat.svg:
        _write_text(master_svg(geom), out)
    elif to is RenderFormat.html:
        _write_text(html_document(geom), out)
    else:  # png / pdf — binary; refuse to spew bytes at a terminal (§6.3)
        if out is None and sys.stdout.isatty():
            raise XwordError(
                f"refusing to write binary {to.value} to a terminal; use --out FILE or redirect"
            )
        svg = master_svg(geom)
        data = raster.to_png(svg) if to is RenderFormat.png else raster.to_pdf(svg)
        if out is not None:
            out.write_bytes(data)
        else:
            sys.stdout.buffer.write(data)


@app.command
def stats(
    file: Annotated[
        Optional[ExistingFile], Parameter(help="Layout file; omit to read stdin.")
    ] = None,
    *,
    from_: Annotated[
        Optional[DataFormat], Parameter(help="Input format; overrides detection.")
    ] = None,
    json_: Annotated[
        bool, Parameter(name="--json", help="Emit deterministic JSON.")
    ] = False,
    out: Annotated[
        Optional[Path], Parameter(name=("--out", "-o"), help="Write to FILE instead of stdout.")
    ] = None,
) -> None:
    """Report deterministic descriptive layout statistics."""
    text = _read_input(file)
    _write_text(stats_text(text, from_.value if from_ else None, json_output=json_), out)


@app.command
def diff(
    left: Annotated[
        str, Parameter(help="First layout file, or '-' for stdin.", allow_leading_hyphen=True)
    ],
    right: Annotated[
        str, Parameter(help="Second layout file, or '-' for stdin.", allow_leading_hyphen=True)
    ],
    *,
    from_left: Annotated[
        Optional[DataFormat], Parameter(name="--from-left", help="First input format override.")
    ] = None,
    from_right: Annotated[
        Optional[DataFormat], Parameter(name="--from-right", help="Second input format override.")
    ] = None,
    json_: Annotated[
        bool, Parameter(name="--json", help="Emit deterministic JSON.")
    ] = False,
    out: Annotated[
        Optional[Path], Parameter(name=("--out", "-o"), help="Write to FILE instead of stdout.")
    ] = None,
) -> None:
    """Compare two layouts by normalised puzzle semantics."""
    try:
        left_text, right_text = _read_diff_inputs(left, right)
        result, different = diff_texts(
            left_text,
            right_text,
            from_left.value if from_left else None,
            from_right.value if from_right else None,
            json_output=json_,
        )
        _write_text(result, out)
    except XwordError as e:
        print(f"xword: error: {e}", file=sys.stderr)
        raise SystemExit(2) from None
    if different:
        raise SystemExit(1)


def main() -> None:
    try:
        app()
    except XwordError as e:
        print(f"xword: error: {e}", file=sys.stderr)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
