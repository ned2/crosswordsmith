"""Terminal renderer: a pure `BoardGeom -> Rich renderable` (S2).

Manual `Text` box-drawing, NOT a Rich `Table` (a Table always reflows to
console width). Cell interior is 3 chars wide x 1 row tall; the 4x2 char
pitch cancels the terminal's ~2:1 aspect so the grid reads square. Blocks
are solid `█` (an SGR reverse-video block would vanish in plain capture).
Corner numbers are superscript digits inline before the letter.
"""

from __future__ import annotations

from rich.console import Console
from rich.text import Text

from .geometry import BoardGeom, ClueEntry, GeomCell

CELL_W = 3
GRID_STYLE = "grey42"
NUMBER_STYLE = "bold cyan"
LETTER_STYLE = "white"
HEADER_STYLE = "bold"

_SUPERSCRIPT = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

_JUNCTIONS = {
    (False, True, False, True): "┌",
    (False, True, True, False): "┐",
    (True, False, False, True): "└",
    (True, False, True, False): "┘",
    (False, True, True, True): "┬",
    (True, False, True, True): "┴",
    (True, True, False, True): "├",
    (True, True, True, False): "┤",
    (True, True, True, True): "┼",
}


def render_grid(geom: BoardGeom) -> Text:
    rows, cols = geom.rows, geom.cols
    out = Text()

    def hborder(row_above: int, row_below: int) -> None:
        up = 0 <= row_above < rows
        down = 0 <= row_below < rows
        for c in range(cols):
            out.append(_JUNCTIONS[(up, down, c > 0, True)], style=GRID_STYLE)
            out.append("─" * CELL_W, style=GRID_STYLE)
        out.append(_JUNCTIONS[(up, down, True, False)], style=GRID_STYLE)
        out.append("\n")

    hborder(-1, 0)
    for r in range(rows):
        out.append("│", style=GRID_STYLE)
        for c in range(cols):
            out.append(_cell_text(geom.cells[r][c]))
            out.append("│", style=GRID_STYLE)
        out.append("\n")
        hborder(r, r + 1)
    if out.plain.endswith("\n"):
        out = out[: len(out.plain) - 1]
    return out


def _cell_text(cell: GeomCell) -> Text:
    if cell.block:
        return Text("█" * CELL_W, style=GRID_STYLE)
    letter = cell.letter or " "
    t = Text()
    if cell.number is not None:
        number = str(cell.number).translate(_SUPERSCRIPT)
        body = (number + letter)[:CELL_W].ljust(CELL_W)
        t.append(body[: len(number)], style=NUMBER_STYLE)
        t.append(body[len(number):], style=LETTER_STYLE if cell.letter else "")
    else:
        t.append(letter.center(CELL_W), style=LETTER_STYLE if cell.letter else "")
    return t


def render_clue_list(header: str, entries: list[ClueEntry]) -> Text:
    out = Text()
    out.append(header, style=HEADER_STYLE)
    num_width = max((len(str(e.number)) for e in entries), default=1)
    for entry in entries:
        out.append("\n")
        out.append(f"{entry.number:>{num_width}} ", style=NUMBER_STYLE)
        out.append(entry.text)
    return out


def print_view(console: Console, geom: BoardGeom) -> None:
    """Print the full `view` layout: title, grid (never wrapped), clue lists."""
    if geom.title:
        console.print(Text(geom.title, style=HEADER_STYLE))
        console.print()
    # Fixed geometry; if wider than the console it overflows, never re-wraps.
    console.print(render_grid(geom), soft_wrap=True)
    for header, entries in (("Across", geom.across), ("Down", geom.down)):
        if entries:
            console.print()
            console.print(render_clue_list(header, entries))
