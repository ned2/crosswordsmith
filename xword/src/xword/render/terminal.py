"""Terminal renderer: a pure `BoardGeom -> Rich renderable` (S2).

Manual `Text` box-drawing, NOT a Rich `Table` (a Table always reflows to
console width). Cell interior is 3 chars wide x 1 row tall; the 4x2 char
pitch cancels the terminal's ~2:1 aspect so the grid reads square. Blocks
are solid `█` (an SGR reverse-video block would vanish in plain capture).
Corner numbers are superscript digits inline before the letter.

Layout mirrors a printed broadsheet (the Guardian): the grid sits top-left
and the clues fill a borderless pane to its right, Across and Down as two
side-by-side columns. Everything is composed line-by-line into one fixed
`Text`, so — like the grid — it overflows a narrow console rather than
reflowing.
"""

from __future__ import annotations

import textwrap

from rich.console import Console
from rich.text import Text

from .geometry import BoardGeom, ClueEntry, GeomCell

CELL_W = 3
GRID_STYLE = "grey42"
NUMBER_STYLE = "bold cyan"
LETTER_STYLE = "white"
HEADER_STYLE = "bold"

# Right-pane geometry. CLUE_TEXT_WIDTH is the wrap width of the clue text
# alone; the number sits in a gutter to its left and wrapped lines hang under
# the text (never under the number), as a printed puzzle sets them.
CLUE_TEXT_WIDTH = 32
GUTTER_GAP = 2  # spaces between a clue number and its text
COLUMN_GAP = 4  # blank columns between the grid and each clue column

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


def _clue_lines(entry: ClueEntry, num_width: int) -> list[Text]:
    """A single clue as styled lines: number in the gutter, text hanging under itself."""
    gutter = f"{entry.number:>{num_width}}{' ' * GUTTER_GAP}"
    indent = " " * len(gutter)
    segments = textwrap.wrap(entry.text, width=CLUE_TEXT_WIDTH, break_on_hyphens=False) or [""]
    lines: list[Text] = []
    for i, segment in enumerate(segments):
        line = Text()
        if i == 0:
            line.append(gutter, style=NUMBER_STYLE)
        else:
            line.append(indent)
        line.append(segment)
        lines.append(line)
    return lines


def render_clue_column(header: str, entries: list[ClueEntry]) -> tuple[list[Text], int]:
    """Lay out one clue column (header + hanging-indent clues); returns lines + its width."""
    num_width = max((len(str(e.number)) for e in entries), default=1)
    lines: list[Text] = [Text(header, style=HEADER_STYLE)]
    for entry in entries:
        lines.append(Text(""))  # a blank line between clues, as a printed puzzle spaces them
        lines.extend(_clue_lines(entry, num_width))
    return lines, num_width + GUTTER_GAP + CLUE_TEXT_WIDTH


def _compose(grid: Text, columns: list[tuple[list[Text], int]]) -> Text:
    """Stitch the grid and each clue column side by side into one fixed-width Text."""
    grid_lines = grid.split("\n")
    grid_width = len(grid_lines[0].plain) if len(grid_lines) else 0
    height = max([len(grid_lines), *(len(lines) for lines, _ in columns)])
    out = Text()
    for i in range(height):
        line = Text()
        cell = grid_lines[i] if i < len(grid_lines) else None
        line.append(cell if cell is not None else Text())
        line.append(" " * (grid_width - (len(cell.plain) if cell is not None else 0)))
        for lines, width in columns:
            line.append(" " * COLUMN_GAP)
            cell = lines[i] if i < len(lines) else None
            line.append(cell if cell is not None else Text())
            line.append(" " * (width - (len(cell.plain) if cell is not None else 0)))
        line.rstrip()  # drop the padding that runs off the final column
        if i:
            out.append("\n")
        out.append(line)
    return out


def print_view(console: Console, geom: BoardGeom) -> None:
    """Print the full `view` layout: title, then the grid with a clue pane to its right."""
    if geom.title:
        console.print(Text(geom.title, style=HEADER_STYLE))
        console.print()
    columns = [
        render_clue_column(header, entries)
        for header, entries in (("Across", geom.across), ("Down", geom.down))
        if entries
    ]
    # Fixed geometry; if wider than the console it overflows, never re-wraps.
    console.print(_compose(render_grid(geom), columns), soft_wrap=True)
