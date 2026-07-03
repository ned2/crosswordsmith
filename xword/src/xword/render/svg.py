"""SVG renderer: a pure, deterministic `BoardGeom -> SVG` (spec §6.3, D4/D6).

Two entry points share one grid-drawing pass:

* ``grid_svg``   — a self-contained ``<svg>`` of *just* the grid, embedded
  inline by the HTML renderer.
* ``master_svg`` — the standalone master: title (if any) + grid + Across/Down
  clue lists laid out in two columns below it. This is the ``svg`` target and
  the source the raster stage converts to PNG/PDF.

Determinism (S3/D6): every element id is a pure function of grid coordinates
(``cell-r{r}-c{c}``, ``word-{n}-{across|down}``) — no uuid, counters, or
timestamps. All coordinates are integers, so serialisation is byte-stable.
Glyphs are ``<text>`` (Q3): host-font-dependent raster is acceptable, and
cairosvg embeds glyph outlines into the PDF regardless.
"""

from __future__ import annotations

from xml.sax.saxutils import escape

from .geometry import BoardGeom, ClueEntry, GeomCell

# --- Geometry constants (all integers -> byte-stable output) -----------------
CELL = 32  # cell side, px
MARGIN = 8  # outer margin so bars/strokes are not clipped
GAP = 20  # vertical gap between grid and clue block

BLOCK_FILL = "#111111"
LIGHT_FILL = "#ffffff"
GRID_LINE = "#333333"
NUMBER_FILL = "#0066cc"
LETTER_FILL = "#111111"
CIRCLE_STROKE = "#444444"
BAR_STROKE = "#000000"
BAR_WIDTH = 3

FONT_FAMILY = "'Helvetica Neue', Arial, sans-serif"
NUMBER_FONT = 9
LETTER_FONT = 18
TITLE_FONT = 20
HEADING_FONT = 14
CLUE_FONT = 13
LINE_H = 18

# Clue block layout.
MIN_CONTENT = 360  # min two-column content width, px
GUTTER = 24  # gap between the Across and Down columns
NUM_COL = 26  # right-aligned number gutter inside a column
NUM_GAP = 7  # gap between the number and the clue text, px
CLUE_ADVANCE = 7  # px per char at CLUE_FONT — fixed => deterministic wrapping


def _esc(text: str) -> str:
    return escape(text)


def _rebus_font(letter: str) -> int:
    """Shrink a multi-char (rebus) glyph so it fits the cell width."""
    if len(letter) <= 1:
        return LETTER_FONT
    return max(8, (CELL - 6) * 2 // (len(letter) + 1))


def _cell_elements(cell: GeomCell, r: int, c: int, ox: int, oy: int) -> list[str]:
    """SVG for one cell: rect + number + letter + circle/bar decorations."""
    x = ox + c * CELL
    y = oy + r * CELL
    fill = BLOCK_FILL if cell.block else LIGHT_FILL
    els = [
        f'<rect id="cell-r{r}-c{c}" x="{x}" y="{y}" width="{CELL}" height="{CELL}" '
        f'fill="{fill}" stroke="{GRID_LINE}" stroke-width="1"/>'
    ]
    if cell.block:
        return els
    if cell.circle:
        els.append(
            f'<circle cx="{x + CELL // 2}" cy="{y + CELL // 2}" r="{CELL // 2 - 3}" '
            f'fill="none" stroke="{CIRCLE_STROKE}" stroke-width="1"/>'
        )
    if cell.number is not None:
        els.append(
            f'<text x="{x + 2}" y="{y + NUMBER_FONT + 1}" '
            f'font-size="{NUMBER_FONT}" fill="{NUMBER_FILL}">{cell.number}</text>'
        )
    if cell.letter:
        size = _rebus_font(cell.letter)
        # Baseline nudged below centre (~0.35em) so the glyph reads centred
        # without relying on dominant-baseline (uneven cairosvg support).
        els.append(
            f'<text x="{x + CELL // 2}" y="{y + CELL // 2 + size * 35 // 100}" '
            f'font-size="{size}" fill="{LETTER_FILL}" text-anchor="middle">'
            f"{_esc(cell.letter)}</text>"
        )
    if cell.bar_right:
        els.append(
            f'<line x1="{x + CELL}" y1="{y}" x2="{x + CELL}" y2="{y + CELL}" '
            f'stroke="{BAR_STROKE}" stroke-width="{BAR_WIDTH}"/>'
        )
    if cell.bar_below:
        els.append(
            f'<line x1="{x}" y1="{y + CELL}" x2="{x + CELL}" y2="{y + CELL}" '
            f'stroke="{BAR_STROKE}" stroke-width="{BAR_WIDTH}"/>'
        )
    return els


def _grid_elements(geom: BoardGeom, ox: int, oy: int) -> list[str]:
    els: list[str] = []
    for r in range(geom.rows):
        for c in range(geom.cols):
            els += _cell_elements(geom.cells[r][c], r, c, ox, oy)
    return els


def _wrap(text: str, budget: int) -> list[str]:
    """Greedy word-wrap to `budget` chars; an over-long word gets its own line."""
    lines: list[str] = []
    cur = ""
    for word in text.split():
        if not cur:
            cur = word
        elif len(cur) + 1 + len(word) <= budget:
            cur += " " + word
        else:
            lines.append(cur)
            cur = word
    lines.append(cur)  # always emit one line (empty clue -> one empty line)
    return lines


def _clue_column(
    heading: str, entries: list[ClueEntry], direction: str, x: int, y: int, width: int
) -> tuple[list[str], int]:
    """Lay out one clue column; return (elements, consumed height in px)."""
    els = [
        f'<text x="{x}" y="{y + HEADING_FONT}" font-size="{HEADING_FONT}" '
        f'font-weight="bold" fill="{LETTER_FILL}">{_esc(heading)}</text>'
    ]
    # Advance a full line below the heading baseline so the first clue clears
    # the heading (matching the inter-clue baseline step); +8 overlapped it.
    cy = y + HEADING_FONT + LINE_H
    budget = max(1, (width - NUM_COL) // CLUE_ADVANCE)
    text_x = x + NUM_COL
    for entry in entries:
        lines = _wrap(entry.text, budget)
        tspans = "".join(
            f'<tspan x="{text_x}" dy="{0 if i == 0 else LINE_H}">{_esc(line)}</tspan>'
            for i, line in enumerate(lines)
        )
        els.append(
            f'<g class="clue" id="word-{entry.number}-{direction}">'
            f'<text x="{x + NUM_COL - NUM_GAP}" y="{cy}" font-size="{CLUE_FONT}" '
            f'text-anchor="end" fill="{NUMBER_FILL}">{entry.number}</text>'
            f'<text y="{cy}" font-size="{CLUE_FONT}" fill="{LETTER_FILL}">{tspans}</text>'
            f"</g>"
        )
        cy += LINE_H * len(lines)
    return els, cy - y


def _svg(width: int, height: int, body: list[str], css_class: str | None = None) -> str:
    cls = f' class="{css_class}"' if css_class else ""
    head = (
        f'<svg xmlns="http://www.w3.org/2000/svg"{cls} '
        f'width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
        f'font-family="{FONT_FAMILY}">'
    )
    return head + "".join(body) + "</svg>"


def grid_svg(geom: BoardGeom, *, css_class: str | None = None) -> str:
    """A self-contained ``<svg>`` of just the grid (embedded inline by HTML).

    The blank/solved choice was already applied when the geometry was built
    (``board_geometry(board, blank=...)``), so the renderer just draws it.
    """
    width = geom.cols * CELL + 2 * MARGIN
    height = geom.rows * CELL + 2 * MARGIN
    body = _grid_elements(geom, MARGIN, MARGIN)
    return _svg(width, height, body, css_class)


def master_svg(geom: BoardGeom) -> str:
    """The standalone master: title (if any) + grid + two-column clue lists.

    Blank/solved is baked into the geometry upstream (``board_geometry``).
    """
    grid_w = geom.cols * CELL
    grid_h = geom.rows * CELL
    body: list[str] = []

    top = MARGIN
    if geom.title:
        body.append(
            f'<text x="{MARGIN}" y="{top + TITLE_FONT}" font-size="{TITLE_FONT}" '
            f'font-weight="bold" fill="{LETTER_FILL}">{_esc(geom.title)}</text>'
        )
        top += TITLE_FONT + 12

    body += _grid_elements(geom, MARGIN, top)
    bottom = top + grid_h

    content_w = max(grid_w, MIN_CONTENT)
    if geom.across or geom.down:
        col_w = (content_w - GUTTER) // 2
        clue_y = bottom + GAP
        left, lh = _clue_column("Across", geom.across, "across", MARGIN, clue_y, col_w)
        right, rh = _clue_column(
            "Down", geom.down, "down", MARGIN + col_w + GUTTER, clue_y, col_w
        )
        body += left + right
        bottom = clue_y + max(lh, rh)

    width = max(grid_w, content_w) + 2 * MARGIN
    height = bottom + MARGIN
    return _svg(width, height, body)
