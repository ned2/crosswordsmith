"""HTML renderer: a semantic page over the shared geometry (spec §6.3).

The grid is the inline SVG from ``svg.grid_svg``; the clue lists are semantic
HTML (``<ol>`` with per-clue ``value`` + coordinate-derived ids). Styling is a
single built-in inline ``<style>`` block (spec §14 Q2 resolution) so the page
is self-contained and deterministic — no external refs, no timestamps.
"""

from __future__ import annotations

from html import escape

from .geometry import BoardGeom, ClueEntry
from .svg import grid_svg

_STYLE = """\
:root { color-scheme: light dark; }
body {
  font-family: 'Helvetica Neue', Arial, sans-serif;
  margin: 2rem; line-height: 1.4;
}
.xword-title { font-size: 1.5rem; font-weight: 700; margin: 0 0 1rem; }
.xword-grid { display: block; height: auto; max-width: 100%; }
.xword-clues {
  display: flex; flex-wrap: wrap; gap: 2rem; margin-top: 1.5rem;
}
.xword-col { flex: 1 1 16rem; min-width: 0; }
.xword-col h2 {
  font-size: 1rem; text-transform: uppercase; letter-spacing: 0.05em;
  border-bottom: 1px solid currentColor; padding-bottom: 0.25rem; margin: 0 0 0.5rem;
}
.xword-col ol { margin: 0; padding-left: 2.5rem; }
.xword-col li { margin: 0.15rem 0; }
.xword-col li::marker { color: #0066cc; font-weight: 600; }
"""


def _clue_list(entries: list[ClueEntry], direction: str) -> str:
    items = "".join(
        f'<li value="{e.number}" id="word-{e.number}-{direction}">{escape(e.text)}</li>'
        for e in entries
    )
    return f"<ol>{items}</ol>"


def html_document(geom: BoardGeom) -> str:
    """A self-contained HTML page: title + inline SVG grid + clue lists.

    Blank/solved is baked into the geometry upstream (``board_geometry``).
    """
    title = geom.title or "Crossword"
    parts = [
        "<!doctype html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{escape(title)}</title>",
        f"<style>\n{_STYLE}</style>",
        "</head>",
        "<body>",
    ]
    if geom.title:
        parts.append(f'<h1 class="xword-title">{escape(geom.title)}</h1>')
    parts.append(grid_svg(geom, css_class="xword-grid"))
    if geom.across or geom.down:
        parts.append('<div class="xword-clues">')
        for heading, entries, direction in (
            ("Across", geom.across, "across"),
            ("Down", geom.down, "down"),
        ):
            if entries:
                parts.append(
                    f'<section class="xword-col"><h2>{heading}</h2>'
                    f"{_clue_list(entries, direction)}</section>"
                )
        parts.append("</div>")
    parts += ["</body>", "</html>", ""]
    return "\n".join(parts)
