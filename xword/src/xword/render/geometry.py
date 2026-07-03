"""BoardGeom — renderer-agnostic grid geometry + clue-list layout (spec §5.2).

Computed once from a Board and shared by every renderer (terminal now; SVG,
HTML, and the Textual widget later). A renderer only decides how a
cell/number/block/clue is *drawn* in its medium.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..board import ACROSS, Board, enumeration_from_cells


@dataclass
class GeomCell:
    block: bool
    letter: str | None = None
    number: int | None = None
    circle: bool = False


@dataclass
class ClueEntry:
    number: int
    text: str  # "clue (enum)" — clue may be empty, never invented


@dataclass
class BoardGeom:
    rows: int
    cols: int
    cells: list[list[GeomCell]]
    title: str | None = None
    across: list[ClueEntry] = field(default_factory=list)
    down: list[ClueEntry] = field(default_factory=list)


def board_geometry(board: Board, *, blank: bool = False) -> BoardGeom:
    """Project a Board into drawable geometry; `blank` strips letters only."""
    cells = [
        [
            GeomCell(block=True)
            if cell is None
            else GeomCell(
                block=False,
                letter=None if blank else cell.letter,
                number=cell.number,
                circle=cell.circle,
            )
            for cell in row
        ]
        for row in board.grid
    ]
    across: list[ClueEntry] = []
    down: list[ClueEntry] = []
    for word in sorted(board.words, key=lambda w: w.number):
        enum = word.enumeration or enumeration_from_cells(word.cells)
        text = f"{word.clue} {enum}" if word.clue else enum
        entry = ClueEntry(number=word.number, text=text)
        (across if word.direction == ACROSS else down).append(entry)
    title = board.meta.get("title")
    return BoardGeom(
        rows=board.height,
        cols=board.width,
        cells=cells,
        title=title if isinstance(title, str) else None,
        across=across,
        down=down,
    )
