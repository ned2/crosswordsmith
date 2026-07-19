"""Cross-format semantic Board comparison (spec §6.5)."""

from __future__ import annotations

import json
from typing import Any

from .board import ACROSS, Board, Cell, Word, enumeration_from_cells
from .formats import parse_board
from .util import dump_json

Difference = dict[str, Any]
WordKey = tuple[str, tuple[tuple[int, int], ...]]

_CELL_FIELDS = (
    ("letter", "letter"),
    ("number", "number"),
    ("across", "across_num"),
    ("down", "down_num"),
    ("circle", "circle"),
    ("barRight", "bar_right"),
    ("barBelow", "bar_below"),
    ("prefilled", "prefilled"),
)


def diff_boards(left: Board, right: Board) -> dict[str, Any]:
    """Return a stable, granular D12 comparison report."""
    differences: list[Difference] = []
    _add(differences, "dimensions.height", left.height, right.height)
    _add(differences, "dimensions.width", left.width, right.width)

    for r in range(min(left.height, right.height)):
        for c in range(min(left.width, right.width)):
            left_cell = left.grid[r][c]
            right_cell = right.grid[r][c]
            path = f"grid[{r}][{c}]"
            if (left_cell is None) != (right_cell is None):
                _add(
                    differences,
                    path,
                    "block" if left_cell is None else "light",
                    "block" if right_cell is None else "light",
                )
                continue
            if left_cell is None or right_cell is None:
                continue
            for field, attribute in _CELL_FIELDS:
                _add(
                    differences,
                    f"{path}.{field}",
                    getattr(left_cell, attribute),
                    getattr(right_cell, attribute),
                )
            _add(
                differences,
                f"{path}.style",
                _residual_style(left_cell),
                _residual_style(right_cell),
            )

    left_words = {_word_key(word): word for word in left.words}
    right_words = {_word_key(word): word for word in right.words}
    for key in sorted(set(left_words) | set(right_words), key=_word_sort_key):
        left_word = left_words.get(key)
        right_word = right_words.get(key)
        path = _word_path(key)
        if left_word is None or right_word is None:
            _add(
                differences,
                path,
                _word_value(left_word) if left_word is not None else None,
                _word_value(right_word) if right_word is not None else None,
            )
            continue
        _add(differences, f"{path}.number", left_word.number, right_word.number)
        _add(differences, f"{path}.clue", left_word.clue or "", right_word.clue or "")
        _add(
            differences,
            f"{path}.enumeration",
            left_word.enumeration or enumeration_from_cells(left_word.cells),
            right_word.enumeration or enumeration_from_cells(right_word.cells),
        )

    _add(differences, "meta.title", left.meta.get("title"), right.meta.get("title"))
    _add(differences, "meta.author", left.meta.get("author"), right.meta.get("author"))
    return {"differences": differences, "equal": not differences}


def _add(differences: list[Difference], path: str, left: Any, right: Any) -> None:
    if left != right:
        differences.append({"left": left, "path": path, "right": right})


def _residual_style(cell: Cell) -> dict[str, Any] | None:
    style = dict(cell.style or {})
    if cell.circle and style.get("shapebg") == "circle":
        style.pop("shapebg")
    barred = style.get("barred")
    if isinstance(barred, str):
        value = barred.upper()
        if value and all(ch in "LTRB" for ch in value):
            style.pop("barred")
    return style or None


def _word_key(word: Word) -> WordKey:
    return (word.direction, tuple(word.cells))


def _word_sort_key(key: WordKey) -> tuple[Any, ...]:
    direction, cells = key
    start = cells[0] if cells else (-1, -1)
    return (start[0], start[1], 0 if direction == ACROSS else 1, cells)


def _word_path(key: WordKey) -> str:
    direction, cells = key
    coords = "-".join(f"{r},{c}" for r, c in cells)
    return f"words[{direction}:{coords}]"


def _word_value(word: Word) -> dict[str, Any]:
    return {
        "clue": word.clue or "",
        "enumeration": word.enumeration or enumeration_from_cells(word.cells),
        "number": word.number,
    }


def format_diff_text(report: dict[str, Any]) -> str:
    """Render the canonical difference records as stable plain text."""
    differences = report["differences"]
    if not differences:
        return "equal: no semantic differences\n"
    n = len(differences)
    lines = [f"different: {n} semantic difference{'s' if n != 1 else ''}"]
    for difference in differences:
        left = json.dumps(difference["left"], ensure_ascii=False, sort_keys=True)
        right = json.dumps(difference["right"], ensure_ascii=False, sort_keys=True)
        lines.append(f"{difference['path']}: {left} -> {right}")
    return "\n".join(lines) + "\n"


def diff_texts(
    left_text: str,
    right_text: str,
    left_fmt: str | None = None,
    right_fmt: str | None = None,
    *,
    json_output: bool = False,
) -> tuple[str, bool]:
    """Parse and compare two layouts, returning (report text, differs)."""
    report = diff_boards(parse_board(left_text, left_fmt), parse_board(right_text, right_fmt))
    output = dump_json(report) if json_output else format_diff_text(report)
    return output, not report["equal"]
