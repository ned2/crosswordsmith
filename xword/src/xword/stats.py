"""Deterministic descriptive statistics over a normalised Board (spec §6.4)."""

from __future__ import annotations

from collections import Counter
from typing import Any

from . import XwordError
from .board import ACROSS, DOWN, Board, Coord, enumeration_from_cells
from .formats import parse_board
from .util import dump_json


def board_stats(board: Board) -> dict[str, Any]:
    """Return the D11 descriptive report for ``board``."""
    lights: set[Coord] = set()
    filled = 0
    for r, row in enumerate(board.grid):
        for c, cell in enumerate(row):
            if cell is None:
                continue
            styled_bar = cell.style.get("barred") if cell.style else None
            if cell.bar_right or cell.bar_below or styled_bar:
                raise XwordError(
                    "stats: barred grids are not supported until entry derivation is bar-aware"
                )
            lights.add((r, c))
            if cell.letter is not None:
                filled += 1

    memberships: dict[Coord, set[str]] = {coord: set() for coord in lights}
    for word in board.words:
        for coord in word.cells:
            if coord in memberships:
                memberships[coord].add(word.direction)

    checked = sum(len(directions) >= 2 for directions in memberships.values())
    unchecked = sum(len(directions) == 1 for directions in memberships.values())
    orphan = sum(not directions for directions in memberships.values())
    light_count = len(lights)
    unfilled = light_count - filled

    if filled == 0:
        fill_state = "blank"
    elif unfilled == 0:
        fill_state = "solved"
    else:
        fill_state = "partial"

    length_counts = Counter(len(word.cells) for word in board.words)
    enumeration_counts = Counter(
        word.enumeration or enumeration_from_cells(word.cells) for word in board.words
    )
    components = _component_count(lights)
    deficit = sum(
        (board.height - 1 - r, board.width - 1 - c) not in lights for r, c in lights
    )

    return {
        "blocks": board.height * board.width - light_count,
        "cells": {
            "checked": checked,
            "filled": filled,
            "light": light_count,
            "orphan": orphan,
            "unchecked": unchecked,
            "unfilled": unfilled,
        },
        "connectivity": {
            "components": components,
            "connected": components <= 1,
        },
        "dimensions": {"height": board.height, "width": board.width},
        "fillState": fill_state,
        "symmetry": {"deficit": deficit, "rot180": deficit == 0},
        "words": {
            "across": sum(word.direction == ACROSS for word in board.words),
            "down": sum(word.direction == DOWN for word in board.words),
            "enumerations": [
                {"count": count, "enumeration": enumeration}
                for enumeration, count in sorted(enumeration_counts.items())
            ],
            "lengths": [
                {"count": count, "length": length}
                for length, count in sorted(length_counts.items())
            ],
            "total": len(board.words),
        },
    }


def _component_count(lights: set[Coord]) -> int:
    unseen = set(lights)
    components = 0
    while unseen:
        components += 1
        pending = [unseen.pop()]
        while pending:
            r, c = pending.pop()
            for neighbour in ((r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)):
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    pending.append(neighbour)
    return components


def format_stats_text(report: dict[str, Any]) -> str:
    """Render a stats report as fixed-order, terminal-independent text."""
    dimensions = report["dimensions"]
    cells = report["cells"]
    words = report["words"]
    connectivity = report["connectivity"]
    symmetry = report["symmetry"]
    lengths = ", ".join(
        f"{item['length']}:{item['count']}" for item in words["lengths"]
    ) or "none"
    enumerations = ", ".join(
        f"{item['enumeration']}:{item['count']}" for item in words["enumerations"]
    ) or "none"
    connection = "connected" if connectivity["connected"] else "disconnected"
    component_word = "component" if connectivity["components"] == 1 else "components"
    symmetry_name = "rot180" if symmetry["rot180"] else "asymmetric"
    return "\n".join(
        [
            f"Dimensions: {dimensions['height']}x{dimensions['width']}",
            f"Cells: {cells['light']} light, {report['blocks']} blocks; "
            f"{cells['filled']} filled, {cells['unfilled']} unfilled",
            f"Checking: {cells['checked']} checked, {cells['unchecked']} unchecked, "
            f"{cells['orphan']} orphan",
            f"Words: {words['total']} total, {words['across']} across, {words['down']} down",
            f"Lengths: {lengths}",
            f"Enumerations: {enumerations}",
            f"Fill: {report['fillState']}",
            f"Connectivity: {connection} ({connectivity['components']} {component_word})",
            f"Symmetry: {symmetry_name} (deficit {symmetry['deficit']})",
        ]
    ) + "\n"


def stats_text(text: str, from_fmt: str | None = None, *, json_output: bool = False) -> str:
    """Parse one layout and return its human or canonical-JSON stats report."""
    report = board_stats(parse_board(text, from_fmt))
    return dump_json(report) if json_output else format_stats_text(report)
