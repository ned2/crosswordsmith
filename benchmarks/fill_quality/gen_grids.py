#!/usr/bin/env python3
"""Emit comparison grids in crosswordsmith JSON and ingrid_core ASCII forms."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

EASY_GRIDS = {
    "open4": ["....", "....", "....", "...."],
    "open5": [".....", ".....", ".....", ".....", "....."],
    "mini7": ["...#...", "...#...", ".......", "###.###",
              ".......", "...#...", "...#..."],
    "mini9": ["....#....", "....#....", ".........", "###...###", ".........",
              "###...###", ".........", "....#....", "....#...."],
}

# Grids that must satisfy the strict American-mask invariants below. The easy
# quality grids are exempt (mini9 has deliberate unchecked cells).
AMERICAN = {"amer11"}


def grid_documents():
    docs = {
        name: {"name": name, "size": len(mask),
               "symmetry": "rot180", "mask": mask}
        for name, mask in EASY_GRIDS.items()
    }
    docs["amer11"] = json.loads((ROOT / "grids/amer11.json").read_text())
    return docs


def check_american(name, mask):
    rows, cols = len(mask), len(mask[0])
    white = {(r, c) for r in range(rows) for c in range(cols)
             if mask[r][c] != '#'}
    assert all(len(row) == cols for row in mask), f"{name}: ragged mask"
    for r in range(rows):
        for c in range(cols):
            sym = mask[rows - 1 - r][cols - 1 - c]
            assert (mask[r][c] == '#') == (sym == '#'), \
                f"{name}: not 180-deg symmetric at {(r, c)}"
    runs = {}
    for horiz in (True, False):
        for a in range(rows if horiz else cols):
            run = []
            for b in range((cols if horiz else rows) + 1):
                cell = (a, b) if horiz else (b, a)
                if b < (cols if horiz else rows) and cell in white:
                    run.append(cell)
                else:
                    assert not run or len(run) >= 3, \
                        f"{name}: run shorter than 3 at {run}"
                    for cell2 in run:
                        runs.setdefault(cell2, 0)
                        runs[cell2] += 1
                    run = []
    assert all(runs.get(cell, 0) == 2 for cell in white), \
        f"{name}: unchecked cell(s) {[c for c in white if runs.get(c, 0) != 2]}"
    seen, todo = set(), [next(iter(white))]
    while todo:
        r, c = todo.pop()
        if (r, c) in seen or (r, c) not in white:
            continue
        seen.add((r, c))
        todo += [(r + 1, c), (r - 1, c), (r, c + 1), (r, c - 1)]
    assert seen == white, f"{name}: mask not connected"


def main(outdir="."):
    out = Path(outdir)
    docs = grid_documents()
    for name in AMERICAN:
        check_american(name, docs[name]["mask"])
    for name, doc in docs.items():
        mask = doc["mask"]
        (out / f"{name}.json").write_text(json.dumps(doc))
        (out / f"{name}.txt").write_text("\n".join(mask) + "\n")
        print(f"wrote {name}.json / {name}.txt ({len(mask)}x{len(mask[0])})")

if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
