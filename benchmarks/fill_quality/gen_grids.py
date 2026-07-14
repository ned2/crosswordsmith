#!/usr/bin/env python3
"""Emit the small benchmark grids in BOTH crosswordsmith (JSON) and ingrid_core
(ASCII #/.) form, so the two fillers run on identical masks.

These are deliberately EASY, short-slot grids: crosswordsmith's fixed-budget MRV
search cannot complete a standard 13x13 with full-length slots (see README), so
the fill-QUALITY comparison uses grids both tools finish, and the hard grid is
reported separately as a completion-rate/search-power data point.
"""
import json
from pathlib import Path

GRIDS = {
    "open4": ["....", "....", "....", "...."],
    "open5": [".....", ".....", ".....", ".....", "....."],
    "mini7": ["...#...", "...#...", ".......", "###.###",
              ".......", "...#...", "...#..."],
    "mini9": ["....#....", "....#....", ".........", "###...###", ".........",
              "###...###", ".........", "....#....", "....#...."],
}

def main(outdir="."):
    out = Path(outdir)
    for name, mask in GRIDS.items():
        (out / f"{name}.json").write_text(json.dumps(
            {"name": name, "size": len(mask), "symmetry": "rot180", "mask": mask}))
        (out / f"{name}.txt").write_text("\n".join(mask) + "\n")
        print(f"wrote {name}.json / {name}.txt ({len(mask)}x{len(mask[0])})")

if __name__ == "__main__":
    import sys
    main(*sys.argv[1:])
