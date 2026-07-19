#!/usr/bin/env python3
"""Layout-quality analyzer (idea: cryptic-layout-spec.md, phase v1a).

Reads canonical layout JSON (the output of `crosswordsmith arrange`) on stdin and
prints layout-quality metrics as a JSON object on stdout. Measures how far a
layout sits from the cryptic-lattice ideal (~0.5 checked) and how renderable it
is (bounding box / aspect). No solver involved — pure post-processing.

A cell is CHECKED iff it carries both an across and a down clue number; UNCH
(for a word) iff it is in that word only.

Usage:  ./crosswordsmith arrange --size N --input FILE | python3 benchmarks/analyze_layout.py
"""
import json
import sys


def cell_at(grid, r, c):
    if 0 <= r < len(grid) and 0 <= c < len(grid[r]):
        return grid[r][c]
    return None


def is_checked(cell):
    return cell is not None and cell.get("across") is not None \
        and cell.get("down") is not None


def main():
    doc = json.load(sys.stdin)
    grid = doc["grid"]
    words = doc["words"]
    n_words = len(words)

    # --- cell-level checking ---
    filled = checked = 0
    rows = [r for r, row in enumerate(grid) for c in row if c is not None]
    cols = [ci for row in grid for ci, c in enumerate(row) if c is not None]
    for row in grid:
        for cell in row:
            if cell is not None:
                filled += 1
                if is_checked(cell):
                    checked += 1
    checked_fraction = checked / filled if filled else 0.0

    # --- per-word checking + unch runs ---
    ratios = []
    longest_unch_run = 0
    words_with_3run = 0
    words_meeting_half = 0
    unchecked_ends = 0
    for w in words:
        cells = w["cells"]
        L = len(cells)
        flags = [is_checked(cell_at(grid, r, c)) for (r, c) in cells]
        cc = sum(flags)
        ratios.append(cc / L if L else 0.0)
        if cc >= -(-L // 2):                       # ceil(L/2)
            words_meeting_half += 1
        if flags and (not flags[0] or not flags[-1]):
            unchecked_ends += 1
        run = best = 0
        for f in flags:
            run = 0 if f else run + 1
            best = max(best, run)
        longest_unch_run = max(longest_unch_run, best)
        if best >= 3:
            words_with_3run += 1

    # --- compactness ---
    if rows:
        h = max(rows) - min(rows) + 1
        w_ = max(cols) - min(cols) + 1
        bbox_area = h * w_
        fill_density = filled / bbox_area
        aspect = max(h, w_) / min(h, w_)
    else:
        h = w_ = bbox_area = 0
        fill_density = aspect = 0.0

    lengths = [len(w["cells"]) for w in words]
    out = {
        "words": n_words,
        "filled_cells": filled,
        "checked_cells": checked,
        "checked_fraction": round(checked_fraction, 3),
        "crossings_per_word": round(checked * 2 / n_words, 2) if n_words else 0,
        "word_checked_ratio_min": round(min(ratios), 3) if ratios else 0,
        "word_checked_ratio_mean": round(sum(ratios) / len(ratios), 3) if ratios else 0,
        "words_meeting_min_half": words_meeting_half,
        "frac_words_meeting_min_half": round(words_meeting_half / n_words, 3) if n_words else 0,
        "unchecked_end_words": unchecked_ends,
        "longest_unch_run": longest_unch_run,
        "words_with_3plus_unch_run": words_with_3run,
        "bbox": [h, w_],
        "bbox_area": bbox_area,
        "fill_density": round(fill_density, 3),
        "aspect": round(aspect, 2),
        "word_len_min": min(lengths) if lengths else 0,
        "word_len_mean": round(sum(lengths) / len(lengths), 1) if lengths else 0,
        "word_len_max": max(lengths) if lengths else 0,
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
