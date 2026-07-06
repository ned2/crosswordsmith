#!/usr/bin/env python3
"""Generate a realistic real-word arrange product fixture (rework plan §6/§9).

Same planted-witness construction as gen_mesh_fixture.py — seed a word at
(0,0) across, then repeatedly attach a perpendicular crossing word by a legal
solver move, so the instance is satisfiable AND reachable by construction —
except candidate words are REAL English words drawn from the frozen
public-domain ENABLE list (fixtures/dict/enable1.txt, see its README) instead
of random small-alphabet strings. This gives the product bench what the
stress-test asked for: word sets a real user would feed `arrange` at the
blocked-standard daily sizes (13/15) — a long anchor plus 4–9 letter fills —
rather than tiny-alphabet synthetic meshes.

Construction differences from gen_mesh_fixture.py:
  * the seed word is a LONG anchor of length G-2 (a blocked daily's marquee
    entry), the attached words sample lengths uniformly in [LMIN, LMAX];
  * a candidate must match the grid at EVERY already-filled cell it covers
    (multi-crossings are accepted, not skipped), looked up via a
    (length, position, letter) index over the dictionary;
  * answers are unique and prefix-pair-free (same heuristic as the mesh
    generator, for the same reasons).

Deterministic: a given (dict bytes, G, N, LMIN, LMAX, SEED) always yields the
same fixture. The committed fixture FILE is the source of truth (baselines are
defined against its exact words); this script documents how it was made.

Usage:  gen_real_fixture.py G N DICT LMIN LMAX SEED OUTFILE
   e.g. gen_real_fixture.py 15 20 fixtures/dict/enable1.txt 4 9 11 fixtures/real_15x15_20w.pl
"""
import sys
import random
from collections import defaultdict


def load_dict(path, lmin, lmax, anchor_len):
    """Words (uppercased) of the fill lengths plus the anchor length, and a
    (length, position, letter) -> [word] index for crossing lookups."""
    wanted = set(range(lmin, lmax + 1)) | {anchor_len}
    words, seen = [], set()
    with open(path) as f:
        for line in f:
            w = line.strip()
            if not w or not w.isascii() or not w.isalpha():
                continue
            w = w.upper()
            if len(w) in wanted and w not in seen:
                seen.add(w)
                words.append(w)
    index = defaultdict(list)
    for w in words:
        for i, ch in enumerate(w):
            index[(len(w), i, ch)].append(w)
    return words, index


def build(G, N, dict_path, lmin, lmax, seed):
    rng = random.Random(seed)
    anchor_len = G - 2
    words_all, index = load_dict(dict_path, lmin, lmax, anchor_len)
    anchors = [w for w in words_all if len(w) == anchor_len]
    if not anchors:
        sys.exit("no dictionary word of anchor length %d" % anchor_len)

    grid = {}                      # (r, c) -> letter
    placed = []                    # (cells, direction, word)
    seen = set()

    def ok(r, c):
        return 0 <= r < G and 0 <= c < G

    def legal(cells, direction, word):
        """True iff placing `word` on `cells` is a legal solver move that
        crosses at least one existing word."""
        (r0, c0), (r1, c1) = cells[0], cells[-1]
        prev = (r0, c0 - 1) if direction == 'across' else (r0 - 1, c0)
        nxt = (r1, c1 + 1) if direction == 'across' else (r1 + 1, c1)
        if ok(*prev) and prev in grid:
            return False
        if ok(*nxt) and nxt in grid:
            return False
        crossings = 0
        for (r, c), ch in zip(cells, word):
            if (r, c) in grid:                 # crossing: letters must agree
                if grid[(r, c)] != ch:
                    return False
                crossings += 1
            else:
                nbs = [(r - 1, c), (r + 1, c)] if direction == 'across' \
                    else [(r, c - 1), (r, c + 1)]
                if any(ok(nr, nc) and (nr, nc) in grid for (nr, nc) in nbs):
                    return False
        return crossings >= 1

    def commit(cells, word):
        for (r, c), ch in zip(cells, word):
            grid[(r, c)] = ch

    def prefix_clash(w):
        return any(s.startswith(w) or w.startswith(s) for s in seen)

    # Long anchor at (0,0) across.
    anchor = rng.choice(sorted(anchors))
    cells = [(0, c) for c in range(anchor_len)]
    commit(cells, anchor)
    placed.append((cells, 'across', anchor))
    seen.add(anchor)

    attempts, max_attempts = 0, N * 4000
    while len(placed) < N and attempts < max_attempts:
        attempts += 1
        pcells, pdir, _ = rng.choice(placed)
        (cr, cc) = rng.choice(pcells)
        ndir = 'down' if pdir == 'across' else 'across'
        wl = rng.randint(lmin, lmax)
        p = rng.randrange(wl)                   # crossing offset within new word
        if ndir == 'across':
            cells = [(cr, cc - p + i) for i in range(wl)]
        else:
            cells = [(cr - p + i, cc) for i in range(wl)]
        if any(not ok(r, c) for (r, c) in cells):
            continue
        pool = index.get((wl, p, grid[(cr, cc)]), [])
        if not pool:
            continue
        for w in rng.sample(pool, min(40, len(pool))):
            if w in seen or prefix_clash(w):
                continue
            if legal(cells, ndir, w):
                commit(cells, w)
                placed.append((cells, ndir, w))
                seen.add(w)
                break

    strings = [p[2] for p in placed]
    rng.shuffle(strings)                        # break the witness order
    return strings, len(grid), placed, grid


def witness_metrics(placements, grid):
    """Same layout-quality summary as gen_mesh_fixture.py (the planted
    witness's known-achievable ceiling)."""
    cover = {}
    for (cells, direction, _s) in placements:
        for (r, c) in cells:
            cover.setdefault((r, c), set()).add(direction)
    filled = len(grid)
    checked = sum(1 for dirs in cover.values() if len(dirs) >= 2)
    meet = 0
    for (cells, _d, _s) in placements:
        cc = sum(1 for rc in cells if len(cover.get(rc, ())) >= 2)
        if cc >= -(-len(cells) // 2):           # ceil(L/2)
            meet += 1
    rows = [r for (r, _c) in grid]
    cols = [c for (_r, c) in grid]
    h = max(rows) - min(rows) + 1
    w = max(cols) - min(cols) + 1
    return {
        "checked_fraction": round(checked / filled, 3) if filled else 0.0,
        "frac_words_meeting_min_half": round(meet / len(placements), 3) if placements else 0.0,
        "fill_density": round(filled / (h * w), 3) if h and w else 0.0,
        "aspect": round(max(h, w) / min(h, w), 2) if h and w else 0.0,
    }


def main():
    G, N, dict_path, lmin, lmax, seed, outfile = (
        int(sys.argv[1]), int(sys.argv[2]), sys.argv[3],
        int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7])
    words, filled, placements, grid = build(G, N, dict_path, lmin, lmax, seed)
    with open(outfile, 'w') as f:
        f.write("%% Real-word product fixture: %d words, grid %dx%d, "
                "lengths %d-%d + anchor %d, seed %d.\n"
                % (len(words), G, G, lmin, lmax, G - 2, seed))
        f.write("%% Generated by benchmarks/gen_real_fixture.py "
                "%d %d %s %d %d %d (rework plan par.6/par.9).\n"
                % (G, N, dict_path, lmin, lmax, seed))
        f.write("%% ENABLE dictionary words planted on a legal witness layout: "
                "satisfiable and\n")
        f.write("%% reachable by construction; a realistic arrange input at a "
                "blocked daily size.\n")
        f.write("clues([\n")
        f.write(",\n".join("       ['%s']" % w for w in words))
        f.write("\n      ]).\n")
    wm = witness_metrics(placements, grid)
    print("wrote %s: %d words, grid %d, %d cells filled"
          % (outfile, len(words), G, filled))
    print("  witness ceiling: checked_fraction=%s frac_words>=half=%s "
          "fill_density=%s aspect=%s"
          % (wm["checked_fraction"], wm["frac_words_meeting_min_half"],
             wm["fill_density"], wm["aspect"]))


if __name__ == '__main__':
    main()
