#!/usr/bin/env python3
"""Generate a hard-but-satisfiable short-word crossword mesh fixture (idea I4).

Unlike the abandoned full-width lattice (I3), this builds a witness layout from
SHORT interlocking words by simulating the solver's own placement rules, so the
result is satisfiable AND reachable by construction:

  * seed a word at (0,0) across, so the puzzle is solvable at topleft_across;
  * repeatedly attach a short perpendicular word crossing an already-placed
    word, accepting it only if it passes the solver's checks (no flush-parallel
    neighbours on newly filled cells; the cells immediately before/after the
    word are empty or off-grid). Crossing cells reuse the existing letter.

Every word is therefore added by a legal solver move, i.e. the build order is a
valid placement order and a full layout exists.

Difficulty is tuned by:
  * a SMALL alphabet K -> many repeated letters -> baseline explores many
    spurious crossings (a branchy search) while MRV stays guided;
  * SHUFFLING the output word order, so input order is not the witness order
    and baseline cannot just replay the build.

Deterministic: a given (G, N, K, Lmin, Lmax, SEED) always yields the same words.

Usage:  gen_mesh_fixture.py G N K LMIN LMAX SEED OUTFILE
"""
import sys
import random


def build(G, N, K, lmin, lmax, seed):
    rng = random.Random(seed)
    alpha = [chr(ord('A') + i) for i in range(K)]
    grid = {}                      # (r, c) -> letter
    words = []                     # (cells, direction, string)
    seen = set()

    def ok(r, c):
        return 0 <= r < G and 0 <= c < G

    def attempt(cells, direction):
        """Return the word string if these cells form a legal placement, else None."""
        (r0, c0), (r1, c1) = cells[0], cells[-1]
        prev = (r0, c0 - 1) if direction == 'across' else (r0 - 1, c0)
        nxt = (r1, c1 + 1) if direction == 'across' else (r1 + 1, c1)
        if ok(*prev) and prev in grid:
            return None
        if ok(*nxt) and nxt in grid:
            return None
        letters, to_fill, crossings = [], [], 0
        for (r, c) in cells:
            if (r, c) in grid:                 # crossing: reuse, no adjacency check
                letters.append(grid[(r, c)])
                crossings += 1
            else:
                nbs = [(r - 1, c), (r + 1, c)] if direction == 'across' \
                    else [(r, c - 1), (r, c + 1)]
                if any(ok(nr, nc) and (nr, nc) in grid for (nr, nc) in nbs):
                    return None
                ll = rng.choice(alpha)
                letters.append(ll)
                to_fill.append((r, c, ll))
        return ''.join(letters), to_fill, crossings

    def commit(to_fill):
        for (r, c, ll) in to_fill:
            grid[(r, c)] = ll

    # Seed at (0,0) across.
    seed_len = rng.randint(lmin, lmax)
    seed_cells = [(0, c) for c in range(seed_len)]
    s, to_fill, _ = attempt(seed_cells, 'across')
    commit(to_fill)
    words.append((seed_cells, 'across', s))
    seen.add(s)

    attempts, max_attempts = 0, N * 4000
    while len(words) < N and attempts < max_attempts:
        attempts += 1
        pcells, pdir, _ = rng.choice(words)
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
        res = attempt(cells, ndir)
        if res is None:
            continue
        s, to_fill, crossings = res
        if crossings == 0 or s in seen:         # must actually cross; no dup answers
            continue
        # Reject prefix pairs: if one word is a prefix of another, the solver can
        # place them collinearly at the same start cell (a "word inside a word"),
        # which assign_clue_numbers/3 cannot number. (This is a latent SOLVER
        # bug the dense mesh surfaced; we sidestep it in the witness so the
        # fixture emits cleanly through the full CLI.)
        if any(w.startswith(s) or s.startswith(w) for w in seen):
            continue
        commit(to_fill)
        words.append((cells, ndir, s))
        seen.add(s)

    strings = [w[2] for w in words]
    rng.shuffle(strings)                        # break the witness order for baseline
    return strings, len(grid)


def main():
    G, N, K, lmin, lmax, seed, outfile = (
        int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]),
        int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7])
    words, filled = build(G, N, K, lmin, lmax, seed)
    with open(outfile, 'w') as f:
        f.write("%% Synthetic mesh fixture: %d words, grid %dx%d, alphabet %d, "
                "lengths %d-%d, seed %d.\n" % (len(words), G, G, K, lmin, lmax, seed))
        f.write("%% Generated by benchmarks/gen_mesh_fixture.py "
                "%d %d %d %d %d %d (idea I4).\n" % (G, N, K, lmin, lmax, seed))
        f.write("%% Short interlocking words read off a witness layout: "
                "satisfiable by construction,\n")
        f.write("%% hard for baseline (small alphabet + shuffled order), "
                "tractable for MRV.\n")
        f.write("clues([\n")
        f.write(",\n".join("       ['%s']" % w for w in words))  # no trailing comma
        f.write("\n      ]).\n")
    print("wrote %s: %d words, grid %d, %d cells filled, alphabet %d"
          % (outfile, len(words), G, filled, K))


if __name__ == '__main__':
    main()
