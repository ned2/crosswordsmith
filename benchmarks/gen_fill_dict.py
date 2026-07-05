#!/usr/bin/env python3
"""benchmarks/gen_fill_dict.py - derive frozen dictionary subsets from ENABLE.

The fill benchmark ladder needs dictionaries at several SIZES (a size-scaling
knob for the dict-load and search cost) that all preserve the crossing-
completion structure of a real wordlist. A subset must NOT be the alphabetical
first-N words (that skews the letter distribution - the front of a sorted list
is all a-words - which destroys crossing realism and can make a grid
artificially infeasible). Instead each subset is a SEEDED SHUFFLE of the full
ENABLE list, truncated to the target size, then re-SORTED so the committed file
is in fill's own canonical order (fill sorts+dedupes on load anyway, but a
sorted-on-disk file is diff-friendly and keeps the byte-freeze stable).

Determinism: the shuffle uses random.Random(SEED) with a FIXED seed recorded
below, so re-running reproduces byte-identical subsets. The committed subset
FILES are the source of truth for the bench (load-layer inference counts depend
on the exact bytes); this script documents how they were made, it is not run by
the bench.

Usage:
    python3 benchmarks/gen_fill_dict.py            # regenerate all subsets
    python3 benchmarks/gen_fill_dict.py --check    # verify committed == regenerated

Provenance: fixtures/dict/enable1.txt (ENABLE, public domain; see
fixtures/dict/README.md). All words lowercase ASCII a-z.
"""
import os
import random
import sys

SEED = 20260705                      # recorded shuffle seed (do not change)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SOURCE = os.path.join(ROOT, "fixtures", "dict", "enable1.txt")

# (target_size, output_basename). Target sizes are approximate labels; the
# actual count equals target_size exactly (truncation of the shuffled list).
SUBSETS = [
    (10000, "enable_10k.txt"),
    (25000, "enable_25k.txt"),
    (50000, "enable_50k.txt"),
]


def read_words(path: str) -> list[str]:
    with open(path, "r", encoding="ascii") as f:
        words = [line.strip() for line in f]
    return [w for w in words if w]


def derive(words: list[str], size: int, seed: int) -> list[str]:
    rnd = random.Random(seed)
    shuffled = words[:]
    rnd.shuffle(shuffled)
    subset = shuffled[:size]
    subset.sort()                    # canonical, diff-friendly on-disk order
    return subset


def write_words(path: str, words: list[str]) -> None:
    # ASCII, one word per line, trailing newline; LF only.
    with open(path, "w", encoding="ascii", newline="\n") as f:
        for w in words:
            _ = f.write(w + "\n")


def main(argv: list[str]) -> int:
    check = "--check" in argv
    words = read_words(SOURCE)
    if len(words) < max(s for s, _ in SUBSETS):
        sys.exit("source ENABLE has fewer words than the largest requested subset")
    rc = 0
    for size, base in SUBSETS:
        subset = derive(words, size, SEED)
        out = os.path.join(ROOT, "fixtures", "dict", base)
        if check:
            existing = read_words(out) if os.path.exists(out) else None
            status = "OK" if existing == subset else "MISMATCH"
            if status != "OK":
                rc = 1
            print(f"{base}: {status} ({len(subset)} words, seed {SEED})")
        else:
            write_words(out, subset)
            print(f"wrote {out}: {len(subset)} words (seed {SEED})")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
