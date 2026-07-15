#!/usr/bin/env python3
"""Score-quality comparison for a filled crossword grid.

Prototype fill-quality benchmark: crosswordsmith `fill` (scoreless MRV) vs
ingrid_core (scored CSP), both filling the SAME grid from the SAME word list
(Spread the Wordlist). Every placed entry is scored *post hoc* against STW, so
the quality delta is visible even though crosswordsmith emits no scores itself.

Metrics per fill: n entries, mean score, min score, #entries scored <50
("below clean"), #entries absent from STW (junk / score 0).

Usage:
    python3 score_fill.py stw_scored.txt open4 open5 mini7 mini9
Expects, per grid <g>: <g>.json (mask), cs_<g>.json (crosswordsmith full-dict
fill), cs50_<g>.json (crosswordsmith score>=50 prefiltered dict fill),
ingrid_<g>.txt (ingrid_core min-score 50 fill). Missing variants are skipped.
"""
import json
import sys
from pathlib import Path


def fold(word):
    """Engine-parity key: uppercase, keep A-Z only (digits/punct squeezed).

    crosswordsmith normalizes every dict word this way at load (README --dict),
    so post-hoc scoring must key the same way or entries like CWL's MP3J;100
    (played by the engine as MPJ) would look absent. STW ships pre-normalized,
    so this is a no-op for it (6 digit-bearing score-20 entries aside).
    """
    return ''.join(c for c in word.upper() if 'A' <= c <= 'Z')


def load_scores(path):
    scores = {}
    for line in Path(path).read_text().splitlines():
        if ';' in line:
            word, raw = line.rsplit(';', 1)
            try:
                score = int(raw)
            except ValueError:
                continue
            key = fold(word)
            if key:
                # duplicate folded keys keep the max, matching the engine
                scores[key] = max(score, scores.get(key, score))
    return scores


def slots_from_mask(mask):
    """Maximal horizontal+vertical white runs (len>=2) as lists of (r,c)."""
    rows, cols = len(mask), len(mask[0])
    blk = lambda r, c: mask[r][c] == '#'
    out = []
    for r in range(rows):
        c = 0
        while c < cols:
            if not blk(r, c):
                run = []
                while c < cols and not blk(r, c):
                    run.append((r, c))
                    c += 1
                if len(run) >= 2:
                    out.append(run)
            else:
                c += 1
    for c in range(cols):
        r = 0
        while r < rows:
            if not blk(r, c):
                run = []
                while r < rows and not blk(r, c):
                    run.append((r, c))
                    r += 1
                if len(run) >= 2:
                    out.append(run)
            else:
                r += 1
    return out


def ingrid_words(grid_txt, mask):
    lines = [l for l in Path(grid_txt).read_text().splitlines() if l]
    letters = {(r, c): ch.upper()
               for r, row in enumerate(lines) for c, ch in enumerate(row)}
    return [''.join(letters[cell] for cell in slot)
            for slot in slots_from_mask(mask)]


def cs_words(cs_json):
    doc = json.loads(Path(cs_json).read_text())
    return [fold(w['answer']) for w in doc['words']]


def stats(words, scores):
    vals = [scores.get(w, 0) for w in words]
    n = len(vals)
    return {
        'n': n,
        'mean': round(sum(vals) / n, 1) if n else 0,
        'min': min(vals) if vals else 0,
        'below50': sum(1 for v in vals if v < 50),
        'junk': sum(1 for w in words if w not in scores),
    }


def check_report(grid, s):
    """Cross-check the engine's --report-json against these post-hoc stats.

    The §8.4a acceptance gate: n/min/belowThreshold must match exactly;
    mean within 0.1 (both sides round to 1 d.p., but SWI rounds half away
    from zero while Python rounds half to even, so an exact .x5 tie makes
    the two 1-d.p. values differ by exactly 0.1 — e.g. 81.25 -> 81.3 vs
    81.2; the next legitimate gap is 0.2, so the 0.15 cutoff rejects any
    real disagreement). Returns True/False, or None if the engine wrote no
    report (older engine / failed fill).
    """
    rp = Path(f"cs_minscore_{grid}.report.json")
    if not rp.exists():
        return None
    rep = json.loads(rp.read_text())
    return (rep['n'] == s['n'] and rep['min'] == s['min']
            and rep['belowThreshold'] == s['below50']
            and abs(rep['mean'] - s['mean']) < 0.15)


def main():
    scores = load_scores(sys.argv[1])
    grids = sys.argv[2:]
    hdr = f"{'grid':7} {'tool':31} {'n':>3} {'mean':>5} {'min':>4} {'<50':>4} {'junk':>4}"
    print(hdr)
    print('-' * len(hdr))
    disagreements = 0
    for g in grids:
        mask = json.loads(Path(f"{g}.json").read_text())['mask']
        variants = [
            ('crosswordsmith (full)', f"cs_{g}.json", 'cs'),
            ('crosswordsmith (>=50 dict)', f"cs50_{g}.json", 'cs'),
            ('crosswordsmith (--min-score 50)', f"cs_minscore_{g}.json", 'cs'),
            ('ingrid_core (min-score 50)', f"ingrid_{g}.txt", 'ingrid'),
        ]
        for label, fpath, kind in variants:
            if not Path(fpath).exists():
                print(f"{g:7} {label:31} (missing)")
                continue
            words = cs_words(fpath) if kind == 'cs' else ingrid_words(fpath, mask)
            s = stats(words, scores)
            note = ''
            if fpath.startswith('cs_minscore_'):
                agree = check_report(g, s)
                if agree is True:
                    note = '  report-json AGREES'
                elif agree is False:
                    note = '  report-json DISAGREES'
                    disagreements += 1
            print(f"{g:7} {label:31} {s['n']:>3} {s['mean']:>5} "
                  f"{s['min']:>4} {s['below50']:>4} {s['junk']:>4}{note}")
        print()
    if disagreements:
        print(f"FAIL: {disagreements} --report-json disagreement(s) with post-hoc scoring")
        sys.exit(1)


if __name__ == '__main__':
    main()
