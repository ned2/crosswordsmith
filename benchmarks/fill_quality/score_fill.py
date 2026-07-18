#!/usr/bin/env python3
"""Independent score-quality comparison for filled crossword grids.

The optional comparison driver runs crosswordsmith and ingrid_core on the same
grid and Spread the Wordlist dictionary. Every placed entry is scored post hoc,
independently of crosswordsmith's quality sidecar.

Metrics per fill: n entries, mean score, min score, #entries scored <50
("below clean"), #entries absent from STW (junk / score 0).

Usage:
    python3 score_fill.py stw_scored.txt open4 open5 mini7 mini9 amer11
Expects, per grid <g>: <g>.json (mask), cs_<g>.json (crosswordsmith full-dict
fill), cs50_<g>.json (crosswordsmith score>=50 prefiltered dict fill),
ingrid_<g>.txt (ingrid_core min-score 50 fill). Missing variants are skipped.
"""
import json
import sys
from pathlib import Path


def fold(word):
    """Benchmark key: uppercase, keep A-Z only (digits/punct squeezed).

    The reference STW snapshot and measured CWL inputs are ASCII. Squeezing is
    still needed for entries like CWL's MP3J;100, which the engine plays as MPJ.
    The engine's product loader has broader Unicode folding; this scorer is
    deliberately scoped to the benchmark's pre-normalized ASCII dictionaries.
    """
    return ''.join(c for c in word.upper() if 'A' <= c <= 'Z')


def load_scores(path):
    scores = {}
    for line_number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        if ';' in line:
            word, raw = line.rsplit(';', 1)
            try:
                score = int(raw)
            except ValueError:
                continue
        else:
            word, score = line, 1
        if not word.isascii():
            raise ValueError(
                f"{path}:{line_number}: post-hoc scorer supports ASCII words only"
            )
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


def report_agrees(report_path, s):
    """Return whether one engine sidecar agrees with post-hoc stats."""
    rep = json.loads(Path(report_path).read_text())
    return (rep.get('threshold') == 50
            and rep.get('n') == s['n']
            and rep.get('min') == s['min']
            and rep.get('belowThreshold') == s['below50']
            and isinstance(rep.get('mean'), (int, float))
            and abs(rep['mean'] - s['mean']) < 0.15)


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
    return report_agrees(rp, s)


def main():
    try:
        scores = load_scores(sys.argv[1])
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
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
                else:
                    note = '  report-json MISSING'
                    disagreements += 1
            print(f"{g:7} {label:31} {s['n']:>3} {s['mean']:>5} "
                  f"{s['min']:>4} {s['below50']:>4} {s['junk']:>4}{note}")
        print()
    if disagreements:
        print(f"FAIL: {disagreements} --report-json failure(s) against post-hoc scoring")
        sys.exit(1)


if __name__ == '__main__':
    main()
