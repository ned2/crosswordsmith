# Implementation plan: scored fill (FS-1 + FS-3(a)/(c))

**Status: DONE (2026-07-15).** All phases (C1–C5) landed; the FS-3(a)
STW/ingrid harness gate PASSED (native `--min-score 50` matches the
`>=50`-dict column — mean 50.0 / min 50 / 0 below-clean — on all four
completable grids, `--report-json` cross-check agrees) and the fill ratchet
including `--heavy` rungs is PASS at +0.00% (search untouched). Built the
LOCKED §8.4a contract
([design-spec](../design-spec.md), DP-4 as amended by DP-5) per
[fill-scoring-uplift.md](./fill-scoring-uplift.md) §5 sequencing: DP-5 landed
first (own commit); this plan is FS-1 plus its FS-3(a)/(c) acceptance gates.
FS-4/5/6/7 are out of scope — **no search-power change of any kind** (budget,
restarts, forward-checking): §8.4a is one filter + one sort key, not a new
engine.

## Architecture (the whole feature is load + emit; the search is untouched)

Candidate order in `fill` is the order of words inside each `DictByLen`
length bucket (`candidates/4` maps index-ordset positions back into the
bucket; MRV counting is order-blind). So §8.4a lands entirely in
`load_dict` + the emit boundary:

1. **Ingestion (`fill.pl` `load_dict`).** Per-line `word;score` parse (split
   on the **last** `;`, mirroring `score_fill.py`): score = non-negative
   integer in native units (DP-5); the word part goes through the existing
   `normalize_word/2` pipeline unchanged. Unscored line → uniform score 1
   (DP-5). Malformed score (non-integer / negative) → whole line dropped +
   **counted, reported unconditionally on stderr** (INV-3, same channel as
   the Unicode drop report). Duplicate word → keep the **max** score
   (deterministic, order-independent; same silent dedupe class as today's
   `sort/2`).
2. **Prune.** `min_score(N)` (default 1, DP-5) filters `Score < N` after
   dedupe, before bucket/index build. Pruned count reported unconditionally
   on stderr when > 0 (INV-3) — so the DP-5 default stays quiet on unscored
   dicts (prunes nothing) and on scored dicts without a score-0 floor.
   If the prune empties the dictionary, the report adds the targeted hint
   (dict word count + observed max score — covers the `--min-score 2` vs
   uniform-score-1 degenerate case). Downstream, empty domains are the
   ordinary AC-FILL-1 infeasible outcome (slots named, exit non-zero).
3. **Ordering.** Each length bucket is sorted **score-descending, then the
   existing lexicographic word order** (msort over `s(NegScore, Letters)`
   keys — a total order, INV-2; ties at equal score ARE the mainline case
   and collapse to today's dictionary order). The index build is unchanged —
   index ordsets point into the reordered buckets, so search-time candidate
   enumeration follows automatically. **`fill_search`/MRV/counting kernels:
   zero edits.**
4. **Plain-dict fast path (AC-FILL-6).** If the file contains no `;` and
   `min_score =< 1`, take the existing load path verbatim (one C-level
   whole-string scan decides). Unscored dicts stay byte-identical end to end
   (goldens + `make fuzz`), and the gated bench seam `load_dict/3` stays
   inference-neutral on every existing rung (verify with
   `make bench-fill-check`; `search_inf` must be +0.00% everywhere).
5. **Scores out.** The loader additionally returns `uniform` or
   `scores(Assoc)` (normalized-word atom → score). `load_dict/3` keeps its
   exported signature (the ratchet seam) and delegates to the new
   options-carrying form (internal; tests reach it white-box — no new export
   for tests).
6. **Report (emit boundary, AC-FILL-7).** After a successful fill compute
   n / mean (1 d.p.) / min / below-threshold over the placed answers
   (normalized like ingestion; absent-from-dict — e.g. a non-dict seed —
   scores 0, mirroring `score_fill.py`). Threshold is **50** — the
   documented clean-floor convention (DP-5 recommends `--min-score 50`), so
   the numbers agree with `score_fill.py`'s `<50` column by construction.
   Channels: `--report-json FILE` writes the sorted-key JSON object
   `{"belowThreshold":K,"mean":X,"min":M,"n":N,"threshold":50}` whenever
   requested; the human summary rides the existing `--verbose` stderr
   summary (§5.1's quiet-on-success stderr contract governs — the report is
   informational, not a compromise, so it is not unconditional; the machine
   channel is `--report-json`). Stdout canonical layout: byte-untouched.
7. **CLI (`crosswordsmith` fill verb).** `--min-score N` (integer, default
   1) + `--report-json FILE`. Both documented in usage text + README
   (§8.4a is their spec home; §5's usage block is a sketch that already
   defers flag detail to the verb sections). **v1 artifact scope:** both new
   flags require the text `--dict` path — combining either with `--index`
   is a clean usage error (artifacts carry no scores; baking scored
   ordering/prune into artifacts is future scope, noted in the error).
   `--save-index` from a scored dict inherits ordering/prune for free (it
   calls `load_dict`) but is out of v1's documented surface for the same
   reason.

## Fixture (INV-4: authored from scratch, never derived from any scored list)

`fixtures/dict_scored_sample.txt` — an ORIGINAL union of the existing
3-letter sample vocabulary with hand-assigned scores, containing: a score-0
entry (blocklist floor, must never appear in a default fill), unscored lines
(uniform 1), one malformed line (reported, dropped), equal-score ties
(mainline), and a score spread that makes the score-descending ordering
observably change the chosen fill vs `wordlist_sample.txt` on
`fill_grid_3.json`, and makes `--min-score 50` change it again.

## Phases / commits (STATUS.md swept in each status-changing commit)

- **C1 `feat(fill)`** — ingestion + ordering + prune in `load_dict` (+ the
  fast path), scored fixture, plunit (parse/uniform/malformed/dedupe-max/
  ordering/prune/hint/score-0-default/byte-identity, AC-FILL-5/6/8).
- **C2 `feat(fill)`** — CLI `--min-score` + `--report-json` + report at the
  emit boundary; plunit for the report numbers; README fill section +
  limitations sweep ("fill is unscored" bullet).
- **C3 `test(fill)`** — FS-3(c) CI subset: goldens `fill_scored.json`
  (default flags), `fill_scored_min50.json` (`--min-score 50`), and the
  report JSON golden, wired into `run_tests.sh` + Makefile `golden`/
  `update-golden`; determinism-fuzz lines for the scored paths.
- **C4 `bench(fill)`** — FS-3(a): `cs_minscore_<g>.json` native variant in
  `benchmarks/fill_quality/run.sh` + `score_fill.py` (STW ingested natively
  via `--dict` + `--min-score 50`; report-json cross-checked against
  `score_fill.py` on the same fill). Needs external STW + ingrid_core: if
  absent here, wire it, run what runs, and record the gap in the harness
  README (best-effort STOP — do not fake numbers).
- **C5** — ratchet check (`make bench-fill-check`): `search_inf` +0.00%
  required; if `load_inf` moves past tolerance despite the fast path,
  re-record the load bucket explicitly with rationale (never silently).

Verification per commit: `make unit` while iterating; `make test` before
every commit; `make fuzz` after C1/C2. Goldens change only where C3 adds new
ones — no existing golden may move.
