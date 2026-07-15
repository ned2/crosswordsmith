# Fill-quality benchmark — crosswordsmith `fill` vs ingrid_core

A head-to-head measuring **fill quality**, not just legality: crosswordsmith's
`fill` (plain-dict scoreless mode AND the §8.4a native scored mode) against
[`ingrid_core`](https://github.com/rf-/ingrid_core) (scored backtracking CSP),
all filling the **same grids** from the **same word list** — Spread the
Wordlist (STW), CC BY-NC-SA 4.0. Every placed entry is scored *post hoc*
against STW, independently of anything the engine reports; the engine's own
`--report-json` is cross-checked against that post-hoc scoring. This harness
is the standing regression for `fill --min-score` (design-spec §8.4a).

This is the measured answer to the "load-bearing open question" in
[`docs/research/setter-tool-landscape-2026.md`](../../docs/research/setter-tool-landscape-2026.md)
§D: *does scoreless fill produce materially worse fill, or just occasionally
worse corners?* **Answer: materially worse — down to nonsense strings — plus a
separate completion-rate gap on hard grids.**

> **Not part of `make bench`.** It needs two external, non-bundled deps:
> `ingrid_core` (`cargo install ingrid_core`, needs a Rust toolchain) and the STW
> `word;score` list (download from <https://www.spreadthewordlist.com/>). Run:
> `STW=/path/to/spread-the-wordlist.txt ./run.sh`. Any `word;score` list works
> as `$STW` — the CWL measurement below reused the same driver.

## Files

- `gen_grids.py` — emits the small benchmark grids in both crosswordsmith JSON
  and ingrid_core ASCII (`#`/`.`) form, so both tools fill identical masks.
- `score_fill.py` — scores a fill's entries against STW (n, mean, min, #below-50,
  #junk). Reconstructs ingrid's entries from its output grid + the mask. Also
  cross-checks the native fill's `--report-json` sidecar against its own
  post-hoc stats (n/min/belowThreshold exact; mean within 0.1 — SWI rounds
  half away from zero, Python half to even, so a .x5 tie differs by exactly
  0.1 in the last decimal) and exits non-zero on disagreement. Scores are
  keyed by the engine's own normalization (A–Z fold, digits/punct squeezed,
  duplicate keys keep the max score) so post-hoc lookups match what the
  engine actually played.
- `run.sh` — driver: builds the plain + `score>=50` dicts, fills each grid four
  ways (crosswordsmith full / `>=50`-prefiltered dict / **native
  `--min-score 50`** / ingrid_core), prints the comparison.

## Results (STW snapshots 2026-07-01 / 2026-07-15; score 50 = the "clean" benchmark)

### Fill quality, on grids both tools complete

The **native `--min-score 50`** column (§8.4a scored fill, added 2026-07-15)
is the FS-1 acceptance gate: it must — and does — match the `>=50`-dict
column (mean/min 50, zero below-clean) on every grid, and its
`--report-json` sidecar agreed with `score_fill.py`'s post-hoc stats on all
four fills. Full-dict numbers were identical across both snapshots.

| grid | crosswordsmith (full dict) | crosswordsmith (`score>=50` dict) | crosswordsmith (`--min-score 50`, native) | ingrid_core (`--min-score 50`) |
|---|---|---|---|---|
| open4 (4×4) | mean **41.2**, min 20, **5**/8 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below ✅ | mean 50, min 50, 0 below |
| open5 (5×5) | mean **27.0**, min 20, **10**/10 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below ✅ | mean 50, min 50, 0 below |
| mini7 (7×7) | mean **42.2**, min **0**, **9**/23 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below ✅ | mean 50, min 50, 0 below |
| mini9 (9×9) | mean **41.8**, min 20, **14**/28 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below ✅ | mean 50, min 50, 0 below |

Scoreless MRV grabs the alphabetically-first fitting strings. The `open5` fill is
the vivid case:

```
crosswordsmith (full)   : AAAAA AAAAQ ALEDA ABAAB AMEAL ASALE ALBMS AEAEA ADAAL QABLE
crosswordsmith (>=50)   : AAHED AIART ANNIE HEDDA EAGER DROSS INEAR ANDGO RIDES TEARS
ingrid_core (min 50)    : ASSTS SHERA SAMUS ERIES SISSY
```

**Two findings:**

1. **Scoreless fill is materially worse** — `AAAAA`, `AAAAQ`, `ABAAB`, `AEAEA`
   are non-words. Mean scores sit at 27–42 with entries at score 0.
2. **A crude `score>=50` *dict prefilter* fully recovers ingrid-parity quality**
   (mean/min 50, zero below-clean) on every completable grid — the cheapest form
   of scored fill (filter `--dict` by score, or an in-engine `--min-score`)
   already closes the quality gap here.

### Completion / search power, on a hard grid (`blocked_13a`, full 13-length slots)

| tool | result |
|---|---|
| crosswordsmith (any dict) | **does not complete** — budget-exhausted at ~20s |
| crosswordsmith `--min-score 50` (native, 2026-07-15) | **does not complete** — budget-exhausted at ~21s (the prune removes 195,725 of 315,903 words and shrinks every domain; scoring never lifts the search ceiling, per §8.4a) |
| ingrid_core `--min-score 20` | completes in **10.7s** (mean 45.0, 17 below-50) |
| ingrid_core `--min-score 30` | completes in **9.7s** (mean 44.4, 22 below-50) |
| ingrid_core `--min-score 50` | **times out** (>90s) — even the scored CSP can't hit "clean" on this grid |

A **second, distinct** competitiveness axis: crosswordsmith's fixed-budget MRV
cannot complete a standard 13×13 with full-length slots that ingrid solves in
~10s. This is orthogonal to scoring — it is search power — and is the harder
engineering gap. (Note even ingrid can't fill `blocked_13a` at the clean-50 bar,
so that grid is genuinely hard, not just hard for crosswordsmith.)

## CWL measurement (2026-07-15) — evidence for the DP-6 bundling candidate

The FS-6(b) research confirmed the [Collaborative Word
List](https://github.com/Crossword-Nexus/collaborative-word-list) (CWL) is
**MIT + fully scored** (567,657 `word;score` entries, integer 0–100), making
it the first bundleable scored candidate (design-spec §8.5). The harness was
run with CWL as `$STW` to convert its two recorded caveats — staleness
(frozen 2023-02) and undocumented band semantics — into data:

| grid | crosswordsmith (`--min-score 50`, native) | ingrid_core (`--min-score 50`) |
|---|---|---|
| open4 | mean **81.2**, min 50, 0 below, 0 junk, report-json agrees | mean 78.8, min 50, 0 below |
| open5 | mean **78.5**, min 50, 0 below, 0 junk, report-json agrees | mean 77.0, min 50, 0 below |
| mini7 | mean **82.4**, min 50, 0 below, 0 junk, report-json agrees | mean 77.6, min 50, 0 below |
| mini9 | mean **78.0**, min 50, 0 below, 0 junk, report-json agrees | mean 77.1, min 50, 0 below |

**Findings, in decreasing order of importance for DP-6:**

1. **Capacity ceiling (engine finding, list-independent):** the full 567k-word
   CWL dict **crashes `build_index`** at SWI's default 1GB stack limit — a raw
   stack-overflow ERROR at dict load, not a clean §5.1-style failure report.
   STW's 315k loads fine, so the ceiling sits somewhere in between. Any list
   this size needs `--min-score` (CWL at 50 → 252,564 words, loads cleanly);
   the unpruned "full dict" column is therefore **unmeasurable for CWL**.
   Worth an engine robustness item: fail cleanly (or stream the index build)
   instead of overflowing.
2. **Quality: drop-in clean.** Native `--min-score 50` over raw CWL fills all
   four grids at min 50 / 0 below-clean, with `--report-json` agreeing with
   post-hoc scoring, and means (78–82) at or above ingrid_core on the same
   list. On these small grids the 2023 freeze does not hurt: CWL's clean band
   behaves like STW's.
3. **Digit-entry hygiene:** CWL keeps digit-bearing entries (`MP3J;100`,
   `0AD;25`, 992 such lines). The engine's documented A–Z fold squeezes
   digits, so `MP3J;100` becomes a 3-letter `MPJ` at score **100** — which
   score-descending ordering then plays *eagerly* (the open4 fill opens with
   `MPJS`). A CWL-derived dict should pre-filter digit-bearing entries (or
   the bundling pass should ship a cleaned derivative). This is also why
   `score_fill.py` keys scores by the engine's fold — raw-keyed lookups
   mislabel such entries as junk.
4. **Search power unchanged (expected):** `blocked_13a` with CWL at
   `--min-score 50` is still not proven within budget (~26s) — the ceiling is
   the search, not the list.

## Caveats / how to extend

- **Small-grid sample.** Four easy grids + one hard grid; a real benchmark wants
  a spread of standard American masks (11×11, 15×15) and several STW snapshots.
- **`junk`=0 by construction** — crosswordsmith's dict here is *derived from* STW,
  so every placeable word has an STW score; the quality signal is mean/min/below-50.
- **Post-hoc scoring is the independent check** — `score_fill.py` scores fills
  without trusting the engine, which is exactly how the native
  `cs_minscore_<g>.json` variant and its `--report-json` are gated (the two
  must agree; `score_fill.py` exits non-zero when they don't).
- **CI subset:** this head-to-head needs the two external deps and stays
  on-demand; the no-deps regression for the same §8.4a behavior (AC-FILL-5/-7)
  lives in `make test` — plunit + goldens over the bundled original scored
  fixture (`fixtures/dict_scored_sample.txt`, `tests/golden/fill_scored*`).
- **Mask spread (FS-3(b)) is still open** — standard 11×11/15×15 masks and a
  completion-rate × min-score matrix remain the recorded next growth step;
  they are a search-power measurement, not part of the scored-fill (FS-1)
  acceptance.
