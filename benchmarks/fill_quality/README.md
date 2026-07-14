# Fill-quality benchmark (prototype) — crosswordsmith `fill` vs ingrid_core

A prototype head-to-head measuring **fill quality**, not just legality:
crosswordsmith's scoreless MRV `fill` against
[`ingrid_core`](https://github.com/rf-/ingrid_core) (scored backtracking CSP),
both filling the **same grids** from the **same word list** — Spread the
Wordlist (STW), CC BY-NC-SA 4.0. Every placed entry is scored *post hoc* against
STW, so the quality gap is measurable even though crosswordsmith emits no scores.

This is the measured answer to the "load-bearing open question" in
[`docs/research/setter-tool-landscape-2026.md`](../../docs/research/setter-tool-landscape-2026.md)
§D: *does scoreless fill produce materially worse fill, or just occasionally
worse corners?* **Answer: materially worse — down to nonsense strings — plus a
separate completion-rate gap on hard grids.**

> **Not part of `make bench`.** It needs two external, non-bundled deps:
> `ingrid_core` (`cargo install ingrid_core`, needs a Rust toolchain) and the STW
> `word;score` list (download from <https://www.spreadthewordlist.com/>). Run:
> `STW=/path/to/spread-the-wordlist.txt ./run.sh`.

## Files

- `gen_grids.py` — emits the small benchmark grids in both crosswordsmith JSON
  and ingrid_core ASCII (`#`/`.`) form, so both tools fill identical masks.
- `score_fill.py` — scores a fill's entries against STW (n, mean, min, #below-50,
  #junk). Reconstructs ingrid's entries from its output grid + the mask.
- `run.sh` — driver: builds the plain + `score>=50` dicts, fills each grid three
  ways, prints the comparison.

## Results (STW snapshot 2026-07-01; score 50 = the "clean" benchmark)

### Fill quality, on grids both tools complete

| grid | crosswordsmith (full dict) | crosswordsmith (`score>=50` dict) | ingrid_core (`--min-score 50`) |
|---|---|---|---|
| open4 (4×4) | mean **41.2**, min 20, **5**/8 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below |
| open5 (5×5) | mean **27.0**, min 20, **10**/10 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below |
| mini7 (7×7) | mean **42.2**, min **0**, **9**/23 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below |
| mini9 (9×9) | mean **41.8**, min 20, **14**/28 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below |

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
| ingrid_core `--min-score 20` | completes in **10.7s** (mean 45.0, 17 below-50) |
| ingrid_core `--min-score 30` | completes in **9.7s** (mean 44.4, 22 below-50) |
| ingrid_core `--min-score 50` | **times out** (>90s) — even the scored CSP can't hit "clean" on this grid |

A **second, distinct** competitiveness axis: crosswordsmith's fixed-budget MRV
cannot complete a standard 13×13 with full-length slots that ingrid solves in
~10s. This is orthogonal to scoring — it is search power — and is the harder
engineering gap. (Note even ingrid can't fill `blocked_13a` at the clean-50 bar,
so that grid is genuinely hard, not just hard for crosswordsmith.)

## Caveats / how to extend

- **Small-grid sample.** Four easy grids + one hard grid; a real benchmark wants
  a spread of standard American masks (11×11, 15×15) and several STW snapshots.
- **`junk`=0 by construction** — crosswordsmith's dict here is *derived from* STW,
  so every placeable word has an STW score; the quality signal is mean/min/below-50.
- **Post-hoc scoring is the whole trick** — it works today, before crosswordsmith
  implements any scoring, and is exactly how to regression-test a future
  `fill --min-score` / score-tiebreak feature against this ingrid baseline.
- **Next step** once scored fill lands: add a `cs_minscore_<g>.json` variant
  (native `--min-score 50`) and confirm it matches the `>=50`-dict column here.
