# Fill-Quality Checks

This directory owns two independent checks for `crosswordsmith fill`:

- the permanent AC-FILL-12 quality gate over the hard `blocked_13a` reference
  rows;
- an optional crosswordsmith/ingrid_core comparison over five smaller masks.

The permanent gate requires an external Spread the Wordlist (STW) `word;score`
file. STW is CC BY-NC-SA 4.0 and is intentionally not bundled. The optional
comparison accepts a pre-normalized ASCII scored list. Both checks independently
rescore emitted answers with `score_fill.py` and fail if that result disagrees
with the engine's quality sidecar.

## Permanent Gate

Run the design-spec section 8.4c gate with:

```sh
make bench-fill-quality-check STW=/path/to/spread-the-wordlist.txt
```

The command first runs the no-data Python tests, then fills
`grids/blocked_13a.json` at `--min-score 30` and `--min-score 1`. It supplies no
budget, seed, or shuffle override, so it exercises the shipped default budget
and deterministic default search path.

| row | enforced floor |
|---|---|
| `blocked_13a @30` | 54 entries, mean >= 45.0, min >= 30, no junk |
| `blocked_13a @1` | 54 entries, mean >= 38.7, min >= 1, no junk |

The checker fails on an incomplete or budget-exhausted fill, malformed output,
a quality regression, or sidecar/post-hoc disagreement. Fresh temporary output
prevents stale files from satisfying the gate.

The required STW snapshot has 315,905 lines and SHA-256:

```text
9aa4563eb82f9af45f3a9c4fbc7fe3c8c3e27da6b76cc48f5e09c499bd17e4ee
```

The checker prints both properties and validates the digest before running; the
digest pins the line count and all bytes because the floors are snapshot-specific.
To run only checks that need no external dictionary:

```sh
make bench-fill-quality-test
```

The external-data gate is intentionally not part of `make test`.

## Optional Comparison

With `ingrid_core` on `PATH`, compare both engines over `open4`, `open5`,
`mini7`, `mini9`, and `amer11`:

```sh
STW=/path/to/spread-the-wordlist.txt benchmarks/fill_quality/run.sh
```

The driver reports crosswordsmith with the full plain dictionary, a prefiltered
score-50 dictionary, native `--min-score 50`, and ingrid_core at 50. Missing
fills remain reportable comparison data, but malformed output or scoring
disagreement fails the run. This comparison is not the AC-FILL-12 gate and may
be used with other pre-normalized ASCII scored wordlists.

## Recorded Results

The living docs cite these durable conclusions; detailed campaign mechanics are
kept in the linked evidence rather than in this operating guide:

- The digit-cleaned full CWL (566,665 words) exceeds SWI's default 1 GB stack
  while building the index. Its 437,400-word score-30 band loads but can exceed
  the stack during MAC propagation on full-length slots. The bundled 252,200-word
  score-50 derivative is inside both measured envelopes.
- On the five comparison masks, CWL score-50 fills recorded means 78.1-83.1
  versus ingrid_core's 77.0-78.8; the STW clean-floor paths produce mean/min 50
  with no below-clean entries. The original scoreless probe demonstrated the
  `AAAAA`/`AAAAQ` quality failure that score-aware fill closes.
- AC-FILL-12 closes `blocked_13a` at score floors 30 and 1. `blocked_13b` and
  `blocked_15a` remain out of reach for both engines and are not promised rows.

## Evidence

- Current quality contract and decisions: [`../../docs/design-spec.md`](../../docs/design-spec.md)
- Shipped MAC/dom-wdeg implementation: [`../../docs/plans/fill-mac-dwd-implementation.md`](../../docs/plans/fill-mac-dwd-implementation.md)
- Fill performance closeout: [`../../docs/research/fill-perf-program-closeout-2026-07.md`](../../docs/research/fill-perf-program-closeout-2026-07.md)
- B0 instrumentation result: [`../results/2026-07-16-fill-b0-mac-instrumentation.md`](../results/2026-07-16-fill-b0-mac-instrumentation.md)
- Retired-rig reconstruction recipes: [`../../docs/research/benchmark-probe-historical-reconstruction.md`](../../docs/research/benchmark-probe-historical-reconstruction.md)
