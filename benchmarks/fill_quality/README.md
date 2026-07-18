# Fill-quality gates

This directory owns two distinct checks for `crosswordsmith fill`:

- the permanent AC-FILL-12 quality gate over the hard `blocked_13a` reference
  rows; and
- an optional crosswordsmith/ingrid_core comparison over five smaller masks.

Both use an external Spread the Wordlist (STW) `word;score` file. STW is CC
BY-NC-SA 4.0 and is intentionally not bundled. The engine's quality sidecar is
always cross-checked against `score_fill.py`, which scores the emitted answers
independently.

## AC-FILL-12 gate

Run the standing design-spec section 8.4c gate with:

```sh
make bench-fill-quality-check STW=/path/to/spread-the-wordlist.txt
```

The command first runs the retained Python unit tests, then invokes the product
CLI on `grids/blocked_13a.json` at `--min-score 30` and `--min-score 1`. It
deliberately supplies no `--budget`, `--seed`, or `--shuffle`, so it exercises
the shipped default budget and deterministic default search path.

| row | enforced floor | reference observation |
|---|---|---|
| `blocked_13a @30` | 54 entries, mean >= 45.0, min >= 30, no junk | mean 45.0, min 30, 22 below 50 |
| `blocked_13a @1` | 54 entries, mean >= 38.7, min >= 1, no junk | mean 38.7, min 20, 33 below 50 |

The checker fails on a non-zero fill exit (including budget exhaustion), a
missing or malformed layout/report, a floor regression, or any disagreement
between the engine sidecar and independent post-hoc scoring. It uses a fresh
temporary directory, so stale output cannot make a failed run pass.

The recorded reference STW snapshot is 315,905 lines with SHA-256:

```text
9aa4563eb82f9af45f3a9c4fbc7fe3c8c3e27da6b76cc48f5e09c499bd17e4ee
```

The gate prints the supplied dictionary's line count and digest and requires
this exact snapshot. The floors are snapshot-specific; accepting a different
dictionary would confound wordlist drift with an engine regression. Use the
optional comparison driver for exploratory runs on other scored lists.

This external-data gate is not part of `make test`. Run its no-data unit layer
alone with:

```sh
make bench-fill-quality-test
```

## Optional comparison

`run.sh` compares crosswordsmith with
[`ingrid_core`](https://github.com/rf-/ingrid_core) on `open4`, `open5`,
`mini7`, `mini9`, and `amer11`:

```sh
STW=/path/to/spread-the-wordlist.txt benchmarks/fill_quality/run.sh
```

It requires `ingrid_core` on `PATH`. For each mask it reports crosswordsmith
with the full plain dictionary, a prefiltered score-50 dictionary, native
`--min-score 50`, and ingrid_core at 50. Missing fills remain comparison data;
sidecar/post-hoc disagreement exits non-zero. This driver is not the hard-row
AC-FILL-12 gate. The retained post-hoc scorer accepts pre-normalized ASCII
dictionaries (the measured STW and CWL inputs) and rejects broader Unicode
inputs rather than approximating the product loader's Unicode fold.

Files retained for the permanent checks:

- `check_ac_fill_12.py`: hard-row completion and quality gate.
- `score_fill.py`: independent post-hoc scorer and sidecar cross-check.
- `gen_grids.py`: emits JSON and ingrid ASCII forms for the comparison masks.
  `grids/amer11.json` is the single source of truth for that shipped mask.
- `run.sh`: optional five-mask ingrid comparison.
- `test_fill_quality.py`: no-data unit checks for the gate, scorer, and grid
  source relationship.

## Current measured facts

### STW clean-floor masks

The five comparison masks complete at native `--min-score 50` with zero entries
below 50. The 2026-07-16 replay recorded `n=8/10/23/28/44` for
`open4/open5/mini7/mini9/amer11`, with mean and minimum 50 on every mask.

### CWL on the section 8.4c engine (2026-07-16) - DP-9 grounding

The pinned, digit-cleaned Collaborative Word List measurement at score 50
recorded crosswordsmith means of 83.1, 78.5, 81.7, 80.9, and 78.1 on those
same five masks, all with minimum 50 and zero below 50. The `amer11` + bundled
`dicts/cwl50.dict` pairing is the supported demo; no loadable CWL floor filled
the UK stock grids in the recorded campaign.

The measured CWL capacity envelope remains:

- the cleaned 566,665-word full list exceeds SWI's default 1 GB stack while
  building the index;
- the 437,400-word score-30 band loads, but can exceed the stack during MAC
  propagation on grids with full-length slots; and
- the bundled 252,200-word score-50 derivative is inside both measured
  envelopes.

Since DP-10/AC-FILL-15, both overflow modes produce the ordinary one-line
capacity failure with remedies rather than a raw SWI stack dump.

## Historical evidence

The completed campaign rigs are intentionally not permanent product tooling.
Their results and exact source recovery points remain available:

| campaign | durable evidence |
|---|---|
| FS-1 scored-fill quality and FS-3(b)/FS-4 frontier | design-spec DP-4/DP-6 and the historical README/source recipes in `docs/research/benchmark-probe-historical-reconstruction.md` |
| DP-7 MAC/dom-wdeg adoption | `docs/plans/fill-mac-dwd-implementation.md` and historical probe commit `4653996` |
| DP-8 shipped hard-row result | this gate's floors plus `docs/plans/fill-mac-dwd-implementation.md` |
| B0 shipped-MAC instrumentation | `benchmarks/results/2026-07-16-fill-b0-mac-instrumentation.md` |
| B1/B2/B3/C1 rejected arms | `docs/experiments.md` and `docs/research/fill-perf-program-closeout-2026-07.md` |

See
[`../../docs/research/benchmark-probe-historical-reconstruction.md`](../../docs/research/benchmark-probe-historical-reconstruction.md)
for disposable-worktree reconstruction commands. Git history is the source
archive; this directory keeps only checks that still gate current behavior.
