# Arrange probe infrastructure

This directory is campaign-only measurement infrastructure. It never modifies
product search code and its rows must not be added to `benchmarks/workloads.pl`.

## Rigs

`authority` calls unchanged product predicates. It is the only rig valid for
500M inference cutoff outcomes, placement probability, censoring, success
inferences, reward, or layout. An operation runs strict's two representatives
(`topleft_across`, `topright`) under one shared product budget. A standalone
corner resets input memos once and installs a seed once; operation runs preserve
the product's single mutable seed stream across both corners. Seed cleanup uses
`setup_call_cleanup/3`, including exception and time-limit exits.

`instrumented` is a benchmark-local exact replay. It preserves seed-word order,
equal-count bucket order, crossing proof order/multiplicity, capped stale
overcounts, product grid/boundary trails, and one-reset memo lifecycle. `lean`
counts nodes, decisions, places, unplaces, wipeouts, and max depth. `full` also
counts support-bucket transitions, crossing-proof versus legal-child duplicates,
exact canonical state revisits, states known dead only after normal exhaustive
failure, repeated failed-subtree nodes, parent-local post-failure capture, and
maximum count-map/letter-grid/boundary-grid sizes. Stable word IDs are input-order
integers. A canonical state is the grid size, sorted remaining IDs, and sorted
absolute `(ID,start,dir)` placements; it is never translation-normalized. Dual
fingerprints select copied-key buckets, but every possible hit is verified by
exact canonical equality. Per-run observation facts are cleaned after each twin.
It runs unbounded or under semantic `nodes`/`decisions` caps, never inside a
product inference cap. Its
`measured_inferences` is overhead evidence, not product success inferences.

`P-D0` is a second benchmark-local exact replay in `d0_support.pl`. It observes
every letter-sharing refresh after unchanged product `inc_counts/8`: previous
and exact buckets, capped proof versus geometry multiplicity, old residues,
newest-source proofs, candidate checks, and letter/boundary/adjacency watch
dirtiness. Its proposed classifier is observational only and is checked against
the full exact recount. Bucket 2 always falls back; no shadow value enters MRV
ordering or grid state.

## Commands

All commands run from the repository root and emit one JSON row on stdout plus
start/end heartbeats on stderr. `SEARCH_SEED` is an integer or `none`; the last
argument is an outer health timeout in seconds.

```sh
swipl -q benchmarks/probe_arrange/run.pl -- \
  authority-operation fixtures/bundled_17_clues.pl 17 6 0 none 500000000 \
  op-1 0 control 120

swipl -q benchmarks/probe_arrange/run.pl -- \
  authority-corner fixtures/bundled_17_clues.pl 17 6 0 7 500000000 \
  topleft_across op-2 0 seeded 120

swipl -q benchmarks/probe_arrange/run.pl -- \
  instrumented fixtures/bundled_17_clues.pl 17 6 0 7 topright lean none 0 \
  op-3 0 seeded 120

python3 benchmarks/probe_arrange/measure_process.py \
  instrumented fixtures/bundled_17_clues.pl 17 6 0 none topleft_across lean \
  none 0 op-4 0 overhead 120

swipl -q benchmarks/probe_arrange/profile.pl -- \
  fixtures/ladder_09x09_08w.pl 9 8 none topleft_across

swipl -q benchmarks/probe_arrange/verify_controls.pl
swipl -q benchmarks/probe_arrange/measure_overhead.pl -- \
  fixtures/bundled_17_clues.pl 17 6 none topleft_across lean 21

swipl -q benchmarks/probe_arrange/test_d0_support.pl
swipl -q benchmarks/probe_arrange/measure_d0.pl > /tmp/p-d0.json
swipl -q benchmarks/probe_arrange/pc0_duplicate_work.pl
```

An outer timeout is batch-health protection only. It emits
`outcome=interrupted`, `termination=interrupted`, and `censored=true`; it is
never reported as infeasible and preserves the configured numeric `cutoff`.
Search budget cutoffs are inference-based.

## Frozen inputs

`fixtures/` contains 12 exact-count cliff subjects. `check_fixtures.py`
regenerates each to a temporary directory, checks requested and committed word
counts, and byte-compares without overwriting committed files.

`search_seeds.json` freezes the first 64 SplitMix64 outputs from
`CROSSWOR_ASCII_U64 = 0x43524f5353574f52`: indices `[0,16)` pilot, `[16,32)`
tuning, and `[32,64)` held out. `check_seeds.py` checks both an independent
Python implementation and product `splitmix64/3`.

## Schema and pooling

`schema.py` validates JSONL. Unavailable fields are JSON `null`, particularly
all mechanism counters on authority rows and `success_inferences` on
instrumented rows. Its grouping guard rejects any aggregation group containing
both rigs. Authority and instrumented rows may be joined for replay checks but
must never be pooled for outcome/cutoff analysis.

`limit_kind` identifies the threshold unit. `cutoff` is that configured integer
threshold, or JSON `null` when `limit_kind=none`; it is independent of observed
termination. `termination` is `ok|budget|exhausted|interrupted`. Thus a placed
authority row still carries `cutoff=500000000, termination=ok`, while an outer
interruption retains the same cutoff with `termination=interrupted`.

```sh
make probe-arrange-fixtures
make probe-arrange-seeds
make probe-arrange-schema-test
python3 benchmarks/probe_arrange/schema.py traces.jsonl
```

## P-R0 pilot

`run_pr0.py` is the base-locked, serial launcher for the registered fixed-instance
pilot. It writes one fsynced authority row at a time, resumes by operation ID,
and stops before exceeding 5,400 child CPU-seconds. Its only early stop is the
registered all-censored condition after the first eight seeds.

```sh
python3 benchmarks/probe_arrange/run_pr0.py
python3 benchmarks/probe_arrange/schema.py \
  benchmarks/results/2026-07-17-p-r0-pilot.jsonl
python3 benchmarks/probe_arrange/analyze_pr0.py \
  benchmarks/results/2026-07-17-p-r0-pilot.jsonl
```

The launcher deliberately runs no instrumented row. The analysis treats
standalone corners as screening distributions only; it does not reproduce or
claim operation-wide memo or mutable-stream accounting.
