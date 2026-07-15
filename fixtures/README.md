# Fixtures

This directory contains Prolog `clues/1` fixtures used by examples,
regressions, and local benchmarking. Each fixture can be passed to the CLI with
`--input`, or to the benchmark harness with `make bench BENCH_FIXTURE=...`.
It also holds the small `fill` test assets: grid masks (`fill_grid_*.json`),
seed fragments (`fill_seed_*.json`), and wordlists (`wordlist_sample.txt`,
`dict_cow*.txt`, and the scored `dict_scored_sample.txt` — an ORIGINAL,
hand-authored `word;score` list per INV-4, never derived from any published
scored wordlist; it deliberately contains a score-0 blocklist entry, unscored
lines, equal-score ties, and one malformed line for the §8.4a ingestion
tests). Larger unscored dictionaries live under [`dict/`](dict/README.md).

## Benchmark Fixtures

The benchmark fixtures are synthetic. Their answers are generated strings, not
human-readable crossword entries, so they are useful for solver workload
comparisons rather than puzzle content.

| fixture | words | grid | shape | expected cost | recommended iterations | use case |
| --- | ---: | ---: | --- | --- | ---: | --- |
| `benchmark_08_words.pl` | 8 | 13 | Sparse mesh | Very low | 30+ | Quick benchmark sanity check beyond the bundled example |
| `benchmark_14_words.pl` | 14 | 17 | Sparse mesh | Low | 30+ | Routine local comparisons with a deeper placement sequence |
| `benchmark_16_dense_words.pl` | 16 | 17 | Dense mesh | Very high | 1 | Stress testing search changes and pruning behavior |
| `benchmark_20_words.pl` | 20 | 37 | Comb | Low | 30+ | Larger word-count overhead without a branchy search tree |
| `benchmark_26_words.pl` | 26 | 49 | Comb | Low to moderate | 30+ | Larger grid and word-count scaling checks |
| `benchmark_70_mesh_words.pl` | 70 | 21 | Dense mesh | Baseline: low; MRV: moderate | 1 | Large dense mesh; baseline beats MRV here (idea I4; see I5 note) |

The "expected cost" column and the dense fixture's reputation describe the
**`baseline`** strategy (input-order search). Under baseline, the dense fixture
is intentionally much slower than the others: on one development machine it took
about 23 seconds and 450 million inferences for a single iteration — an
order-of-magnitude guide, not a portable baseline.

`make bench` now benchmarks the **production default strategy** (`mrv_inc`; see
`docs/experiments.md`), which solves every fixture here in well under 40 ms,
including the dense one. To reproduce the heavy baseline search, pass
`BENCH_STRATEGY=baseline`. The `recommended iterations` column likewise assumes
baseline (use 1 iteration for the dense fixture only under baseline).

Example runs:

```sh
make bench BENCH_FIXTURE=fixtures/benchmark_14_words.pl BENCH_GRID=17
# dense under the fast default strategy (mrv_inc):
make bench BENCH_FIXTURE=fixtures/benchmark_16_dense_words.pl BENCH_GRID=17
# dense under baseline (slow; one iteration, no warmup):
make bench BENCH_FIXTURE=fixtures/benchmark_16_dense_words.pl BENCH_GRID=17 BENCH_ITERATIONS=1 BENCH_WARMUP=0 BENCH_STRATEGY=baseline
make bench BENCH_FIXTURE=fixtures/benchmark_26_words.pl BENCH_GRID=49
# whole strategy x fixture matrix:
make bench-matrix
```

## Other Fixtures

`bundled_17_clues.pl` is the small human-readable clue set used by README
examples, golden-output regression, and default benchmark runs. It solves
almost instantly and is mainly useful as a compatibility check.

`benchmark_70_mesh_words.pl` is a large dense short-word mesh. It was added (I4)
as a baseline-hard fixture, but the I5 fix (collinear-overlap maximality)
revealed that its hardness was largely the bug: with the bug fixed, `baseline`
solves it in ~14 ms while the MRV strategies pay their per-node tax (so baseline
*beats* MRV here). It is kept as a large-grid mesh data point. Reproducible via
`benchmarks/gen_mesh_fixture.py 21 70 8 3 5 1 <out>` (deterministic; see
`docs/experiments.md`, ideas I4 and I5). `benchmark_16_dense_words.pl` remains
the suite's genuine baseline-hard fixture.

## Quality Fixtures

For *layout-quality* work (`docs/cryptic-layout-spec.md`) rather than
search-speed. The existing benchmark fixtures were built for search difficulty;
their interlock ceilings are incidental (combs are structurally low-ceiling).
These come with a **known-achievable checking ceiling** — the planted witness's
own quality, reported by `gen_mesh_fixture.py` — so a quality engine can be
measured against a real target, not just an unknown set. Analyze any layout with
`benchmarks/analyze_layout.py` (reads emitted JSON).

| fixture | words | grid | witness ceiling (checked / ≥half) | use |
| --- | ---: | ---: | --- | --- |
| `quality_22_mesh.pl` | 22 | 11 | 0.435 / 0.955 | small, solvable, high ceiling — fast iteration & baseline gap (current solver: 0.389 / 0.636) |
| `quality_61_mesh.pl` | 61 | 17 | 0.470 / 0.967 | dense; the backtracking solver can't pack it at all — a pure greedy-engine target |
| `toc_demo.pl` | 16 | 25 | low (real words; current ≈0.135) | realistic web-TOC input — graceful-degradation / honest-report end |

Reproducible: `gen_mesh_fixture.py 11 22 5 3 4 7 quality_22_mesh.pl` and
`... 17 70 6 3 4 3 quality_61_mesh.pl`. The cryptic ideal is ~0.5 checked with
~all words ≥half-checked; even maximally-dense witnesses top out ~0.47, so for an
arbitrary closed set the realistic target is "approach the witness ceiling",
not 0.5 (see the spec, §15).
