# Fixtures

This directory contains Prolog `clues/1` fixtures used by examples,
regressions, and local benchmarking. Each fixture can be passed to the CLI with
`--input`, or to the benchmark harness with `make bench BENCH_FIXTURE=...`.

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
| `benchmark_70_mesh_words.pl` | 70 | 21 | Dense mesh | Baseline: timeout; MRV: low | 1 | Hard search where baseline times out but MRV solves (idea I4) |

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

`benchmark_70_mesh_words.pl` is a hard short-word mesh: `baseline` does not
finish (the `bench-matrix` harness caps each cell at 60 s and records
`timeout`), while `mrv_inc` solves it quickly — the suite's case for search
pruning on a larger grid. It is reproducible via
`benchmarks/gen_mesh_fixture.py 21 70 8 3 5 1 <out>` (deterministic; see
`docs/experiments.md`, idea I4).
