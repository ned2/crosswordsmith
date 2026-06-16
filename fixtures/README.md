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

The dense fixture is intentionally much slower than the others. On one
development machine it took about 23 seconds and 450 million inferences for a
single measured iteration; treat that as an order-of-magnitude guide rather
than a portable baseline.

Example runs:

```sh
make bench BENCH_FIXTURE=fixtures/benchmark_14_words.pl BENCH_GRID=17
make bench BENCH_FIXTURE=fixtures/benchmark_16_dense_words.pl BENCH_GRID=17 BENCH_ITERATIONS=1 BENCH_WARMUP=0
make bench BENCH_FIXTURE=fixtures/benchmark_26_words.pl BENCH_GRID=49
swipl -q benchmarks/run_benchmarks.pl -- --grid 49 fixtures/benchmark_26_words.pl
```

## Other Fixtures

`bundled_17_clues.pl` is the small human-readable clue set used by README
examples, golden-output regression, and default benchmark runs. It solves
almost instantly and is mainly useful as a compatibility check.
