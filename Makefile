# Convenience targets for crosswordsmith.

SHELL := /bin/bash

.PHONY: test unit golden update-golden bench bench-matrix

BENCH_FIXTURE ?= fixtures/bundled_17_clues.pl
BENCH_GRID ?= 17
BENCH_START_LOC ?= topleft_across
BENCH_ITERATIONS ?= 30
BENCH_WARMUP ?= 3
BENCH_FORMAT ?= text
BENCH_STRATEGY ?=

# Full suite: plunit tests + golden-output regression.
test:
	./run_tests.sh

# Just the plunit unit/integration tests.
unit:
	swipl -q tests/run_tests.pl

# Just the golden-output regression checks (legacy crossword + arrange engine).
golden:
	@diff -u tests/golden/grid_17_topleft_across.txt \
		<(./crossword.pl --input fixtures/bundled_17_clues.pl 17 topleft_across 2>/dev/null) \
		&& echo "golden (crossword): OK"
	@diff -u tests/golden/arrange_bundled_17_fixed.json \
		<(swipl -q -g 'arrange_run("fixtures/bundled_17_clues.pl",17,fixed),halt' -t 'halt(1)' arrange.pl 2>/dev/null) \
		&& echo "golden (arrange fixed): OK"
	@diff -u tests/golden/arrange_toc_demo_max.json \
		<(swipl -q -g 'arrange_run("fixtures/toc_demo.pl",25,max),halt' -t 'halt(1)' arrange.pl 2>/dev/null) \
		&& echo "golden (arrange max): OK"

# Regenerate the golden files. Use only after an INTENTIONAL output change,
# and review the diff before committing.
update-golden:
	./crossword.pl --input fixtures/bundled_17_clues.pl 17 topleft_across 2>/dev/null > tests/golden/grid_17_topleft_across.txt
	swipl -q -g 'arrange_run("fixtures/bundled_17_clues.pl",17,fixed),halt' -t 'halt(1)' arrange.pl 2>/dev/null > tests/golden/arrange_bundled_17_fixed.json
	swipl -q -g 'arrange_run("fixtures/toc_demo.pl",25,max),halt' -t 'halt(1)' arrange.pl 2>/dev/null > tests/golden/arrange_toc_demo_max.json
	@echo "Regenerated golden files (crossword + arrange fixed/max)"

# Local performance baselines. Results are machine-specific and reporting-only.
# Benchmarks the production default strategy unless BENCH_STRATEGY is set, e.g.
#   make bench BENCH_FIXTURE=fixtures/benchmark_16_dense_words.pl BENCH_ITERATIONS=1 BENCH_WARMUP=0 BENCH_STRATEGY=baseline
bench:
	swipl -q benchmarks/run_benchmarks.pl -- \
		--grid $(BENCH_GRID) \
		--start-loc $(BENCH_START_LOC) \
		--iterations $(BENCH_ITERATIONS) \
		--warmup $(BENCH_WARMUP) \
		--format $(BENCH_FORMAT) \
		$(if $(BENCH_STRATEGY),--strategy $(BENCH_STRATEGY),) \
		$(BENCH_FIXTURE)

# Strategy x fixture comparison matrix (CSV on stdout). Each fixture runs on
# its manifest grid (benchmarks/fixtures.pl). Optionally restrict strategies:
#   make bench-matrix BENCH_STRATEGIES="baseline mrv_capped"
# Compare on inferences (machine-independent); wall time is reporting-only.
BENCH_STRATEGIES ?=
bench-matrix:
	swipl -q benchmarks/run_matrix.pl -- $(BENCH_STRATEGIES)
