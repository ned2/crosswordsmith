# Convenience targets for crosswordsmith.

SHELL := /bin/bash

.PHONY: test unit golden update-golden fuzz bench bench-matrix

BENCH_FORMAT ?= text
BENCH_ARGS ?=

# Full suite: plunit tests + golden-output regression.
test:
	./run_tests.sh

# Just the plunit unit/integration tests.
unit:
	swipl -q tests/run_tests.pl

# Determinism + degenerate-input fuzz (INV-2 / AC-EMIT-1 / AC-X-2 byte-identity,
# the §5.2 --out partial-write contract, and graceful handling of pathological
# inputs). On-demand: it spawns many CLI subprocesses, so it is NOT part of the
# fast `make test`.
fuzz:
	./tests/determinism_fuzz.sh

# Just the golden-output regression checks (the crosswordsmith CLI, end to end).
golden:
	@diff -u tests/golden/arrange_bundled_17_fixed.json \
		<(./crosswordsmith arrange --strict --size-mode fixed --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null) \
		&& echo "golden (arrange fixed): OK"
	@diff -u tests/golden/arrange_toc_demo_max.json \
		<(./crosswordsmith arrange --strict --size-mode max --size 25 --input fixtures/toc_demo.pl 2>/dev/null) \
		&& echo "golden (arrange max): OK"
	@diff -u tests/golden/arrange_bundled_17_fragment.json \
		<(./crosswordsmith arrange --strict --size-mode fixed --fragment fixtures/bundled_17_fragment.json --input fixtures/bundled_17_clues.pl 2>/dev/null) \
		&& echo "golden (arrange fragment): OK"
	@diff -u tests/golden/arrange_bundled_17_candidates.json \
		<(./crosswordsmith arrange --strict --size-mode fixed --candidates 3 --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null) \
		&& echo "golden (arrange candidates): OK"
	@diff -u tests/golden/lint_bundled_17_toc.json \
		<(./crosswordsmith lint --profile toc tests/golden/arrange_bundled_17_fixed.json 2>/dev/null) \
		&& echo "golden (lint toc): OK"
	@diff -u tests/golden/export_bundled_17.ipuz \
		<(./crosswordsmith export --to ipuz tests/golden/arrange_bundled_17_fixed.json 2>/dev/null) \
		&& echo "golden (export ipuz): OK"
	@diff -u tests/golden/export_bundled_17.exolve \
		<(./crosswordsmith export --to exolve tests/golden/arrange_bundled_17_fixed.json 2>/dev/null) \
		&& echo "golden (export exolve): OK"
	@diff -u tests/golden/fill_3.json \
		<(./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt 2>/dev/null) \
		&& echo "golden (fill 3x3): OK"

# Regenerate the golden files. Use only after an INTENTIONAL output change,
# and review the diff before committing.
update-golden:
	./crosswordsmith arrange --strict --size-mode fixed --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null > tests/golden/arrange_bundled_17_fixed.json
	./crosswordsmith arrange --strict --size-mode max --size 25 --input fixtures/toc_demo.pl 2>/dev/null > tests/golden/arrange_toc_demo_max.json
	./crosswordsmith arrange --strict --size-mode fixed --fragment fixtures/bundled_17_fragment.json --input fixtures/bundled_17_clues.pl 2>/dev/null > tests/golden/arrange_bundled_17_fragment.json
	./crosswordsmith arrange --strict --size-mode fixed --candidates 3 --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null > tests/golden/arrange_bundled_17_candidates.json
	./crosswordsmith lint --profile toc tests/golden/arrange_bundled_17_fixed.json 2>/dev/null > tests/golden/lint_bundled_17_toc.json
	./crosswordsmith export --to ipuz tests/golden/arrange_bundled_17_fixed.json 2>/dev/null > tests/golden/export_bundled_17.ipuz
	./crosswordsmith export --to exolve tests/golden/arrange_bundled_17_fixed.json 2>/dev/null > tests/golden/export_bundled_17.exolve
	./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt 2>/dev/null > tests/golden/fill_3.json
	@echo "Regenerated golden files (arrange fixed/max/fragment/candidates + lint toc + export ipuz/exolve + fill)"

# Product benchmark for `arrange`: end-to-end command latency, the in-process
# search alone, and the CLI-wrapper overhead between them (rest = command -
# search), over the workloads in benchmarks/workloads.pl. Core workloads only by
# default; results are machine-specific and reporting-only. Compare on the search
# inference counts (machine-independent); wall/rss are reporting-only.
#   make bench
#   make bench BENCH_ARGS=--heavy                       # + budget-saturating probes (~26s each)
#   make bench BENCH_FORMAT=csv BENCH_ARGS="--fixture bundled"
bench:
	swipl -q benchmarks/run_arrange.pl -- --format $(BENCH_FORMAT) $(BENCH_ARGS)

# Strategy x fixture comparison matrix (CSV on stdout). Each fixture runs on
# its manifest grid (benchmarks/fixtures.pl). Optionally restrict strategies:
#   make bench-matrix BENCH_STRATEGIES="baseline mrv_capped"
# Compare on inferences (machine-independent); wall time is reporting-only.
BENCH_STRATEGIES ?=
bench-matrix:
	swipl -q benchmarks/run_matrix.pl -- $(BENCH_STRATEGIES)
