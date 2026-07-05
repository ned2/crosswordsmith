# Convenience targets for crosswordsmith.

SHELL := /bin/bash

.PHONY: test unit golden update-golden fuzz bench bench-check bench-record bench-log bench-history bench-matrix

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
		<(./crosswordsmith arrange --strict --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null) \
		&& echo "golden (arrange fixed): OK"
	@diff -u tests/golden/arrange_toc_demo_max.json \
		<(./crosswordsmith arrange --strict --max-size 25 --input fixtures/toc_demo.pl 2>/dev/null) \
		&& echo "golden (arrange max): OK"
	@diff -u tests/golden/arrange_bundled_17_fragment.json \
		<(./crosswordsmith arrange --strict --fragment fixtures/bundled_17_fragment.json --input fixtures/bundled_17_clues.pl 2>/dev/null) \
		&& echo "golden (arrange fragment): OK"
	@diff -u tests/golden/arrange_bundled_17_candidates.json \
		<(./crosswordsmith arrange --strict --candidates 3 --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null) \
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
	./crosswordsmith arrange --strict --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null > tests/golden/arrange_bundled_17_fixed.json
	./crosswordsmith arrange --strict --max-size 25 --input fixtures/toc_demo.pl 2>/dev/null > tests/golden/arrange_toc_demo_max.json
	./crosswordsmith arrange --strict --fragment fixtures/bundled_17_fragment.json --input fixtures/bundled_17_clues.pl 2>/dev/null > tests/golden/arrange_bundled_17_fragment.json
	./crosswordsmith arrange --strict --candidates 3 --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null > tests/golden/arrange_bundled_17_candidates.json
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

# Performance ratchet: run the 15x15 ladder and diff each rung's search-inference
# count against benchmarks/baseline.json. A DROP is a win; a RISE past the
# baseline's tolerance is a regression (nonzero exit). Deterministic + machine-
# independent, so it hill-climbs the arrange algorithm (a -X% here predicts ~X%
# under WASM). wall/rss are reported but never gated. Core rungs only by default;
# add the hard tail with BENCH_ARGS=--heavy. NOT on the `make test` path.
#   make bench-check
#   make bench-check BENCH_ARGS=--heavy
bench-check:
	swipl -q benchmarks/check_baseline.pl $(BENCH_ARGS)

# Ratchet the baseline DOWN to the currently-measured numbers - run this to accept
# an improvement (or after an intentional algorithm change) and review the diff.
# Also appends the run to benchmarks/history.jsonl (the over-time ledger). Add
# BENCH_ARGS=--heavy to also re-record the heavy tail rungs.
#   make bench-record
#   make bench-record BENCH_ARGS=--heavy
bench-record:
	swipl -q benchmarks/check_baseline.pl --record $(BENCH_ARGS)

# Append the current measurement to benchmarks/history.jsonl WITHOUT moving the
# baseline - use it to log a data point (e.g. a mid-experiment reading, or the same
# commit on another host) while leaving the ratchet reference untouched. Stamps the
# git commit + timestamp so rungs stay comparable over time.
#   make bench-log
#   make bench-log BENCH_ARGS=--heavy
bench-log:
	swipl -q benchmarks/check_baseline.pl --log $(BENCH_ARGS)

# Render the recorded history (benchmarks/history.jsonl) as a per-rung trend:
# latest search_inf, the last step's delta, and the cumulative delta since the first
# entry. Reads the ledger only - runs no benchmark.
#   make bench-history
bench-history:
	swipl -q benchmarks/check_baseline.pl --history

# Strategy x fixture comparison matrix (CSV on stdout). Each fixture runs on
# its manifest grid (benchmarks/fixtures.pl). Optionally restrict strategies:
#   make bench-matrix BENCH_STRATEGIES="baseline mrv_capped"
# Compare on inferences (machine-independent); wall time is reporting-only.
BENCH_STRATEGIES ?=
bench-matrix:
	swipl -q benchmarks/run_matrix.pl -- $(BENCH_STRATEGIES)

# --- fill product bench + ratchet (campaign/fill Phase 0) ---------------------
# Product benchmark for `fill`: four attribution buckets per ladder rung -
# command (end-to-end CLI), dict_load (load_dict/3), grid (fill_grid/4), and
# search (fill_attempt/8, FRESH slots per sample) - over the rungs in
# benchmarks/fill_workloads.pl. Compare on search_inf/load_inf (deterministic,
# machine-independent); wall/rss are reporting-only. Core rungs by default.
#   make bench-fill
#   make bench-fill BENCH_ARGS=--heavy
#   make bench-fill BENCH_FORMAT=csv BENCH_ARGS="--fixture g11"
.PHONY: bench-fill bench-fill-check bench-fill-record bench-fill-log bench-fill-history
bench-fill:
	swipl -q benchmarks/run_fill.pl -- --format $(BENCH_FORMAT) $(BENCH_ARGS)

# Fill performance ratchet: diff each rung's search_inf against
# benchmarks/fill_baseline.json. search_inf GATES (a rise past tolerance fails,
# exit 1); load_inf is REPORTED but informational until Phase 3 decides
# otherwise. NOT on the `make test` path.
#   make bench-fill-check
#   make bench-fill-check BENCH_ARGS=--heavy
bench-fill-check:
	swipl -q benchmarks/check_fill_baseline.pl $(BENCH_ARGS)

# Ratchet the fill baseline to the currently-measured numbers (accept a win /
# intentional change) and append to benchmarks/fill_history.jsonl. After ANY
# record, read fill_baseline.json back and verify every rung is present.
#   make bench-fill-record BENCH_ARGS=--heavy
bench-fill-record:
	swipl -q benchmarks/check_fill_baseline.pl --record $(BENCH_ARGS)

# Append the current fill measurement to fill_history.jsonl WITHOUT moving the
# baseline.
bench-fill-log:
	swipl -q benchmarks/check_fill_baseline.pl --log $(BENCH_ARGS)

# Render the recorded fill history as a per-rung trend (reads the ledger only).
bench-fill-history:
	swipl -q benchmarks/check_fill_baseline.pl --history
