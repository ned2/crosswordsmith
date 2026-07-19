# Convenience targets for crosswordsmith.

SHELL := /bin/bash

.PHONY: test test-wasm test-xword xword-parity unit golden update-golden fuzz bench bench-check bench-exact bench-record bench-log bench-history bench-matrix bench-arrange-verify bench-arrange-promote bench-greedy bench-greedy-check bench-greedy-exact bench-greedy-record bench-greedy-log bench-greedy-history bench-greedy-identity bench-greedy-verify bench-greedy-promote

BENCH_FORMAT ?= text
BENCH_ARGS ?=

# Full suite: plunit tests + golden-output regression.
test:
	./run_tests.sh

# One-command wasm battery (wasm/test/run_all.sh): stage swipl-web artifacts +
# qcompile the app qlf with the WASM swipl (both skipped when fresh), then
# value-golden, the type lock, and the headless-Chrome probes behind a
# self-managed static server on a free port. Prereqs: the swipl-devel wasm
# build tree (wasm/build/build-wasm.sh; override WASM_BUILD/WASM_SWIPL) and a
# one-time `( cd wasm/test && PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 npm install )`.
# The ~60s gate #1 probe is opt-in:
#   make test-wasm
#   make test-wasm WASM_TEST_ARGS=--yield
WASM_TEST_ARGS ?=
test-wasm:
	wasm/test/run_all.sh $(WASM_TEST_ARGS)

# xword companion (Python conversion toolchain) tests. NOT part of `make test`.
test-xword:
	cd xword && uv run pytest -q

# Engine<->xword export parity probe: run both toolchains over one native
# layout and byte-diff the ipuz + exolve outputs. Byte-parity is BEST-EFFORT
# policy, not a contract (docs/xword-spec.md §14) - a nonzero exit is
# information, not a broken build. NOT part of `make test`.
#   make xword-parity
#   make xword-parity XWORD_PARITY_FIXTURE=fixtures/titled_layout.json
XWORD_PARITY_FIXTURE ?= tests/golden/arrange_bundled_17_fixed.json
xword-parity:
	scripts/xword-parity.sh $(XWORD_PARITY_FIXTURE)

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
	@diff -u tests/golden/arrange_bundled_17_fragment.json \
		<(./crosswordsmith arrange --strict --fragment fixtures/bundled_17_fragment_thin.json --size 17 --input fixtures/bundled_17_clues.pl 2>/dev/null) \
		&& echo "golden (arrange fragment, thin form): OK"
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
	@diff -u tests/golden/export_titled.ipuz \
		<(./crosswordsmith export --to ipuz fixtures/titled_layout.json 2>/dev/null) \
		&& echo "golden (export ipuz, titled+authored): OK"
	@diff -u tests/golden/export_titled.exolve \
		<(./crosswordsmith export --to exolve fixtures/titled_layout.json 2>/dev/null) \
		&& echo "golden (export exolve, titled+authored): OK"
	@diff -u tests/golden/fill_3.json \
		<(./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt 2>/dev/null) \
		&& echo "golden (fill 3x3): OK"
	@diff -u tests/golden/fill_3_seeded.json \
		<(./crosswordsmith fill --grid fixtures/fill_grid_split3.json --seeds fixtures/fill_seed_cow_top.json --dict fixtures/dict_cow_pig.txt 2>/dev/null) \
		&& echo "golden (fill seeded, canonical): OK"
	@diff -u tests/golden/fill_3_seeded.json \
		<(./crosswordsmith fill --grid fixtures/fill_grid_split3.json --seeds fixtures/fill_seed_cow_top_thin.json --dict fixtures/dict_cow_pig.txt 2>/dev/null) \
		&& echo "golden (fill seeded, thin form): OK"
	@diff -u tests/golden/fill_3_rng_seed7.json \
		<(./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt --seed 7 2>/dev/null) \
		&& echo "golden (fill --seed 7): OK"
	@diff -u tests/golden/fill_scored.json \
		<(./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/dict_scored_sample.txt 2>/dev/null) \
		&& echo "golden (fill scored, default prune): OK"
	@tmp=$$(mktemp) && \
		diff -u tests/golden/fill_scored_min50.json \
			<(./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/dict_scored_sample.txt --min-score 50 --report-json $$tmp 2>/dev/null) \
		&& echo "golden (fill scored, min-score 50): OK" \
		&& diff -u tests/golden/fill_scored_min50_report.json $$tmp \
		&& echo "golden (fill scored, report json): OK"; \
		rc=$$?; rm -f $$tmp; exit $$rc

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
	./crosswordsmith export --to ipuz fixtures/titled_layout.json 2>/dev/null > tests/golden/export_titled.ipuz
	./crosswordsmith export --to exolve fixtures/titled_layout.json 2>/dev/null > tests/golden/export_titled.exolve
	./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt 2>/dev/null > tests/golden/fill_3.json
	./crosswordsmith fill --grid fixtures/fill_grid_split3.json --seeds fixtures/fill_seed_cow_top.json --dict fixtures/dict_cow_pig.txt 2>/dev/null > tests/golden/fill_3_seeded.json
	./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt --seed 7 2>/dev/null > tests/golden/fill_3_rng_seed7.json
	./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/dict_scored_sample.txt 2>/dev/null > tests/golden/fill_scored.json
	./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/dict_scored_sample.txt --min-score 50 --report-json tests/golden/fill_scored_min50_report.json 2>/dev/null > tests/golden/fill_scored_min50.json
	./crosswordsmith fill --grid fixtures/fill_grid_15a.json --dict fixtures/dict/enable1.txt 2>/dev/null > tests/golden/fill_15_bench.json
	./crosswordsmith fill --grid grids/amer11.json --dict dicts/cwl50.dict 2>/dev/null > tests/golden/fill_cwl_amer11.json
	@echo "Regenerated golden files (arrange fixed/max/fragment/candidates + lint toc + export ipuz/exolve + export titled ipuz/exolve + fill plain/seeded/scored/15-bench/cwl-amer11)"

# Product benchmark for `arrange`: end-to-end command latency, the in-process
# search alone, and the CLI-wrapper overhead between them (rest = command -
# search), over the workloads in benchmarks/workloads.pl. Core workloads only by
# default; results are machine-specific and reporting-only. Search inference
# counts are regression signals only against a same-SWI baseline.
#   make bench
#   make bench BENCH_ARGS=--heavy                       # + gated tail and one ~26s latency probe
#   make bench BENCH_FORMAT=csv BENCH_ARGS="--fixture real"
bench:
	swipl -q benchmarks/run_arrange.pl -- --format $(BENCH_FORMAT) $(BENCH_ARGS)

# Performance ratchet: run the 9x9/13x13/15x15/21x21 strict workload set and
# diff each rung's search-inference count against benchmarks/baseline.json. A
# DROP is a win; a RISE past the baseline's tolerance is a regression (nonzero
# exit). Same-SWI inference deltas hill-climb the arrange algorithm; wall/rss are
# reported but never gated. Core rungs only by default;
# add the hard tail with BENCH_ARGS=--heavy. NOT on the `make test` path.
#   make bench-check
#   make bench-check BENCH_ARGS=--heavy
bench-check:
	swipl -q benchmarks/check_baseline.pl $(BENCH_ARGS)

# Refactor gate: require exact same-SWI inference counts over every core and
# heavy rung. Unlike bench-check, both increases and decreases fail.
bench-exact:
	swipl -q benchmarks/check_baseline.pl --exact

# One-command strict arrange adjudication: full native tests, diagnostics-bearing
# ladder identity, then the inference ratchet. Selection mirrors bench-check:
# core by default, core+heavy with BENCH_ARGS=--heavy.
bench-arrange-verify:
	./run_tests.sh
	benchmarks/check_arrange_identity.sh $(BENCH_ARGS)
	swipl -q benchmarks/check_baseline.pl $(BENCH_ARGS)

# Check and promote ONE strict-ladder measurement. The recorder rejects
# regressions before writing, appends history only on success, then reads the
# baseline back and verifies every measured rung/value persisted.
bench-arrange-promote:
	swipl -q benchmarks/check_baseline.pl --promote $(BENCH_ARGS)

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

# Independent greedy best-effort/candidates substrate. Construction, the two
# directly searched blocks plus transpose partners, postprocess, and command
# layers are separate; sweep inferences are primary. This ratchet never reads or
# writes baseline.json or history.jsonl.
bench-greedy:
	swipl -q benchmarks/run_arrange_greedy.pl -- --format $(BENCH_FORMAT) $(BENCH_ARGS)

bench-greedy-check:
	swipl -q benchmarks/check_greedy_baseline.pl $(BENCH_ARGS)

bench-greedy-exact:
	swipl -q benchmarks/check_greedy_baseline.pl --exact

bench-greedy-record:
	swipl -q benchmarks/check_greedy_baseline.pl --record $(BENCH_ARGS)

bench-greedy-log:
	swipl -q benchmarks/check_greedy_baseline.pl --log $(BENCH_ARGS)

bench-greedy-history:
	swipl -q benchmarks/check_greedy_baseline.pl --history

bench-greedy-identity:
	benchmarks/check_greedy_identity.sh

bench-greedy-verify:
	./run_tests.sh
	benchmarks/check_greedy_identity.sh
	swipl -q benchmarks/check_greedy_baseline.pl $(BENCH_ARGS)

# Check and record the exact same one measurement, with atomic write plus
# read-back completeness verification before replacing the baseline.
bench-greedy-promote:
	swipl -q benchmarks/check_greedy_baseline.pl --promote $(BENCH_ARGS)

# Strategy x fixture comparison matrix (CSV on stdout). Each fixture runs on
# its manifest grid (benchmarks/fixtures.pl). Optionally restrict strategies:
#   make bench-matrix BENCH_STRATEGIES="baseline mrv_capped"
# Compare inference counts only within one SWI version; wall time is reporting-only.
BENCH_STRATEGIES ?=
bench-matrix:
	swipl -q benchmarks/run_matrix.pl -- $(BENCH_STRATEGIES)

# --- browser payload/startup bench (payload plan Phase 0) ---------------------
# Payload/request/startup benchmark for the browser bundle
# (docs/plans/wasm-payload-performance.md §4): artifact raw/gzip/brotli sizes,
# per-scenario cold request counts through a cache-aware counting server, and
# fresh-Chrome-per-sample readiness/first-arrange timing. Wall-clock numbers
# are machine-specific and reporting-only; only the deterministic artifact
# sizes gate (bench-wasm-payload-check). Prereqs: staged + stamped wasm/client/
# (run `make test-wasm` or the build first) and the wasm/test npm toolchain.
#   make bench-wasm-payload                             # measure + report
#   make bench-wasm-payload BENCH_ARGS="--samples 20"
#   make bench-wasm-payload-record                      # + rewrite the committed baseline
#   make bench-wasm-payload-check                       # sizes vs baseline (no browser)
.PHONY: bench-wasm-payload bench-wasm-payload-record bench-wasm-payload-check
bench-wasm-payload:
	node wasm/test/payload_bench.mjs $(BENCH_ARGS)

bench-wasm-payload-record:
	node wasm/test/payload_bench.mjs --record $(BENCH_ARGS)

bench-wasm-payload-check:
	node wasm/test/payload_bench.mjs --check

# --- fill product bench + ratchet ---------------------------------------------
# Product benchmark for `fill`: four attribution buckets per ladder rung -
# command (end-to-end CLI), dict_load (load_dict/3), grid (fill_grid/4), and
# search (fill_attempt/8, FRESH slots per sample) - over the rungs in
# benchmarks/fill_workloads.pl. Compare search_inf/load_inf only against a
# same-SWI baseline; wall/rss are reporting-only. Core rungs by default.
#   make bench-fill
#   make bench-fill BENCH_ARGS=--heavy
#   make bench-fill BENCH_FORMAT=csv BENCH_ARGS="--fixture g11"
.PHONY: bench-fill bench-fill-check bench-fill-exact bench-fill-record bench-fill-log bench-fill-history bench-fill-verify bench-fill-promote bench-fill-quality-test bench-fill-quality-check
bench-fill:
	swipl -q benchmarks/run_fill.pl -- --format $(BENCH_FORMAT) $(BENCH_ARGS)

# Fill performance ratchet: diff each rung's search_inf and load_inf against
# benchmarks/fill_baseline.json. Either metric rising past tolerance fails with
# exit 1. NOT on the `make test` path.
#   make bench-fill-check
#   make bench-fill-check BENCH_ARGS=--heavy
bench-fill-check:
	swipl -q benchmarks/check_fill_baseline.pl $(BENCH_ARGS)

bench-fill-exact:
	swipl -q benchmarks/check_fill_baseline.pl --exact

# Ratchet the fill baseline to the currently-measured numbers (accept a win /
# intentional change) and append to benchmarks/fill_history.jsonl. After ANY
# record, read fill_baseline.json back and verify every rung is present.
#   make bench-fill-record BENCH_ARGS=--heavy
bench-fill-record:
	swipl -q benchmarks/check_fill_baseline.pl --record $(BENCH_ARGS)

# One-command fill adjudication battery: full test suite, CLI byte-identity,
# artifact byte-identity, then the inference ratchet - the four checks every
# fill-experiment accept runs, in one invocation (fail-fast, in that order).
#   make bench-fill-verify
#   make bench-fill-verify BENCH_ARGS=--heavy
bench-fill-verify:
	./run_tests.sh
	benchmarks/check_fill_identity.sh
	benchmarks/check_fill_identity_artifact.sh
	swipl -q benchmarks/check_fill_baseline.pl $(BENCH_ARGS)

# Check-then-ratchet from ONE measurement: run the ladder once, diff against
# the baseline, and only if clean (zero regressions) record that same run as
# the new baseline + history entry. Replaces the bench-fill-check ->
# bench-fill-record double measurement when accepting a win.
#   make bench-fill-promote
#   make bench-fill-promote BENCH_ARGS=--heavy
bench-fill-promote:
	swipl -q benchmarks/check_fill_baseline.pl --promote $(BENCH_ARGS)

# Append the current fill measurement to fill_history.jsonl WITHOUT moving the
# baseline.
bench-fill-log:
	swipl -q benchmarks/check_fill_baseline.pl --log $(BENCH_ARGS)

# Render the recorded fill history as a per-rung trend (reads the ledger only).
bench-fill-history:
	swipl -q benchmarks/check_fill_baseline.pl --history

# Python-only checks for the retained independent quality scorer and gate.
bench-fill-quality-test:
	python3 -m unittest discover -s benchmarks/fill_quality -p 'test_*.py'

# AC-FILL-12 external-dictionary gate. This is intentionally not part of
# `make test`: Spread the Wordlist is CC BY-NC-SA and is never bundled.
#   make bench-fill-quality-check STW=/path/to/spread-the-wordlist.txt
bench-fill-quality-check: bench-fill-quality-test
	python3 benchmarks/fill_quality/check_ac_fill_12.py --dict "$(STW)"
