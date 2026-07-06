#!/usr/bin/env bash
#
# run_tests.sh - run the full crosswordsmith test suite.
#
# Three layers:
#   1. plunit unit/integration tests (tests/*.plt) via tests/run_tests.pl
#   2. golden-output regression: deterministic crosswordsmith CLI runs whose
#      stdout must match tests/golden/ byte-for-byte
#   3. CLI exit-code contract checks (e.g. AC-LINT-2: FAIL -> non-zero)
#
# Exits non-zero if any layer fails. Run from anywhere.

set -u
cd "$(dirname "$0")"

status=0

echo "=== plunit tests ==="
if swipl -q tests/run_tests.pl; then
    echo "plunit: OK"
else
    echo "plunit: FAILED"
    status=1
fi

# check_golden <name> <golden-file> -- <command...>: run the command, compare
# its stdout (stderr discarded) byte-for-byte against the golden file.
check_golden() {
    local name="$1" golden="$2"; shift 2
    local actual; actual="$(mktemp)"
    "$@" 2>/dev/null > "$actual"
    if diff -u "$golden" "$actual" >/dev/null; then
        echo "golden ($name): OK"
    else
        echo "golden ($name): FAILED (output differs from $golden)"
        diff -u "$golden" "$actual" | head -40
        status=1
    fi
    rm -f "$actual"
}

# check_exit <name> <expected-code> <command...>: run the command (stdout +
# stderr discarded) and assert its exit status. Covers the AC-LINT-2 exit-code
# contract (FAIL -> non-zero), which check_golden cannot (it discards the code).
check_exit() {
    local name="$1" want="$2"; shift 2
    "$@" >/dev/null 2>&1
    local got=$?
    if [ "$got" -eq "$want" ]; then
        echo "exit ($name): OK (exit $got)"
    else
        echo "exit ($name): FAILED (got exit $got, want $want)"
        status=1
    fi
}

# check_fail_report <name> <command...>: assert a clean fill-failure (exit 1,
# a `fill:`-framed stderr line, and NO raw internal throw). Distinguishes the
# clean report from an uncaught engine error, which check_exit (code only)
# cannot - both exit 1. The fill: match is unanchored because hooked seed
# errors render behind print_message's "ERROR: [Thread main] " prefix.
check_fail_report() {
    local name="$1"; shift
    local err; err="$("$@" 2>&1 >/dev/null)"; local got=$?
    if [ "$got" -eq 1 ] \
       && printf '%s' "$err" | grep -q 'fill:' \
       && ! printf '%s' "$err" | grep -qi 'Domain error\|unique_key_pairs'; then
        echo "fail-report ($name): OK"
    else
        echo "fail-report ($name): FAILED (exit $got)"; printf '%s\n' "$err" | head -3; status=1
    fi
}

echo
echo "=== golden output regression (the crosswordsmith CLI, end to end) ==="
check_golden "arrange fixed" \
    tests/golden/arrange_bundled_17_fixed.json \
    ./crosswordsmith arrange --strict --size 17 --input fixtures/bundled_17_clues.pl
check_golden "arrange max" \
    tests/golden/arrange_toc_demo_max.json \
    ./crosswordsmith arrange --strict --max-size 25 --input fixtures/toc_demo.pl
check_golden "arrange fragment" \
    tests/golden/arrange_bundled_17_fragment.json \
    ./crosswordsmith arrange --strict --fragment fixtures/bundled_17_fragment.json --input fixtures/bundled_17_clues.pl
# The thin spelling of the same fragment against the SAME golden file:
# byte-identity between the two forms is AC-FRAG-4 itself.
check_golden "arrange fragment (thin form)" \
    tests/golden/arrange_bundled_17_fragment.json \
    ./crosswordsmith arrange --strict --fragment fixtures/bundled_17_fragment_thin.json --size 17 --input fixtures/bundled_17_clues.pl
check_golden "arrange candidates" \
    tests/golden/arrange_bundled_17_candidates.json \
    ./crosswordsmith arrange --strict --candidates 3 --size 17 --input fixtures/bundled_17_clues.pl
check_golden "lint toc" \
    tests/golden/lint_bundled_17_toc.json \
    ./crosswordsmith lint --profile toc tests/golden/arrange_bundled_17_fixed.json
check_golden "export ipuz" \
    tests/golden/export_bundled_17.ipuz \
    ./crosswordsmith export --to ipuz tests/golden/arrange_bundled_17_fixed.json
check_golden "export exolve" \
    tests/golden/export_bundled_17.exolve \
    ./crosswordsmith export --to exolve tests/golden/arrange_bundled_17_fixed.json
check_golden "fill 3x3" \
    tests/golden/fill_3.json \
    ./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt
# Scale golden (fill bench Phase 0): a 15x15 stock grid filled from full ENABLE.
# Locks byte-identity of a realistic-scale fill (~10.7M search inferences, a few
# seconds CLI). Stdout byte-compared; stderr (the --verbose summary path) discarded.
check_golden "fill 15 bench" \
    tests/golden/fill_15_bench.json \
    ./crosswordsmith fill --grid fixtures/fill_grid_15a.json --dict fixtures/dict/enable1.txt

echo
echo "=== CLI exit-code contract (AC-LINT-2) ==="
# A FAIL-severity verdict under an enforcing profile must exit non-zero...
check_exit "lint blocked-uk FAIL -> nonzero" 1 \
    ./crosswordsmith lint --profile blocked-uk fixtures/lint_fail_layout.json
# ...while PASS/WARN (advisory toc never FAILs) exits zero.
check_exit "lint toc PASS/WARN -> zero" 0 \
    ./crosswordsmith lint --profile toc tests/golden/arrange_bundled_17_fixed.json
# A --no-<bool> negation (optparse auto-negates booleans) is an undocumented
# flag name and must be rejected, not silently accepted (AC-CLI-2, R9). Without
# the guard this run would be accepted (best_effort(false) -> still strict) -> 0.
check_exit "arrange --no- flag rejected" 1 \
    ./crosswordsmith arrange --no-best-effort --size 17 --input fixtures/bundled_17_clues.pl
# --size (exact NxN) and --max-size (ceiling + tight-square crop) are mutually
# exclusive framings; supplying both is a usage error, not a silent precedence.
check_exit "arrange --size + --max-size rejected" 1 \
    ./crosswordsmith arrange --size 17 --max-size 20 --input fixtures/bundled_17_clues.pl
# A seed answer the dictionary can also place must report cleanly, not throw a
# raw unique_key_pairs domain error (regression: fill-seed-pin-crash-fix).
check_fail_report "fill seed-reused unsolvable -> clean fail" \
    ./crosswordsmith fill --grid fixtures/fill_grid_split3.json \
        --seeds fixtures/fill_seed_cow_top.json --dict fixtures/dict_cow.txt
# Two identical seed answers must be rejected before searching (fill_seed_duplicate).
check_fail_report "fill duplicate seeds -> clean fail" \
    ./crosswordsmith fill --grid fixtures/fill_grid_split3.json \
        --seeds fixtures/fill_seed_cow_both.json --dict fixtures/dict_cow.txt

# check_stderr <name> <grep-pattern> <want: present|absent> <command...>: run
# the command (stdout discarded) and assert whether its stderr matches.
check_stderr() {
    local name="$1" pattern="$2" want="$3"; shift 3
    local err; err="$("$@" 2>&1 >/dev/null)"
    local matched=absent
    if printf '%s' "$err" | grep -q "$pattern"; then matched=present; fi
    if [ "$matched" = "$want" ]; then
        echo "stderr ($name): OK ($pattern $want)"
    else
        echo "stderr ($name): FAILED (want $pattern $want, got: $err)"
        status=1
    fi
}

# check_stdout <name> <grep-pattern> <want: present|absent> <command...>: same
# as check_stderr but against stdout (stderr discarded) — payload assertions.
check_stdout() {
    local name="$1" pattern="$2" want="$3"; shift 3
    local out; out="$("$@" 2>/dev/null)"
    local matched=absent
    if printf '%s' "$out" | grep -q "$pattern"; then matched=present; fi
    if [ "$matched" = "$want" ]; then
        echo "stdout ($name): OK ($pattern $want)"
    else
        echo "stdout ($name): FAILED (want $pattern $want)"
        status=1
    fi
}

echo
echo "=== stderr contract (design-spec §5.1: quiet success, --verbose summaries) ==="
# Default: a clean success prints NOTHING on stderr - quality caveats (cap
# inert, dropped words) ride the payload's diagnostics, not the terminal.
check_stderr "arrange quiet by default" "placed" absent \
    ./crosswordsmith arrange --size 17 --input fixtures/bundled_17_clues.pl
check_stderr "cap-inert off stderr" "cap inert" absent \
    ./crosswordsmith arrange --size 17 --input fixtures/bundled_17_clues.pl
# ...the cap-inert compromise is reported in the payload instead (§7.2, INV-3,
# json-output-spec §6.4), alongside the dropped set.
check_stdout "capInert in diagnostics" '"capInert":true' present \
    ./crosswordsmith arrange --size 17 --input fixtures/bundled_17_clues.pl
check_stdout "dropped set in diagnostics" '"OMEGA POINT"' present \
    ./crosswordsmith arrange --best-effort --size 5 --input fixtures/bundled_17_clues.pl
# fill's payload must NOT grow arrange diagnostics.
check_stdout "fill payload has no diagnostics" '"diagnostics"' absent \
    ./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt
# --verbose opts into the summary (cap-status note included).
check_stderr "arrange --verbose summary" "placed 6, reward 60 (cap inert" present \
    ./crosswordsmith arrange --verbose --size 17 --input fixtures/bundled_17_clues.pl
# fill: quiet by default, summary under --verbose.
check_stderr "fill quiet by default" "filled" absent \
    ./crosswordsmith fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt
check_stderr "fill --verbose summary" "filled" present \
    ./crosswordsmith fill --verbose --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt

echo
if [ "$status" -eq 0 ]; then
    echo "ALL TESTS PASSED"
else
    echo "TESTS FAILED"
fi
exit "$status"
