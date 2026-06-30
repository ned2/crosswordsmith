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

echo
echo "=== golden output regression (the crosswordsmith CLI, end to end) ==="
check_golden "arrange fixed" \
    tests/golden/arrange_bundled_17_fixed.json \
    ./crosswordsmith arrange --strict --size-mode fixed --size 17 --input fixtures/bundled_17_clues.pl
check_golden "arrange max" \
    tests/golden/arrange_toc_demo_max.json \
    ./crosswordsmith arrange --strict --size-mode max --size 25 --input fixtures/toc_demo.pl
check_golden "arrange fragment" \
    tests/golden/arrange_bundled_17_fragment.json \
    ./crosswordsmith arrange --strict --size-mode fixed --fragment fixtures/bundled_17_fragment.json --input fixtures/bundled_17_clues.pl
check_golden "arrange candidates" \
    tests/golden/arrange_bundled_17_candidates.json \
    ./crosswordsmith arrange --strict --size-mode fixed --candidates 3 --size 17 --input fixtures/bundled_17_clues.pl
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

echo
if [ "$status" -eq 0 ]; then
    echo "ALL TESTS PASSED"
else
    echo "TESTS FAILED"
fi
exit "$status"
