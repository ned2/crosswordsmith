#!/usr/bin/env bash
#
# run_tests.sh - run the full crosswordsmith test suite.
#
# Two layers:
#   1. plunit unit/integration tests (tests/crossword.plt, tests/arrange.plt)
#      via tests/run_tests.pl
#   2. golden-output regression tests: deterministic runs whose stdout must
#      match tests/golden/ (the legacy crossword CLI + the arrange engine).
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

echo
echo "=== golden output regression ==="
check_golden "crossword 17 topleft_across" \
    tests/golden/grid_17_topleft_across.txt \
    ./crossword.pl --input fixtures/bundled_17_clues.pl 17 topleft_across
check_golden "arrange bundled_17 fixed" \
    tests/golden/arrange_bundled_17_fixed.json \
    swipl -q -g 'arrange_run("fixtures/bundled_17_clues.pl",17,fixed),halt' -t 'halt(1)' arrange.pl
check_golden "arrange toc_demo max" \
    tests/golden/arrange_toc_demo_max.json \
    swipl -q -g 'arrange_run("fixtures/toc_demo.pl",25,max),halt' -t 'halt(1)' arrange.pl

echo
if [ "$status" -eq 0 ]; then
    echo "ALL TESTS PASSED"
else
    echo "TESTS FAILED"
fi
exit "$status"
