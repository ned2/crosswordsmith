#!/usr/bin/env bash
#
# run_tests.sh - run the full crosswordsmith test suite.
#
# Two layers:
#   1. plunit unit/integration tests (tests/crossword.plt) via tests/run_tests.pl
#   2. a golden-output regression test: the deterministic CLI run
#      `crossword.pl 17 topleft_across` must match tests/golden/.
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

echo
echo "=== golden output regression (17 topleft_across) ==="
golden="tests/golden/grid_17_topleft_across.txt"
actual="$(mktemp)"
trap 'rm -f "$actual"' EXIT
# Compare files directly (not via $(...), which strips trailing newlines).
./crossword.pl 17 topleft_across 2>/dev/null > "$actual"
if diff -u "$golden" "$actual" >/dev/null; then
    echo "golden: OK"
else
    echo "golden: FAILED (output differs from $golden)"
    diff -u "$golden" "$actual" | head -40
    status=1
fi

echo
if [ "$status" -eq 0 ]; then
    echo "ALL TESTS PASSED"
else
    echo "TESTS FAILED"
fi
exit "$status"
