#!/usr/bin/env bash
# golden_parity.sh - prove the browser entry produces byte-identical output to the
# CLI's `arrange` verb for the same request. Native only: no wasm build, no
# browser, no server. This locks the "the web app really uses the arrange verb"
# contract as a regression test - the CLI (real binary) is the golden reference.
#
# For each case we feed the SAME clue set two ways: a JSON --input file to the CLI
# (its load_clues path + flag resolution) and an equivalent JSON payload to the
# browser entry (doc_to_words + solve_browser_json). Both dispatch the identical
# arrange_solve/4, so byte-identical output ⇔ the two input adapters agree.
#
# Covers what the spike actually exposes: fixed grid + strict/best-effort. (max
# size / seed / fragment are not wired in the spike - see plan §9.1.)
#
# Run: wasm/test/golden_parity.sh   (exit 0 = all cases byte-identical)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

CLUES='[{"answer":"CAT"},{"answer":"CAR"},{"answer":"ARC"},{"answer":"RAT"},{"answer":"TAR"}]'
printf '{"clues":%s}\n' "$CLUES" > "$TMP/clues.json"

fail=0
# check NAME SIZE CLI_FLAG BEST_EFFORT_BOOL
check() {
  local name="$1" size="$2" flag="$3" best="$4"
  printf '{"clues":%s,"size":%s,"mode":"fixed","bestEffort":%s}' "$CLUES" "$size" "$best" > "$TMP/payload.json"
  ./crosswordsmith arrange --input "$TMP/clues.json" --size "$size" "$flag" > "$TMP/golden.json" 2>"$TMP/cli.err" || {
    echo "FAIL  $name — CLI errored"; cat "$TMP/cli.err"; fail=1; return; }
  swipl -q "$HERE/golden_parity.pl" -- "$TMP/payload.json" > "$TMP/browser.json" 2>"$TMP/br.err" || {
    echo "FAIL  $name — browser entry errored"; cat "$TMP/br.err"; fail=1; return; }
  if diff -u "$TMP/golden.json" "$TMP/browser.json" > "$TMP/diff.txt"; then
    echo "PASS  $name (size $size $flag) — $(wc -c < "$TMP/golden.json") bytes identical"
  else
    echo "FAIL  $name (size $size $flag) — outputs differ:"; head -40 "$TMP/diff.txt"; fail=1
  fi
}

check toy-strict       5 --strict      false
check toy-best-effort  5 --best-effort true

[ "$fail" = 0 ] && echo "golden-parity OK" || echo "golden-parity FAILED"
exit "$fail"
