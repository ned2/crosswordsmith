#!/usr/bin/env bash
#
# scripts/xword-parity.sh - engine<->xword export parity probe.
#
# Runs the SWI-Prolog engine (./crosswordsmith export) and the Python
# companion (xword convert) over the SAME native layout JSON and byte-diffs
# the ipuz and exolve outputs. One command instead of hand-rolled
# export/convert/diff probes.
#
# Byte-parity is BEST-EFFORT policy, not a contract (docs/xword-spec.md §14):
# engine JSON goes through SWI json_write_dict (type-dependent colon spacing)
# while xword uses json.dumps. Exit 0 = byte-identical on every format;
# exit 1 = differences (diff shown) - treat that as information, not a broken
# build. Exit 2 = a toolchain failed outright.
#
#   scripts/xword-parity.sh [layout.json]
#   make xword-parity [XWORD_PARITY_FIXTURE=<layout.json>]
#
# Default fixture: tests/golden/arrange_bundled_17_fixed.json

set -u
cd "$(dirname "$0")/.."

fixture=${1:-tests/golden/arrange_bundled_17_fixed.json}
if [[ ! -f $fixture ]]; then
    echo "xword-parity: no such layout: $fixture" >&2
    exit 2
fi
fixture_abs=$(realpath "$fixture")

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

status=0
for fmt in ipuz exolve; do
    if ! ./crosswordsmith export --to "$fmt" "$fixture_abs" \
            >"$tmp/engine.$fmt" 2>"$tmp/engine.$fmt.err"; then
        echo "xword-parity: engine export --to $fmt failed:" >&2
        cat "$tmp/engine.$fmt.err" >&2
        exit 2
    fi
    if ! (cd xword && uv run --quiet xword convert --to "$fmt" "$fixture_abs") \
            >"$tmp/xword.$fmt" 2>"$tmp/xword.$fmt.err"; then
        echo "xword-parity: xword convert --to $fmt failed:" >&2
        cat "$tmp/xword.$fmt.err" >&2
        exit 2
    fi
    if diff -u --label "engine.$fmt" --label "xword.$fmt" \
            "$tmp/engine.$fmt" "$tmp/xword.$fmt" >"$tmp/diff.$fmt"; then
        echo "parity ($fmt): byte-identical"
    else
        changed=$(grep -c '^[+-][^+-]' "$tmp/diff.$fmt" || true)
        echo "parity ($fmt): DIFFERS ($changed changed lines; best-effort per xword-spec §14)"
        cat "$tmp/diff.$fmt"
        status=1
    fi
done
exit $status
