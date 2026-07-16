#!/usr/bin/env bash
# B0-I benchmark-only runner. STW words remain external and are never copied.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
PROBE="$HERE/probe_mac_b0.pl"
: "${STW:=/tmp/opencode/stw-lowercase-crossfire.txt}"

run_probe() {
    swipl -q -l "$ROOT/load.pl" -l "$PROBE" -g probe_mac_b0:run_cli -t halt -- "$@"
}

case "${1:-}" in
compare)
    run_probe compare "$2" "$3" "$4" "${5:-none}" "${6:-2000000000}"
    ;;
run)
    run_probe run "$2" "$3" "$4" "$5" "${6:-none}" "${7:-800000000}" "${8:-0}" "${9:-counters}"
    ;;
profile)
    run_probe profile "$2" "$3" "$4" "$5" "${6:-none}" "${7:-800000000}"
    ;;
*)
    echo "usage: $0 compare GRID DICT MIN [SEEDS] [BUDGET]" >&2
    echo "       $0 run LABEL GRID DICT MIN [SEEDS] [BUDGET] [TIMEOUT] [off|counters|timing]" >&2
    echo "       $0 profile LABEL GRID DICT MIN [SEEDS] [BUDGET]" >&2
    exit 2
    ;;
esac
