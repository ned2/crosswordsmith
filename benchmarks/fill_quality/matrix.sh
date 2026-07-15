#!/usr/bin/env bash
# FS-3(b): the completion x min-score frontier — the SEARCH-POWER axis.
#
# For every mask (the four easy quality grids, the American 11x11, and the
# three bundled blocked stock grids) x every --min-score threshold, record
# whether crosswordsmith's fixed-budget MRV completes and whether ingrid_core
# does (90s wall cap, the README's timeout convention). Quality per completed
# native fill rides along for free from --report-json.
#
# This measures where the budget ceiling bites and whether pruning moves it —
# evidence for the FS-4 search-power decision pass. It never gates: unlike
# run.sh (the FS-1 quality gate, which exits non-zero on report-json
# disagreement), an incomplete fill here is a DATA POINT, not a failure.
#
# Same external deps as run.sh: ingrid_core on PATH + a word;score list as
# $STW. Usage: STW=/path/to/list.txt ./matrix.sh [workdir]
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CS="$ROOT/crosswordsmith"
WORK="${1:-$HERE/work_matrix}"
: "${STW:?set STW to a word;score wordlist file}"
command -v ingrid_core >/dev/null || { echo "ingrid_core not on PATH (cargo install ingrid_core)"; exit 1; }

GRIDS="open4 open5 mini7 mini9 amer11 blocked_13a blocked_13b blocked_15a"
THRESHOLDS="1 30 50"
INGRID_CAP=90   # seconds; matches the README's ">90s = timeout" convention

mkdir -p "$WORK"; cd "$WORK"
python3 "$HERE/gen_grids.py" .
for b in blocked_13a blocked_13b blocked_15a; do
  cp "$ROOT/grids/$b.json" "$b.json"
  python3 -c "
import json, sys
mask = json.load(open('$b.json'))['mask']
open('$b.txt', 'w').write('\n'.join(mask) + '\n')
"
done

printf '%-12s %9s  %-34s %s\n' grid minscore crosswordsmith ingrid_core
printf '%s\n' "-----------------------------------------------------------------------------"
for g in $GRIDS; do
  for t in $THRESHOLDS; do
    if "$CS" fill --grid "$g.json" --dict "$STW" --min-score "$t" \
        --report-json "cs_${g}_ms${t}.report.json" \
        --out "cs_${g}_ms${t}.json" >/dev/null 2>&1; then
      cs=$(python3 -c "
import json
r = json.load(open('cs_${g}_ms${t}.report.json'))
print(f\"ok    mean={r['mean']} min={r['min']} below50={r['belowThreshold']}\")
")
    else
      cs="NOT-COMPLETED"
    fi
    if timeout "$INGRID_CAP" ingrid_core --wordlist "$STW" --min-score "$t" \
        "$g.txt" > "ingrid_${g}_ms${t}.txt" 2>/dev/null; then
      ig="ok"
    else
      ig="NOT-COMPLETED(<=${INGRID_CAP}s)"
    fi
    printf '%-12s %9s  %-34s %s\n' "$g" "$t" "$cs" "$ig"
  done
done
