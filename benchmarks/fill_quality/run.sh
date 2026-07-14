#!/usr/bin/env bash
# Prototype fill-quality benchmark driver: crosswordsmith `fill` (scoreless MRV)
# vs ingrid_core (scored CSP) on identical grids + the same word list.
#
# NOT part of `make bench` — it needs two external, non-bundled dependencies:
#   1. ingrid_core   (`cargo install ingrid_core`; needs a Rust toolchain)
#   2. Spread the Wordlist scored list, as $STW (word;score, one per line).
#      Download from https://www.spreadthewordlist.com/ (CC BY-NC-SA 4.0).
#
# Usage: STW=/path/to/spread-the-wordlist.txt ./run.sh [workdir]
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CS="$ROOT/crosswordsmith"
WORK="${1:-$HERE/work}"
: "${STW:?set STW to the Spread-the-Wordlist word;score file}"
command -v ingrid_core >/dev/null || { echo "ingrid_core not on PATH (cargo install ingrid_core)"; exit 1; }

mkdir -p "$WORK"; cd "$WORK"
cut -d';' -f1 "$STW" > stw_plain.txt                 # scoreless dict for crosswordsmith
awk -F';' '$2>=50{print $1}' "$STW" > stw_ge50.txt   # score>=50 prefiltered dict
python3 "$HERE/gen_grids.py" .

for g in open4 open5 mini7 mini9; do
  "$CS" fill --grid "$g.json" --dict stw_plain.txt --out "cs_$g.json"  >/dev/null 2>&1 || echo "cs full  $g: no fill within budget"
  "$CS" fill --grid "$g.json" --dict stw_ge50.txt  --out "cs50_$g.json" >/dev/null 2>&1 || echo "cs >=50  $g: no fill within budget"
  ingrid_core --wordlist "$STW" --min-score 50 "$g.txt" > "ingrid_$g.txt" 2>/dev/null   || echo "ingrid   $g: no fill"
done

cp "$HERE/score_fill.py" .
python3 score_fill.py "$STW" open4 open5 mini7 mini9
