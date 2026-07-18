#!/usr/bin/env bash
# Optional fill-quality comparison: crosswordsmith vs ingrid_core on identical
# grids and the same scored dictionary. This reports comparisons; the permanent
# AC-FILL-12 gate is `make bench-fill-quality-check`.
#
# NOT part of `make test` — it needs two external, non-bundled dependencies:
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
[[ -r "$STW" ]] || { echo "STW is not readable: $STW"; exit 1; }
STW="$(realpath "$STW")"
command -v ingrid_core >/dev/null || { echo "ingrid_core not on PATH (cargo install ingrid_core)"; exit 1; }

mkdir -p "$WORK"; cd "$WORK"
cut -d';' -f1 "$STW" > stw_plain.txt                 # scoreless dict for crosswordsmith
awk -F';' '$2>=50{print $1}' "$STW" > stw_ge50.txt   # score>=50 prefiltered dict
python3 "$HERE/gen_grids.py" .

grids=(open4 open5 mini7 mini9 amer11)
for g in "${grids[@]}"; do
  rm -f "cs_$g.json" "cs50_$g.json" "cs_minscore_$g.json" \
      "cs_minscore_$g.report.json" "ingrid_$g.txt"
  # A non-zero exit can be budget exhaustion, a capacity report, or bad input.
  # Rerun without the stderr redirect to distinguish them.
  "$CS" fill --grid "$g.json" --dict stw_plain.txt --out "cs_$g.json"  >/dev/null 2>&1 || echo "cs full  $g: no fill (budget or dict-load error)"
  "$CS" fill --grid "$g.json" --dict stw_ge50.txt  --out "cs50_$g.json" >/dev/null 2>&1 || echo "cs >=50  $g: no fill (budget or dict-load error)"
  # Native scored fill (§8.4a): STW ingested directly as word;score, hard
  # prune at the clean floor. The FS-1 acceptance gate: this column MUST
  # match the >=50-prefiltered-dict column (mean/min 50, 0 below-clean),
  # and its --report-json numbers must agree with score_fill.py's post-hoc
  # stats on the same fill (checked by score_fill.py below).
  "$CS" fill --grid "$g.json" --dict "$STW" --min-score 50 \
      --report-json "cs_minscore_$g.report.json" --out "cs_minscore_$g.json" >/dev/null 2>&1 \
      || echo "cs min50 $g: no fill (budget or dict-load error)"
  ingrid_core --wordlist "$STW" --min-score 50 "$g.txt" > "ingrid_$g.txt" 2>/dev/null   || echo "ingrid   $g: no fill"
done

python3 "$HERE/score_fill.py" "$STW" "${grids[@]}"
