#!/usr/bin/env bash
# Fetch the pinned Collaborative Word List snapshot and emit a crosswordsmith
# scored dict (word;score) at a chosen --min-score floor.
#
# The pin is DP-9's (design-spec §10): a fixed upstream commit + sha256, so
# the output is a pure function of the flags — an unpinned fetch would break
# the determinism claims. The default invocation regenerates the bundled
# dicts/cwl50.dict byte-identically (see dicts/README.md for both hashes).
#
# Usage: scripts/fetch-cwl.sh [--min-score N] [--out FILE]
#   --min-score N   score floor to bake in (default 50)
#   --out FILE      output path (default dicts/cwl50.dict)
#
# Engine capacity envelope (measured, benchmarks/fill_quality/README.md):
#   N >= 50 : inside both envelopes (the bundled artifact).
#   30-49   : loads (>=30 is 437,400 words) but §8.4c search can blow the
#             global stack on grids with full-length slots (blocked_13a
#             crashed at ~17s) — usable on small/mid grids only.
#   1-29    : approaches the untested load ceiling (the full 566,665-word
#             cleaned list crashes dict load at SWI's default 1GB stack).
set -euo pipefail

PIN_COMMIT=2efe76e11ef311315e76d59700752733d69733d7   # xwordlist.dict last changed 2023-02-12
RAW_SHA256=a945a839a5f1e6f48caf9c8de446e5cd85f3567d7f62afcf54c6b738e8906ff4
RAW_URL="https://raw.githubusercontent.com/Crossword-Nexus/collaborative-word-list/$PIN_COMMIT/xwordlist.dict"

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MIN_SCORE=50
OUT="$ROOT/dicts/cwl50.dict"
while [ $# -gt 0 ]; do
  case "$1" in
    --min-score) MIN_SCORE="$2"; shift 2 ;;
    --out)       OUT="$2"; shift 2 ;;
    *) echo "usage: $0 [--min-score N] [--out FILE]" >&2; exit 2 ;;
  esac
done

if [ "$MIN_SCORE" -lt 50 ]; then
  echo "fetch-cwl: WARNING: floors below 50 exceed measured engine capacity" >&2
  echo "fetch-cwl: on some grids (see the envelope note in this script's header)" >&2
fi

TMP="$(mktemp)"; trap 'rm -f "$TMP"' EXIT
echo "fetch-cwl: fetching pinned snapshot $PIN_COMMIT ..." >&2
curl -sfL "$RAW_URL" -o "$TMP"
echo "$RAW_SHA256  $TMP" | sha256sum -c --quiet - \
  || { echo "fetch-cwl: sha256 MISMATCH — refusing to emit (pin drift or corrupt fetch)" >&2; exit 1; }

# Digit-bearing words (992 lines, e.g. MP3J;100) fold to wrong short entries
# under the engine's A-Z squeeze — always filtered (DP-9).
awk -F';' -v t="$MIN_SCORE" '$1 !~ /[0-9]/ && ($2 + 0) >= (t + 0)' "$TMP" > "$OUT"
echo "fetch-cwl: wrote $OUT ($(wc -l < "$OUT") words, score >= $MIN_SCORE)" >&2
