#!/usr/bin/env bash
#
# determinism_fuzz.sh - fuzz INV-2 / AC-EMIT-1 / AC-X-2 (deterministic,
# byte-identical output) and degenerate-input robustness (INV-3) across a broad
# matrix of verbs x flags x inputs, including pathological inputs (empty,
# single-letter, unicode, duplicate, hyphen/space/punct/digit, lowercase).
#
# For each case it runs the CLI THREE times as separate processes and asserts:
#   - identical exit code across runs,
#   - byte-identical stdout across runs (the INV-2 contract),
#   - a DEFINED exit (0 or 1: success or a clean reported failure - never a
#     hang, signal, or other code),
# plus targeted --out partial-write checks (spec §5.2).
#
# This is an ON-DEMAND harness (many subprocess runs); it is not part of the
# default `make test`. Run via `make fuzz` or `tests/determinism_fuzz.sh`.
# Deterministic itself: the matrix is enumerated, not random.
set -u
cd "$(dirname "$0")/.."
CS=./crosswordsmith
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
fails=0
cases=0

# --- degenerate / edge-case input sets ---------------------------------------
printf '{"clues":[]}\n'                                              > "$WORK/empty.json"
printf 'clues([]).\n'                                               > "$WORK/empty.pl"
printf '{"clues":[{"answer":"A"},{"answer":"I"},{"answer":"AI"}]}\n' > "$WORK/single_letters.json"
printf '{"clues":[{"answer":"CAF\xc3\x89"},{"answer":"NA\xc3\x8fVE"},{"answer":"\xc3\x89CLAIR"}]}\n' > "$WORK/unicode.json"
printf '{"clues":[{"answer":"CAT"},{"answer":"CAT"}]}\n'            > "$WORK/dup.json"
printf '{"clues":[{"answer":"WELL-BEING"},{"answer":"BELOW"}]}\n'   > "$WORK/hyphen.json"
printf '{"clues":[{"answer":"DON%sT"},{"answer":"NOTE"}]}\n' "'"    > "$WORK/punct.json"
printf '{"clues":[{"answer":"R2D2"},{"answer":"C3PO"}]}\n'          > "$WORK/digits.json"
printf '{"clues":[{"answer":"cat"},{"answer":"car"}]}\n'            > "$WORK/lower.json"
printf '{"clues":[{"answer":"TT "},{"answer":"AT"}]}\n'             > "$WORK/trailing.json"
printf '{"clues":[{"answer":"AB"},{"answer":"BA"}]}\n'             > "$WORK/two_short.json"
printf '{ this is not json '                                       > "$WORK/malformed.json"
# a valid layout to feed lint/export:
$CS arrange --strict --size 17 \
    --input fixtures/bundled_17_clues.pl > "$WORK/layout.json" 2>/dev/null

# det "desc" cmd... : run 3x, require identical stdout + exit, defined exit (0/1)
det() {
    local desc="$1"; shift
    cases=$((cases + 1))
    timeout 90 "$@" >"$WORK/o1" 2>"$WORK/e1"; local e1=$?
    timeout 90 "$@" >"$WORK/o2" 2>/dev/null;  local e2=$?
    timeout 90 "$@" >"$WORK/o3" 2>/dev/null;  local e3=$?
    if [ "$e1" = 124 ] || [ "$e2" = 124 ] || [ "$e3" = 124 ]; then
        echo "  HANG          | $desc"; fails=$((fails + 1)); return
    fi
    if [ "$e1" != "$e2" ] || [ "$e1" != "$e3" ]; then
        echo "  NONDET-EXIT   | $desc (exits $e1/$e2/$e3)"; fails=$((fails + 1)); return
    fi
    if ! cmp -s "$WORK/o1" "$WORK/o2" || ! cmp -s "$WORK/o1" "$WORK/o3"; then
        echo "  NONDET-OUTPUT | $desc"; fails=$((fails + 1)); return
    fi
    if [ "$e1" != 0 ] && [ "$e1" != 1 ]; then
        echo "  BAD-EXIT($e1)  | $desc"; fails=$((fails + 1)); return
    fi
}

# --out behaviour. Two contracts:
#   nofile_on_fail : a genuine production failure (arrange unplaceable, fill
#                    unfillable) exits non-zero and MUST leave NO file - never a
#                    partial one (spec §5.2). On success the file is non-empty.
#   report_always  : lint always writes a COMPLETE report; its exit code is the
#                    VERDICT (AC-LINT-2), so a FAIL verdict legitimately pairs a
#                    complete (not partial) --out file with exit 1.
outcheck() {
    local mode="$1" desc="$2"; shift 2
    cases=$((cases + 1))
    local of="$WORK/out_$cases.json"
    rm -f "$of"
    timeout 90 "$@" --out "$of" >/dev/null 2>&1; local e=$?
    case "$mode" in
      nofile_on_fail)
        if [ "$e" = 0 ]; then
            [ -s "$of" ] || { echo "  OUT-EMPTY-OK0 | $desc (exit 0 but no file)"; fails=$((fails + 1)); }
        else
            [ -e "$of" ] && { echo "  PARTIAL-WRITE | $desc (exit $e left $of)"; fails=$((fails + 1)); }
        fi ;;
      report_always)
        [ -s "$of" ] || { echo "  NO-REPORT     | $desc (no complete report written)"; fails=$((fails + 1)); } ;;
    esac
}

echo "== determinism fuzz (3x byte-identity per case) =="

# --- arrange: flag matrix on a normal input ----------------------------------
N=fixtures/bundled_17_clues.pl
det "arrange strict fixed 17"        $CS arrange --strict --size 17 --input $N
det "arrange strict max 17"          $CS arrange --strict --max-size 17 --input $N
det "arrange strict default size"    $CS arrange --strict --input $N
det "arrange best-effort max 9"      $CS arrange --best-effort --max-size 9 --input $N
det "arrange best-effort fixed 11"   $CS arrange --best-effort --size 11 --input $N
det "arrange candidates 3"           $CS arrange --strict --size 17 --candidates 3 --input $N
det "arrange check-target 1"         $CS arrange --strict --size 17 --check-target 1 --input $N
det "arrange fragment"               $CS arrange --strict --fragment fixtures/bundled_17_fragment.json --input $N
det "arrange flag=value form"        $CS arrange --strict --size=17 --input=$N
det "arrange size+max-size (reject)" $CS arrange --strict --size 17 --max-size 20 --input $N
det "arrange enumerate small"        $CS arrange --enumerate --size 15 --input "$WORK/hyphen.json"

# --- arrange: degenerate inputs (graceful + deterministic) -------------------
for f in empty single_letters unicode dup hyphen punct digits lower trailing two_short malformed; do
    det "arrange degen $f (max)" $CS arrange --strict --max-size 15 --input "$WORK/$f.json"
    det "arrange degen $f (be)"  $CS arrange --best-effort --max-size 15 --input "$WORK/$f.json"
done
det "arrange empty .pl"          $CS arrange --strict --size 15 --input "$WORK/empty.pl"

# --- lint: every profile x allow-asymmetry x {valid, failing} ----------------
for p in toc blocked-uk american barred-ximenean; do
    det "lint $p valid"     $CS lint --profile $p "$WORK/layout.json"
    det "lint $p allow-asym" $CS lint --profile $p --allow-asymmetry "$WORK/layout.json"
    det "lint $p failing"   $CS lint --profile $p fixtures/lint_fail_layout.json
done

# --- export ------------------------------------------------------------------
det "export ipuz"   $CS export --to ipuz "$WORK/layout.json"
det "export exolve" $CS export --to exolve "$WORK/layout.json"

# --- fill --------------------------------------------------------------------
det "fill 3x3"        $CS fill --grid fixtures/fill_grid_3.json --dict fixtures/wordlist_sample.txt
det "fill 3x3 seed"   $CS fill --grid fixtures/fill_grid_3.json --seeds fixtures/fill_seed_3.json --dict fixtures/wordlist_sample.txt
det "fill 13a (infeasible w/ sample dict)" $CS fill --grid grids/blocked_13a.json --dict fixtures/wordlist_sample.txt
# scored fill (§8.4a): the scored ingestion + prune + score-desc ordering paths
det "fill scored (default prune)" $CS fill --grid fixtures/fill_grid_3.json --dict fixtures/dict_scored_sample.txt
det "fill scored --min-score 50"  $CS fill --grid fixtures/fill_grid_3.json --dict fixtures/dict_scored_sample.txt --min-score 50
det "fill scored min-score prunes all (infeasible)" $CS fill --grid fixtures/fill_grid_3.json --dict fixtures/dict_scored_sample.txt --min-score 99

# --- --out partial-write contract (§5.2) -------------------------------------
echo "== --out partial-write contract =="
outcheck nofile_on_fail "arrange ok --out"          $CS arrange --strict --size 17 --input $N
outcheck nofile_on_fail "arrange infeasible --out"  $CS arrange --strict --size 3 --input $N
outcheck report_always  "lint FAIL --out (complete report)" $CS lint --profile blocked-uk fixtures/lint_fail_layout.json
outcheck nofile_on_fail "fill infeasible --out"     $CS fill --grid grids/blocked_13a.json --dict fixtures/wordlist_sample.txt
outcheck nofile_on_fail "fill scored prune-all --out" $CS fill --grid fixtures/fill_grid_3.json --dict fixtures/dict_scored_sample.txt --min-score 99

echo "== summary: $cases cases, $fails failure(s) =="
[ "$fails" -eq 0 ] && { echo "DETERMINISM FUZZ PASSED"; exit 0; } || { echo "DETERMINISM FUZZ FAILED"; exit 1; }
