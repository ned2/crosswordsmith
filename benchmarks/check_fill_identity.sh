#!/usr/bin/env bash
#
# benchmarks/check_fill_identity.sh - CLI byte-identity oracle for the fill ladder.
#
# For every rung in benchmarks/fill_workloads.pl, run the product CLI
#
#     ./crosswordsmith fill --grid G --dict D [--seeds S]
#
# and sha256 its STDOUT (the canonical filled layout; stderr discarded). Each
# digest is compared to the committed manifest benchmarks/fill_identity.sha256.
#
# This is Phase 2's EQUIVALENCE GATE: fill's output is byte-identical and
# deterministic, so any search/index optimisation must leave every rung's stdout
# digest unchanged. A mismatch means the change altered the fill (order, choice,
# or content) and must be escalated to the orchestrator, never silently accepted.
# It is the CLI-level companion to the inference ratchet (check_fill_baseline.pl):
# the ratchet proves the work didn't grow, this proves the answer didn't move.
#
#   benchmarks/check_fill_identity.sh            # CHECK: diff digests, exit 0/1
#   benchmarks/check_fill_identity.sh --record   # (re)generate the manifest
#
# The manifest is COMMITTED; regenerating it is an explicit, reviewed act (the
# same discipline as the goldens). The rung list is read from fill_workloads.pl
# (single source of truth) so this never drifts from the ladder.
#
# NB every ladder rung completes as `filled` under the CLI's shipped 800M budget
# (all rungs top out at ~34.9M search inferences); a non-zero CLI exit here is a
# real failure, reported per rung.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 2          # repo (worktree) root

# Overrides let the safety test exercise recording without touching committed data.
MANIFEST="${FILL_IDENTITY_MANIFEST:-benchmarks/fill_identity.sha256}"
WORKLOADS="${FILL_IDENTITY_WORKLOADS:-benchmarks/fill_workloads.pl}"
CLI="${FILL_IDENTITY_CLI:-./crosswordsmith}"

RECORD=0
for arg in "$@"; do
    case "$arg" in
        --record) RECORD=1 ;;
        *) echo "check_fill_identity: unknown option $arg" >&2; exit 2 ;;
    esac
done

# Enumerate rungs as "rungid<TAB>grid<TAB>dict<TAB>seeds" (seeds = the atom `none`
# or a fragment path) straight from the manifest of record. </dev/null so a load
# error can never drop swipl into an stdin-blocking toplevel.
if ! rungs="$(swipl -q -l "$WORKLOADS" \
  -g 'forall(fill_workload(Id,G,D,S,_,_,_,_,_), format("~w\t~w\t~w\t~w~n",[Id,G,D,S])), halt' \
  -g 'halt(1)' </dev/null 2>/dev/null)"; then
    echo "check_fill_identity: failed to enumerate rungs from $WORKLOADS" >&2
    exit 2
fi
if [ -z "$rungs" ]; then
    echo "check_fill_identity: could not enumerate rungs from $WORKLOADS" >&2
    exit 2
fi

declare -A expected
if [ "$RECORD" -eq 0 ] && [ -f "$MANIFEST" ]; then
    while IFS=$'\t' read -r id digest extra || [ -n "$id$digest$extra" ]; do
        [ -z "$id" ] && continue
        if [ -n "${expected[$id]+present}" ]; then
            echo "check_fill_identity: duplicate manifest id $id" >&2
            exit 2
        fi
        if [[ ! "$digest" =~ ^[0-9a-f]{64}$ ]] || [ -n "$extra" ]; then
            echo "check_fill_identity: malformed manifest entry $id" >&2
            exit 2
        fi
        expected["$id"]="$digest"
    done <"$MANIFEST"
fi

status=0
if [ "$RECORD" -eq 1 ]; then
    if ! tmp_manifest="$(mktemp "${MANIFEST}.tmp.XXXXXX")"; then
        echo "check_fill_identity: could not allocate manifest temporary" >&2
        exit 2
    fi
else
    if ! tmp_manifest="$(mktemp)"; then
        echo "check_fill_identity: could not allocate digest temporary" >&2
        exit 2
    fi
fi
if ! seen="$(mktemp)"; then
    echo "check_fill_identity: could not allocate seen-set temporary" >&2
    rm -f "$tmp_manifest"
    exit 2
fi
if ! out="$(mktemp)"; then
    echo "check_fill_identity: could not allocate output temporary" >&2
    rm -f "$tmp_manifest" "$seen"
    exit 2
fi
trap 'rm -f "$tmp_manifest" "$seen" "$out"' EXIT
declare -A known
known_count=0
seen_count=0

while IFS=$'\t' read -r id grid dict seeds; do
    [ -z "$id" ] && continue
    if [ -n "${known[$id]:-}" ]; then
        echo "identity ($id): DUPLICATE WORKLOAD ID"; status=1
        continue
    fi
    known["$id"]=1
    known_count=$((known_count + 1))
    args=(fill --grid "$grid" --dict "$dict")
    [ "$seeds" != none ] && args+=(--seeds "$seeds")
    if "$CLI" "${args[@]}" >"$out" 2>/dev/null; then
        hash_output=''
        if ! hash_output="$(sha256sum "$out")" \
           || ! read -r digest _ <<<"$hash_output" \
           || [[ ! "$digest" =~ ^[0-9a-f]{64}$ ]]; then
            echo "identity ($id): SHA-256 FAILED"; status=1
            continue
        fi
        if ! printf '%s\t%s\n' "$id" "$digest" >>"$tmp_manifest"; then
            echo "identity ($id): MANIFEST TEMPORARY WRITE FAILED"; status=1
            continue
        fi
        seen_count=$((seen_count + 1))
        if [ "$RECORD" -eq 0 ]; then
            if ! printf '%s\n' "$id" >>"$seen"; then
                echo "identity ($id): SEEN-SET WRITE FAILED"; status=1
                continue
            fi
            if [ -z "${expected[$id]+present}" ]; then
                echo "identity ($id): NO BASELINE (rung absent from $MANIFEST)"; status=1
            elif [ "${expected[$id]}" = "$digest" ]; then
                echo "identity ($id): OK"
            else
                echo "identity ($id): MISMATCH (want ${expected[$id]} got $digest)"; status=1
            fi
        fi
    else
        echo "identity ($id): CLI FAILED (fill did not exit 0)"; status=1
    fi
done <<<"$rungs"

if [ "$RECORD" -eq 1 ]; then
    if [ "$status" -eq 0 ] && [ "$seen_count" -eq "$known_count" ] \
       && sort -o "$tmp_manifest" "$tmp_manifest" \
       && mv "$tmp_manifest" "$MANIFEST"; then
        echo "recorded $seen_count complete rung digest(s) -> $MANIFEST"
        exit 0
    fi
    echo "IDENTITY: record FAILED; $MANIFEST was not changed"
    exit 1
fi

# A manifest rung that this pass never measured is a silent drop - fail it (the
# same discipline as the baseline recorder's new-rung guard).
for mid in "${!expected[@]}"; do
    if ! grep -qxF "$mid" "$seen"; then
        echo "identity ($mid): IN MANIFEST BUT NOT MEASURED (rung dropped?)"; status=1
    fi
done

if [ "$status" -eq 0 ]; then
    echo "IDENTITY: all rung digests match $MANIFEST"
else
    echo "IDENTITY: FAILED"
fi
exit "$status"
