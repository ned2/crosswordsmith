#!/usr/bin/env bash
#
# benchmarks/check_fill_identity_artifact.sh - ARTIFACT-MODE byte-identity oracle
# for the fill ladder (F-L2 companion to check_fill_identity.sh).
#
# For every rung in benchmarks/fill_workloads.pl, run the product CLI in
# ARTIFACT-CONSUMING mode:
#
#     ./crosswordsmith fill --save-index A --dict D      # build once per dict
#     ./crosswordsmith fill --grid G --index A [--seeds S]
#
# and sha256 its STDOUT. Each digest is compared to the SAME committed manifest
# the raw-mode oracle pins (benchmarks/fill_identity.sha256). This is the F-L2
# gate: a precomputed index artifact must produce output byte-identical to the
# raw text-dict path on every rung, or the artifact is not equivalent and the
# change must be escalated (never silently accepted). Artifacts are built into a
# scratch dir and deleted on exit (nothing binary is committed - the browser
# milestone builds its own).
#
#   benchmarks/check_fill_identity_artifact.sh          # CHECK: diff, exit 0/1
#
# </dev/null so a load error can never park swipl at an stdin-blocking toplevel.

set -u
cd "$(dirname "$0")/.." || exit 2          # repo (worktree) root

MANIFEST="benchmarks/fill_identity.sha256"
WORKLOADS="benchmarks/fill_workloads.pl"
CLI="./crosswordsmith"

rungs="$(swipl -q -l "$WORKLOADS" \
  -g 'forall(fill_workload(Id,G,D,S,_,_,_,_,_), format("~w\t~w\t~w\t~w~n",[Id,G,D,S])), halt' \
  -g 'halt(1)' </dev/null 2>/dev/null)"
if [ -z "$rungs" ]; then
    echo "check_fill_identity_artifact: could not enumerate rungs from $WORKLOADS" >&2
    exit 2
fi

status=0
artdir="$(mktemp -d)"
seen="$(mktemp)"
trap 'rm -rf "$artdir" "$seen"' EXIT

# Build one artifact per distinct dict (cache keyed by a slugged dict path).
artifact_for() {
    local dict="$1"
    local slug art
    slug="$(printf '%s' "$dict" | tr -c 'A-Za-z0-9' '_')"
    art="$artdir/$slug.idx"
    if [ ! -f "$art" ]; then
        if ! "$CLI" fill --dict "$dict" --save-index "$art" </dev/null 2>/dev/null; then
            echo "check_fill_identity_artifact: FAILED to build artifact for $dict" >&2
            return 1
        fi
    fi
    printf '%s' "$art"
}

while IFS=$'\t' read -r id grid dict seeds; do
    [ -z "$id" ] && continue
    if ! art="$(artifact_for "$dict")"; then status=1; continue; fi
    args=(fill --grid "$grid" --index "$art")
    [ "$seeds" != none ] && args+=(--seeds "$seeds")
    out="$(mktemp)"
    if "$CLI" "${args[@]}" >"$out" 2>/dev/null; then
        digest="$(sha256sum "$out" | cut -d' ' -f1)"
        expected="$(awk -F'\t' -v r="$id" '$1==r{print $2}' "$MANIFEST" 2>/dev/null)"
        printf '%s\n' "$id" >>"$seen"
        if [ -z "$expected" ]; then
            echo "identity-artifact ($id): NO BASELINE (rung absent from $MANIFEST)"; status=1
        elif [ "$expected" = "$digest" ]; then
            echo "identity-artifact ($id): OK"
        else
            echo "identity-artifact ($id): MISMATCH (want $expected got $digest)"; status=1
        fi
    else
        echo "identity-artifact ($id): CLI FAILED (fill --index did not exit 0)"; status=1
    fi
    rm -f "$out"
done <<<"$rungs"

# A manifest rung this pass never measured is a silent drop - fail it.
if [ -f "$MANIFEST" ]; then
    while IFS=$'\t' read -r mid _; do
        [ -z "$mid" ] && continue
        if ! grep -qxF "$mid" "$seen"; then
            echo "identity-artifact ($mid): IN MANIFEST BUT NOT MEASURED (rung dropped?)"; status=1
        fi
    done <"$MANIFEST"
fi

if [ "$status" -eq 0 ]; then
    echo "IDENTITY (artifact mode): all rung digests match $MANIFEST"
else
    echo "IDENTITY (artifact mode): FAILED"
fi
exit "$status"
