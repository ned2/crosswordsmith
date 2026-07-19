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
# FILL_SAVE_INDEX_FLAGS (env, optional): extra flags for the artifact build.
# FILL_SAVE_INDEX_FLAGS=--masks proves identity with the F-H2 bitset counting
# masks ACTIVE (the bignum count kernel); unset/empty exercises the default
# no-masks v2 artifact (the ordset kernel). Both must match the same manifest.
#
# </dev/null so a load error can never park swipl at an stdin-blocking toplevel.

set -uo pipefail
cd "$(dirname "$0")/.." || exit 2          # repo (worktree) root

MANIFEST="${FILL_IDENTITY_MANIFEST:-benchmarks/fill_identity.sha256}"
WORKLOADS="${FILL_IDENTITY_WORKLOADS:-benchmarks/fill_workloads.pl}"
CLI="${FILL_IDENTITY_CLI:-./crosswordsmith}"

if [ "$#" -ne 0 ]; then
    echo "check_fill_identity_artifact: unknown option $1" >&2
    exit 2
fi

if ! rungs="$(swipl -q -l "$WORKLOADS" \
  -g 'forall(fill_workload(Id,G,D,S,_,_,_,_,_), format("~w\t~w\t~w\t~w~n",[Id,G,D,S])), halt' \
  -g 'halt(1)' </dev/null 2>/dev/null)"; then
    echo "check_fill_identity_artifact: failed to enumerate rungs from $WORKLOADS" >&2
    exit 2
fi
if [ -z "$rungs" ]; then
    echo "check_fill_identity_artifact: could not enumerate rungs from $WORKLOADS" >&2
    exit 2
fi

declare -A expected
if [ -f "$MANIFEST" ]; then
    while IFS=$'\t' read -r id digest extra || [ -n "$id$digest$extra" ]; do
        [ -z "$id" ] && continue
        if [ -n "${expected[$id]+present}" ]; then
            echo "check_fill_identity_artifact: duplicate manifest id $id" >&2
            exit 2
        fi
        if [[ ! "$digest" =~ ^[0-9a-f]{64}$ ]] || [ -n "$extra" ]; then
            echo "check_fill_identity_artifact: malformed manifest entry $id" >&2
            exit 2
        fi
        expected["$id"]="$digest"
    done <"$MANIFEST"
fi

status=0
if ! artdir="$(mktemp -d)"; then
    echo "check_fill_identity_artifact: could not allocate artifact directory" >&2
    exit 2
fi
if ! seen="$(mktemp)"; then
    echo "check_fill_identity_artifact: could not allocate seen-set temporary" >&2
    rm -rf "$artdir"
    exit 2
fi
trap 'rm -rf "$artdir" "$seen"' EXIT
out="$artdir/out"
declare -A known

# Build one artifact per distinct dict (cache keyed by a slugged dict path).
# FILL_SAVE_INDEX_FLAGS (word-split on purpose) selects the artifact flavour.
artifact_for() {
    local dict="$1"
    local slug art
    slug="$(printf '%s' "$dict" | tr -c 'A-Za-z0-9' '_')"
    art="$artdir/$slug.idx"
    if [ ! -f "$art" ]; then
        # shellcheck disable=SC2086
        if ! "$CLI" fill --dict "$dict" --save-index "$art" ${FILL_SAVE_INDEX_FLAGS:-} </dev/null 2>/dev/null; then
            echo "check_fill_identity_artifact: FAILED to build artifact for $dict" >&2
            return 1
        fi
    fi
    printf '%s' "$art"
}

while IFS=$'\t' read -r id grid dict seeds; do
    [ -z "$id" ] && continue
    if [ -n "${known[$id]:-}" ]; then
        echo "identity-artifact ($id): DUPLICATE WORKLOAD ID"; status=1
        continue
    fi
    known["$id"]=1
    if ! art="$(artifact_for "$dict")"; then status=1; continue; fi
    args=(fill --grid "$grid" --index "$art")
    [ "$seeds" != none ] && args+=(--seeds "$seeds")
    if "$CLI" "${args[@]}" >"$out" 2>/dev/null; then
        hash_output=''
        if ! hash_output="$(sha256sum "$out")" \
           || ! read -r digest _ <<<"$hash_output" \
           || [[ ! "$digest" =~ ^[0-9a-f]{64}$ ]]; then
            echo "identity-artifact ($id): SHA-256 FAILED"; status=1
            continue
        fi
        if ! printf '%s\n' "$id" >>"$seen"; then
            echo "identity-artifact ($id): SEEN-SET WRITE FAILED"; status=1
            continue
        fi
        if [ -z "${expected[$id]+present}" ]; then
            echo "identity-artifact ($id): NO BASELINE (rung absent from $MANIFEST)"; status=1
        elif [ "${expected[$id]}" = "$digest" ]; then
            echo "identity-artifact ($id): OK"
        else
            echo "identity-artifact ($id): MISMATCH (want ${expected[$id]} got $digest)"; status=1
        fi
    else
        echo "identity-artifact ($id): CLI FAILED (fill --index did not exit 0)"; status=1
    fi
done <<<"$rungs"

# A manifest rung this pass never measured is a silent drop - fail it.
for mid in "${!expected[@]}"; do
    if ! grep -qxF "$mid" "$seen"; then
        echo "identity-artifact ($mid): IN MANIFEST BUT NOT MEASURED (rung dropped?)"; status=1
    fi
done

if [ "$status" -eq 0 ]; then
    echo "IDENTITY (artifact mode): all rung digests match $MANIFEST"
else
    echo "IDENTITY (artifact mode): FAILED"
fi
exit "$status"
