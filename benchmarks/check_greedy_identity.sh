#!/usr/bin/env bash
# Semantic identity gate for the complete greedy ladder. The generated JSON
# includes ordered direct attempts, exact mode-specific raw pools, selected
# semantics, normalized assocs, distances, and CLI stdout/stderr hashes; the
# committed file stores its SHA-256.

set -euo pipefail
cd "$(dirname "$0")/.."

manifest="benchmarks/greedy_identity.sha256"
record=0
tmp_manifest=""

for arg in "$@"; do
    case "$arg" in
        --record) record=1 ;;
        *) echo "check_greedy_identity: unknown option $arg" >&2; exit 2 ;;
    esac
done

actual="$(mktemp)"
if [[ "$record" -eq 1 ]]; then
    tmp_manifest="$(mktemp "${manifest}.tmp.XXXXXX")"
fi
trap 'rm -f "$actual"; [[ -z "$tmp_manifest" ]] || rm -f "$tmp_manifest"' EXIT

swipl -q benchmarks/run_arrange_greedy.pl -- --identity --heavy >"$actual"
digest="$(sha256sum "$actual" | cut -d' ' -f1)"

if [[ "$record" -eq 1 ]]; then
    printf '%s  greedy-identity.json\n' "$digest" >"$tmp_manifest"
    mv "$tmp_manifest" "$manifest"
    printf 'recorded greedy identity: %s\n' "$digest"
    exit 0
fi

expected="$(cut -d' ' -f1 "$manifest" 2>/dev/null || true)"
if [[ -n "$expected" && "$digest" == "$expected" ]]; then
    printf 'GREEDY IDENTITY: OK (%s)\n' "$digest"
else
    printf 'GREEDY IDENTITY: FAILED (want %s got %s)\n' "${expected:-missing}" "$digest" >&2
    exit 1
fi
