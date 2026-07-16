#!/usr/bin/env bash
# Semantic identity gate for the complete greedy ladder. The generated JSON
# includes ordered direct attempts, exact mode-specific raw pools, selected
# semantics, normalized assocs, distances, and CLI stdout/stderr hashes; the
# committed file stores its SHA-256.

set -euo pipefail
cd "$(dirname "$0")/.."

manifest="benchmarks/greedy_identity.sha256"
actual="$(mktemp)"
trap 'rm -f "$actual"' EXIT

swipl -q benchmarks/run_arrange_greedy.pl -- --identity --heavy >"$actual"
digest="$(sha256sum "$actual" | cut -d' ' -f1)"

if [[ "${1:-}" == "--record" ]]; then
    printf '%s  greedy-identity.json\n' "$digest" >"$manifest"
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
