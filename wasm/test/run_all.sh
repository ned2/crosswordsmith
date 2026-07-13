#!/usr/bin/env bash
# run_all.sh - the one-command wasm test battery (`make test-wasm`).
#
# From a cold checkout with the swipl-devel wasm build tree present, this
# stages everything and runs every wasm-side gate that used to be hand-run
# (deployment plan §10.4 "CI entry point"):
#
#   0. preflight   the build tree (node swipl.js + swipl-web.{js,wasm,data})
#                  and wasm/test/node_modules must exist - fail loudly with
#                  the fix, never skip silently
#   1. stage       copy swipl-web.{js,wasm,data} into wasm/client/ (skipped
#                  when already fresh)
#   2. qcompile    wasm/client/solve_browser.pl -> crosswordsmith.qlf with the
#                  WASM node swipl - NEVER native swipl: qlf is word-size-
#                  specific, a native x86-64 qlf would load-and-misbehave
#                  under wasm32 (wasm/README.md "Why the app .qlf must be
#                  wasm-produced"). Skipped when no input is newer.
#   2.5 stamp      build-manifest.json + content-hashed artifact names
#                  (wasm/build/stamp-manifest.sh) so the battery drives the
#                  manifest-preferring load path in worker.js
#   3. battery     value_golden.sh (21 checks), golden_type_check.mjs (4
#                  goldens), then a self-managed static server on a free port
#                  for headless.mjs (16 scenarios) + error_probe.mjs (4 cases).
#                  yield_probe.mjs (gate #1, ~60s) is opt-in via --yield.
#                  Every step runs even after a failure; exit is non-zero if
#                  ANY step failed.
#
# The server's web root is wasm/ - NOT wasm/client/: sdk/ and client/ are
# siblings and the harness page imports ../sdk/crosswordsmith.mjs.
#
# Env overrides:
#   WASM_BUILD=/path/to/build.wasm     the swipl-devel wasm build tree
#   WASM_SWIPL=/path/to/swipl.js       the wasm node swipl (default: under
#                                      $WASM_BUILD/src; setting only this
#                                      derives WASM_BUILD from it)
# Run: wasm/test/run_all.sh [--yield]   (exit 0 = all steps pass)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
CLIENT="$ROOT/wasm/client"
cd "$ROOT"

RUN_YIELD=0
for arg in "$@"; do
  case "$arg" in
    --yield) RUN_YIELD=1 ;;
    *) echo "usage: $0 [--yield]" >&2; exit 2 ;;
  esac
done

# --- 0. preflight -------------------------------------------------------------
if [ -n "${WASM_SWIPL:-}" ] && [ -z "${WASM_BUILD:-}" ]; then
  WASM_BUILD="$(dirname "$(dirname "$WASM_SWIPL")")"
fi
WASM_BUILD="${WASM_BUILD:-$HOME/src/swipl-devel/build.wasm}"
WASM_SWIPL="${WASM_SWIPL:-$WASM_BUILD/src/swipl.js}"
export WASM_SWIPL   # value_golden.sh dispatches through the same wasm swipl

for f in "$WASM_SWIPL" "$WASM_BUILD/src/swipl-web.js" \
         "$WASM_BUILD/src/swipl-web.wasm" "$WASM_BUILD/src/swipl-web.data"; do
  if [ ! -f "$f" ]; then
    echo "test-wasm FAILED: missing $f" >&2
    echo "  The swipl-devel wasm build tree is required. Build it per wasm/README.md:" >&2
    echo "    wasm/build/build-wasm.sh" >&2
    echo "  or point WASM_BUILD (or WASM_SWIPL) at an existing build.wasm tree." >&2
    exit 1
  fi
done

if [ ! -d "$HERE/node_modules" ]; then
  echo "test-wasm FAILED: wasm/test/node_modules missing (typescript + playwright)." >&2
  echo "  Install once from the committed lockfile - the tests drive system Chrome" >&2
  echo "  (Playwright channel:'chrome'), so skip the browser download:" >&2
  echo "    ( cd wasm/test && PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 npm ci )" >&2
  exit 1
fi

# --- 1. stage the engine artifacts (skip when fresh) ---------------------------
# Prefer the crosswordsmith-web profile image when the build tree carries one
# (build-wasm.sh step 2.6, payload plan Phase 2) — that is what production
# ships, so it is what the battery must gate. The plain src/swipl-web.* set is
# SWI's full default image (kept for diagnostics/hand-staged trees).
ENGINE_SRC="$WASM_BUILD/src"
PRELOAD_PROFILE=""
if [ -f "$WASM_BUILD/src/profile-crosswordsmith-web/swipl-web.js" ]; then
  ENGINE_SRC="$WASM_BUILD/src/profile-crosswordsmith-web"
  PRELOAD_PROFILE="crosswordsmith-web"
fi
echo "== staging (build tree: $WASM_BUILD; engine: ${PRELOAD_PROFILE:-full image}) =="
for a in swipl-web.js swipl-web.wasm swipl-web.data; do
  if [ ! -f "$CLIENT/$a" ] || [ "$ENGINE_SRC/$a" -nt "$CLIENT/$a" ]; then
    cp "$ENGINE_SRC/$a" "$CLIENT/$a"
    echo "  staged $a"
  else
    echo "  fresh  $a - skip"
  fi
done

# --- 2. qcompile the app qlf with the WASM swipl (skip when fresh) -------------
# Inputs: the load root + every project module include(user) folds in + the
# staged engine (a new swipl-web.js means a new VM signature - requalify).
QLF="$CLIENT/crosswordsmith.qlf"
stale=""
if [ ! -f "$QLF" ]; then
  stale="no qlf yet"
elif [ -n "$(find "$CLIENT/solve_browser.pl" "$ROOT/load.pl" "$ROOT/prolog" \
                  -name '*.pl' -newer "$QLF" -print -quit)" ]; then
  stale="a .pl source is newer"
elif [ "$CLIENT/swipl-web.js" -nt "$QLF" ]; then
  stale="the engine is newer"
fi
if [ -n "$stale" ]; then
  echo "  qcompile -> crosswordsmith.qlf ($stale; wasm swipl - qlf is word-size-specific)"
  # qcompile writes the .qlf ALONGSIDE ITS SOURCE (basename from the source,
  # independent of CWD) - rename it in place, as build-wasm.sh step 4 does.
  node "$WASM_SWIPL" \
    -g "qcompile('$CLIENT/solve_browser.pl', [include(user)]), halt" \
    -t 'halt(1)'
  mv "$CLIENT/solve_browser.qlf" "$QLF"
else
  echo "  fresh  crosswordsmith.qlf - skip qcompile"
fi

# --- 2.1 qlf dependency inventory (payload plan Phase 1) ------------------------
# The browser load root (solve_browser.pl) imports exactly the browser surface:
# core/metrics/arrange/lint/export/browser. fill and stockgrid must NOT be in
# the qlf's folded closure - they drag fill-only baggage (library(sha),
# library(fastrw), dictionary/stock-grid handling) into every browser cold
# load, and the reduced preload/package profiles (Phases 2-3) are derived from
# THIS closure. qlf files embed their source paths (provenance for error
# messages), so grep -a is a faithful inventory. The browser.pl positive probe
# keeps this guard honest: if qlf provenance strings ever vanish, the check
# fails loudly instead of silently passing.
if ! grep -aq 'crosswordsmith/browser\.pl' "$QLF"; then
  echo "test-wasm FAILED: qlf inventory probe broken - no browser.pl provenance in $QLF" >&2
  exit 1
fi
LEAKED="$(grep -aoE 'crosswordsmith/(fill|stockgrid)\.pl|library\(sha\)|fastrw' "$QLF" | sort -u || true)"
if [ -n "$LEAKED" ]; then
  echo "test-wasm FAILED: fill-only dependencies leaked into the browser qlf:" >&2
  echo "$LEAKED" | sed 's/^/    /' >&2
  echo "  the browser load root (wasm/client/solve_browser.pl) must import only" >&2
  echo "  core/metrics/arrange/lint/export/browser (payload plan Phase 1)." >&2
  exit 1
fi
echo "  qlf inventory OK (no fill/stockgrid/sha/fastrw in the closure)"

# --- 2.5 stamp provenance (build-manifest.json + content-hashed names) ---------
# Steps 1-2 may have replaced artifacts, so re-stamp; worker.js prefers the
# manifest's hashed names when present, so the battery below exercises the
# real manifest-driven load path (supply-chain Batch 2). Deterministic and ~1s.
CLIENT_DIR="$CLIENT" SWIPL_SRC="$(dirname "$WASM_BUILD")" \
  PRELOAD_PROFILE="$PRELOAD_PROFILE" \
  "$ROOT/wasm/build/stamp-manifest.sh"

# --- 3. the battery ------------------------------------------------------------
status=0
results=()
step() {
  local name="$1"; shift
  echo
  echo "=== $name ==="
  if "$@"; then
    results+=("PASS  $name")
  else
    results+=("FAIL  $name")
    status=1
  fi
}

step "value golden (wasm deep-equals CLI)" wasm/test/value_golden.sh
step "golden type lock (goldens satisfy the .d.ts)" node wasm/test/golden_type_check.mjs

# Static server for the browser probes, on an OS-assigned free port. -u: python
# block-buffers the "Serving HTTP on ... port N" line into a pipe, and we need
# to read N back.
SERVER_LOG="$(mktemp)"
SERVER_PID=""
cleanup() {
  if [ -n "$SERVER_PID" ]; then
    kill "$SERVER_PID" 2>/dev/null || true
    wait "$SERVER_PID" 2>/dev/null || true
  fi
  rm -f "$SERVER_LOG"
}
trap cleanup EXIT

( cd wasm && exec python3 -u -m http.server 0 --bind 127.0.0.1 ) >"$SERVER_LOG" 2>&1 &
SERVER_PID=$!
PORT=""
for _ in $(seq 1 50); do
  PORT="$(sed -n 's/.*port \([0-9][0-9]*\).*/\1/p' "$SERVER_LOG" | head -n1)"
  [ -n "$PORT" ] && break
  if ! kill -0 "$SERVER_PID" 2>/dev/null; then break; fi
  sleep 0.1
done
if [ -z "$PORT" ]; then
  echo "test-wasm FAILED: static server did not come up" >&2
  cat "$SERVER_LOG" >&2
  exit 1
fi
HARNESS="http://127.0.0.1:$PORT/client/harness.html"
echo
echo "== serving wasm/ at http://127.0.0.1:$PORT (pid $SERVER_PID) =="

step "headless SDK end-to-end" env URL="$HARNESS" node wasm/test/headless.mjs
step "worker error paths" env URL="$HARNESS?noauto=1" node wasm/test/error_probe.mjs
step "spare-worker policies" env URL="$HARNESS?noauto=1" node wasm/test/spare_policy_probe.mjs
if [ "$RUN_YIELD" = 1 ]; then
  step "gate #1 responsiveness + cancel (~60s)" env URL="$HARNESS?noauto=1" node wasm/test/yield_probe.mjs
else
  echo
  echo "(yield_probe.mjs skipped - opt in: make test-wasm WASM_TEST_ARGS=--yield)"
fi

echo
echo "=== summary ==="
printf '%s\n' "${results[@]}"
if [ "$status" -eq 0 ]; then
  echo "test-wasm OK"
else
  echo "test-wasm FAILED"
fi
exit "$status"
