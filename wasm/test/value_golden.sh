#!/usr/bin/env bash
# value_golden.sh - the wasm value-correctness golden + same-instance
# determinism, under the REAL wasm VM (node build.wasm/src/swipl.js). Native
# only otherwise: no browser, no server. Supersedes the spike's
# golden_parity.sh (byte-parity died with the envelope — DEC-8: `result` is a
# live dict, compared by value after parsing).
#
# For each case the SAME clue fixture feeds two paths:
#   reference : the real `crosswordsmith` CLI (its load_clues + flag resolvers)
#   browser   : a v1 request envelope through browser_dispatch/3 under wasm
# and the envelope's `result` must deep-equal the CLI's JSON. Cases cover
# seedless, the max crop, strict and best-effort. All bundled-17/seed/seedless
# requests run through ONE wasm Prolog instance, interleaved, so the
# per-request state reset is proven where it matters (a reused Worker
# instance) — the same-instance determinism debt from the deployment plan §10.4.
#
# SEEDED requests are the one deliberate exception to cross-VM equality
# (finding 2026-07-06, strategy §10): set_random(seed(N)) drives GMP's RNG
# natively but SWI's builtin RNG under the USE_GMP=OFF wasm build — same seed,
# different permutation, so native CLI seed:42 and wasm seed:42 CANNOT agree
# by design. Seeded reproducibility is engine-build-scoped. What we lock
# instead: (a) NATIVE dispatch of the seeded request == the CLI (the seam
# honours the seed identically within one VM); under wasm: (b) seed:42 twice
# is identical, (c) the seed genuinely perturbs the layout (≠ seedless — an
# ignored seed can't fake this), (d) provenance in diagnostics.
#
# Prereq: the swipl-devel wasm build tree (wasm/README.md). Override via
#   WASM_SWIPL=/path/to/build.wasm/src/swipl.js
# Run: wasm/test/value_golden.sh   (exit 0 = all checks pass)
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
cd "$ROOT"

WASM_SWIPL="${WASM_SWIPL:-$HOME/src/swipl-devel/build.wasm/src/swipl.js}"
if [ ! -f "$WASM_SWIPL" ]; then
  echo "value-golden FAILED: wasm swipl not found at $WASM_SWIPL" >&2
  echo "  build it per wasm/README.md, or set WASM_SWIPL" >&2
  exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

printf '{"clues":[{"answer":"CAT"},{"answer":"CAR"},{"answer":"ARC"},{"answer":"RAT"},{"answer":"TAR"}]}\n' \
  > "$TMP/toy_clues.json"

emit() { swipl -q wasm/test/value_golden.pl -- emit-request "$@"; }

echo "== building requests (native) =="
emit "$TMP/req_b17.json"        fixtures/bundled_17_clues.pl '{"size":17}'
emit "$TMP/req_b17_seed.json"   fixtures/bundled_17_clues.pl '{"size":17,"seed":42}'
emit "$TMP/req_toc_max.json"    fixtures/toc_demo.pl         '{"size":25,"mode":"max"}'
emit "$TMP/req_toy_strict.json" "$TMP/toy_clues.json"        '{"size":5}'
emit "$TMP/req_toy_be.json"     "$TMP/toy_clues.json"        '{"size":5,"bestEffort":true}'

echo "== building CLI references (native) =="
./crosswordsmith arrange --input fixtures/bundled_17_clues.pl --size 17            > "$TMP/ref_b17.json"      2>/dev/null
./crosswordsmith arrange --input fixtures/bundled_17_clues.pl --size 17 --seed 42  > "$TMP/ref_b17_seed.json" 2>/dev/null
./crosswordsmith arrange --input fixtures/toc_demo.pl --max-size 25                > "$TMP/ref_toc_max.json"  2>/dev/null
./crosswordsmith arrange --input "$TMP/toy_clues.json" --size 5                    > "$TMP/ref_toy_strict.json" 2>/dev/null
./crosswordsmith arrange --input "$TMP/toy_clues.json" --size 5 --best-effort      > "$TMP/ref_toy_be.json"   2>/dev/null

echo "== dispatching the seeded request natively (seed seam parity vs CLI) =="
swipl -q wasm/test/value_golden.pl -- dispatch "$TMP/req_b17_seed.json" \
  > "$TMP/native_seed.jsonl" 2>/dev/null

echo "== dispatching through ONE wasm instance =="
# Order is load-bearing: seedless / seeded / seedless-again / seeded-again on
# one instance is the determinism lock (a leaked seed would flip line 3 or 4).
node "$WASM_SWIPL" -q wasm/test/value_golden.pl -- dispatch \
  "$TMP/req_b17.json" "$TMP/req_b17_seed.json" \
  "$TMP/req_b17.json" "$TMP/req_b17_seed.json" \
  "$TMP/req_toy_strict.json" "$TMP/req_toy_be.json" "$TMP/req_toc_max.json" \
  > "$TMP/envelopes.jsonl" 2>"$TMP/wasm.err" || {
    echo "value-golden FAILED: wasm dispatch errored" >&2
    cat "$TMP/wasm.err" >&2
    exit 1
  }

echo "== comparing (deep value-equality after parse) =="
if python3 - "$TMP" <<'PY'
import json, sys
tmp = sys.argv[1]

envs = [json.loads(l) for l in open(f"{tmp}/envelopes.jsonl") if l.strip()]
assert len(envs) == 7, f"expected 7 envelopes, got {len(envs)}"

def ref(name):
    return json.load(open(f"{tmp}/ref_{name}.json"))

fails = []
def check(cond, name):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond: fails.append(name)

for e in envs:
    check(e.get("status") == "success" and e.get("v") == 1 and e.get("verb") == "arrange",
          f"envelope well-formed ({e.get('status')})")

b17, b17_seed, b17_again, b17_seed_again, toy_strict, toy_be, toc_max = envs

check(b17["result"] == ref("b17"),               "wasm result == CLI (bundled-17, seedless)")
check(toc_max["result"] == ref("toc_max"),       "wasm result == CLI (toc-demo, max crop)")
check(toy_strict["result"] == ref("toy_strict"), "wasm result == CLI (toy, strict)")
check(toy_be["result"] == ref("toy_be"),         "wasm result == CLI (toy, best-effort)")

# Seeded: cross-VM equality is impossible (GMP vs builtin RNG — see header);
# lock the seam natively and the semantics under wasm.
native_seed = json.loads(open(f"{tmp}/native_seed.jsonl").read())
check(native_seed["result"] == ref("b17_seed"),
      "NATIVE dispatch == CLI (bundled-17, seed 42): seed seam parity")
check(b17_seed["result"] != b17["result"],
      "wasm: seed 42 perturbs the layout (seed not silently ignored)")
check(b17_again["result"] == b17["result"],
      "same-instance determinism: seedless after seeded == seedless")
check(b17_seed_again["result"] == b17_seed["result"],
      "same-instance determinism: seed 42 twice identical")
check(b17_seed["result"]["diagnostics"]["arrange"].get("seed") == 42,
      "seed provenance in diagnostics")
check("seed" not in b17["result"]["diagnostics"]["arrange"],
      "no seed provenance on the seedless result")

sys.exit(1 if fails else 0)
PY
then
  echo "value-golden OK"
else
  echo "value-golden FAILED"
  exit 1
fi
