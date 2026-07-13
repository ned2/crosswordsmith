#!/usr/bin/env bash
# stamp-manifest.sh — artifact provenance for the wasm bundle (hardening plan
# Batch 2, docs/plans/wasm-supply-chain-batch2.md §1).
#
# For each of the four build outputs in $CLIENT_DIR (swipl-web.{js,wasm,data} +
# crosswordsmith.qlf) this records a sha256, copies it to a content-hashed
# sibling <stem>.<sha256:12><ext>, and writes build-manifest.json tying the set
# together with the swipl-devel commit, the compiled submodule SHAs
# (pins.sh $WASM_SUBMODULES), and the toolchain pins. worker.js prefers the hashed names when the manifest is
# present (falls back to unhashed when absent), so a CDN deploy ships the
# hashed set + worker.js + this manifest — a partial deploy can no longer pair
# a new js with a long-cached stale wasm.
#
# Deterministic on purpose: NO timestamps; buildId is derived from the four
# artifact hashes, so identical inputs stamp a byte-identical manifest.
# Provenance degrades explicitly: without a git tree at $SWIPL_SRC the swipl
# fields are null — visible, never fabricated.
#
# Standalone and build-free (~1s): called by build-wasm.sh step 5 and by
# wasm/test/run_all.sh after it re-stages/re-qcompiles.
#
# Usage:
#   wasm/build/stamp-manifest.sh                # stamp the in-repo wasm/client/
#   CLIENT_DIR=… SWIPL_SRC=… stamp-manifest.sh  # override locations
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
: "${CLIENT_DIR:=$REPO_ROOT/wasm/client}"
: "${SWIPL_SRC:=$HOME/src/swipl-devel}"

# shellcheck disable=SC1091
source "$SCRIPT_DIR/pins.sh"

ARTIFACTS="swipl-web.js swipl-web.wasm swipl-web.data crosswordsmith.qlf"

for a in $ARTIFACTS; do
  if [ ! -f "$CLIENT_DIR/$a" ]; then
    echo "stamp-manifest: missing $CLIENT_DIR/$a — build first (wasm/build/build-wasm.sh)" >&2
    exit 1
  fi
done

# profile provenance (payload plan Phases 2+3) — recorded ONLY when the
# caller states which profile produced the staged artifacts (build-wasm.sh
# step 2.6 and run_all.sh's profile-preferring staging both pass it). Null
# otherwise, matching the swipl-provenance policy: visible, never fabricated —
# this script cannot introspect a .data/.wasm to tell profile from full.
# packageList is the explicit -DSWIPL_PACKAGE_LIST configure input;
# staticExtensions is the ACTUAL foreign-plugin surface linked into the
# shipped swipl-web (everything else was dropped from the link); preload pins
# the .data keep-list. All from pins.sh — the single source of truth.
: "${PRELOAD_PROFILE:=}"
profile_json="null"
if [ -n "$PRELOAD_PROFILE" ]; then
  profile_file="$SCRIPT_DIR/preload-profile.txt"
  if [ ! -f "$profile_file" ]; then
    echo "stamp-manifest: PRELOAD_PROFILE=$PRELOAD_PROFILE but $profile_file is missing" >&2
    exit 1
  fi
  profile_sha="$(sha256sum "$profile_file" | cut -d' ' -f1)"
  profile_entries="$(grep -cvE '^\s*(#|$)' "$profile_file")"
  static_exts_json="$(printf '"%s", ' $WASM_STATIC_EXTENSIONS)"
  static_exts_json="[ ${static_exts_json%, } ]"
  profile_json="$(printf '{ "name": "%s", "packageList": "%s", "staticExtensions": %s, "preload": { "manifestSha256": "%s", "files": %s } }' \
    "$PRELOAD_PROFILE" "$WASM_PACKAGE_LIST" "$static_exts_json" "$profile_sha" "$profile_entries")"
fi

# swipl provenance — null (not guessed) when $SWIPL_SRC is not a git tree.
swipl_commit="null"
submodules_json="null"
if git -C "$SWIPL_SRC" rev-parse HEAD >/dev/null 2>&1; then
  swipl_commit="\"$(git -C "$SWIPL_SRC" rev-parse HEAD)\""
  # the WASM-compiled submodules (pins.sh $WASM_SUBMODULES — the same set
  # verify-pin.sh asserts on)
  # shellcheck disable=SC2086
  submodules_json="$(git -C "$SWIPL_SRC" submodule status $WASM_SUBMODULES 2>/dev/null \
    | awk '{ sha=$1; sub(/^[+\-U]/, "", sha); printf "%s\"%s\": \"%s\"", sep, $2, sha; sep=", " } END { if (NR==0) exit 1 }' \
    | sed 's/^/{ /; s/$/ }/')" || submodules_json="null"
fi

# licence notice: the deployable bundle must be self-contained — copy the
# notice in next to the artifacts so "deploy wasm/client/" ships it, and the
# manifest's "licenses" field can name a SIBLING file a consumer actually has
# (a repo-relative path is meaningless once the bundle leaves the repo).
cp "$REPO_ROOT/wasm/THIRD_PARTY_NOTICES.md" "$CLIENT_DIR/THIRD_PARTY_NOTICES.md"

# hash + content-named copies (clean stale hashed copies first — the glob
# cannot match the unhashed originals: they have no infix between stem and ext)
rm -f "$CLIENT_DIR"/swipl-web.*.js "$CLIENT_DIR"/swipl-web.*.wasm \
      "$CLIENT_DIR"/swipl-web.*.data "$CLIENT_DIR"/crosswordsmith.*.qlf
entry_lines=()
hashcat=""
for a in $ARTIFACTS; do
  sha="$(sha256sum "$CLIENT_DIR/$a" | cut -d' ' -f1)"
  bytes="$(wc -c < "$CLIENT_DIR/$a" | tr -d ' ')"
  stem="${a%.*}"; ext="${a##*.}"
  hashed="$stem.${sha:0:12}.$ext"
  cp "$CLIENT_DIR/$a" "$CLIENT_DIR/$hashed"
  entry_lines+=("$(printf '    "%s": { "hashed": "%s", "sha256": "%s", "bytes": %s }' \
                   "$a" "$hashed" "$sha" "$bytes")")
  hashcat="$hashcat$sha"
done
entries="$(printf '%s,\n' "${entry_lines[@]}")"
entries="${entries%,}"          # drop the trailing comma (keep valid JSON)
build_id="$(printf '%s' "$hashcat" | sha256sum | cut -d' ' -f1)"

cat > "$CLIENT_DIR/build-manifest.json" <<EOF
{
  "v": 1,
  "buildId": "$build_id",
  "profile": $profile_json,
  "swipl": { "commit": $swipl_commit, "submodules": $submodules_json },
  "toolchain": {
    "emsdk": "$EMSDK_VERSION",
    "zlib": "$ZLIB_VERSION",
    "zlibSha256": "$ZLIB_SHA256"
  },
  "licenses": "THIRD_PARTY_NOTICES.md",
  "artifacts": {
$entries
  }
}
EOF

# fail loudly on malformed JSON rather than shipping a manifest the worker
# would silently fall back from
python3 -c "import json,sys; json.load(open(sys.argv[1]))" "$CLIENT_DIR/build-manifest.json"
echo "stamped $CLIENT_DIR/build-manifest.json (buildId ${build_id:0:12}…)"
