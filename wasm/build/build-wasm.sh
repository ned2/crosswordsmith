#!/usr/bin/env bash
# build-wasm.sh — reproducible build of the crosswordsmith browser bundle.
#
# Turns the prose recipe in docs/plans/wasm-browser-deployment.md §2 into one
# runnable script. It (0) activates emsdk, (1) stages zlib into $WASM_HOME,
# (2) builds swipl-web from our pinned swipl-devel commit (cloning it first on
# a cold runner) with the crosswordsmith-web package profile (pins.sh
# $WASM_PACKAGE_LIST, payload plan Phase 3), (2.6) relinks the browser image
# against the minimal preload keep-list + static-extension list (Phases 2+3),
# (3) copies the three web artifacts into wasm/client/, (4) qcompiles the app
# into wasm/client/crosswordsmith.qlf using the *wasm* swipl (correct pointer
# size), and (5) stamps build-manifest.json + content-hashed artifact names
# (provenance — see wasm/build/stamp-manifest.sh).
#
# Everything it produces under wasm/client/ is gitignored (see .gitignore) — the
# repo tracks this recipe, not the ~1.8 MB of binaries.
#
# First run took a couple of minutes for configure (hundreds of cross-compile
# probes) + ~600 compile steps. Re-runs skip already-staged deps.
#
# Usage:
#   wasm/build/build-wasm.sh            # full build (verifies the swipl-devel pin)
#   SWIPL_ALLOW_CHECKOUT=1 …            # allow the script to `git checkout` the pin
#   SWIPL_ALLOW_DIRTY=1 …               # build despite a dirty swipl-devel tree
#   WASM_HOME=… SWIPL_SRC=… …          # override locations
#   SWIPL_REPO_URL=…                    # remote for the cold-runner clone
set -euo pipefail

# --- configuration (override via env) ---------------------------------------
: "${WASM_HOME:=$HOME/wasm}"                 # staging prefix for cross-compiled deps
: "${SWIPL_SRC:=$HOME/src/swipl-devel}"      # the swipl-devel checkout to build from

# repo root = two levels up from this script (wasm/build/ -> repo root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# every supply-chain pin (SWIPL_COMMIT, EMSDK/ZLIB versions + hashes) and the
# crosswordsmith-web package/static-extension profile live in pins.sh — one
# source of truth shared with stamp-manifest.sh.
# shellcheck disable=SC1091
source "$SCRIPT_DIR/pins.sh"
CLIENT_DIR="$REPO_ROOT/wasm/client"
BUILD_DIR="$SWIPL_SRC/build.wasm"

log() { printf '\n\033[1;36m== %s\033[0m\n' "$*"; }

# supply-chain guards (submodule/dirty-tree/emcc/pin asserts) live in a
# sourced helper so each is unit-testable in ~1s without the ~600-step ninja
# build. See docs/plans/wasm-supply-chain-hardening.md.
# shellcheck disable=SC1091
source "$SCRIPT_DIR/verify-pin.sh"

# --- 0. toolchain: emsdk -----------------------------------------------------
log "0. emsdk $EMSDK_VERSION"
if [ ! -d "$WASM_HOME/emsdk" ]; then
  # --branch pins the repo tag too (emsdk tags per SDK version) — cosmetic
  # defense-in-depth; verify_emcc_version below is the load-bearing guard.
  git clone --branch "$EMSDK_VERSION" --depth 1 \
      https://github.com/emscripten-core/emsdk "$WASM_HOME/emsdk"
fi
( cd "$WASM_HOME/emsdk" && ./emsdk install "$EMSDK_VERSION" && ./emsdk activate "$EMSDK_VERSION" )
# shellcheck disable=SC1091
source "$WASM_HOME/emsdk/emsdk_env.sh"
# assert the ACTIVE compiler is the pinned one (the load-bearing emsdk guard —
# activation resolves the toolchain from the version number, not the repo HEAD).
verify_emcc_version "$EMSDK_VERSION"

# --- 1. stage wasm-side deps into $WASM_HOME --------------------------------
# zlib only: the SWI core links it (saved states / qlf are zip-packaged).
# pcre2 is gone with the pcre package (payload plan Phase 3) — nothing in the
# crosswordsmith-web package profile links it.
log "1. deps: zlib $ZLIB_VERSION -> $WASM_HOME"
if [ ! -f "$WASM_HOME/lib/libz.a" ]; then
  # fossils/ is the DURABLE path (the bare /zlib-… path serves only the current
  # release and 404s once superseded); -f fails on 4xx/5xx instead of saving an
  # HTML error body; download → sha256sum -c → extract (can't checksum a stream
  # already piped into tar); rm -rf first so we never extract over a dirty dir.
  ( cd "$WASM_HOME"
    rm -rf "zlib-$ZLIB_VERSION"
    curl -fsSL "https://zlib.net/fossils/zlib-$ZLIB_VERSION.tar.gz" \
      -o "zlib-$ZLIB_VERSION.tar.gz"
    echo "$ZLIB_SHA256  zlib-$ZLIB_VERSION.tar.gz" | sha256sum -c -
    tar xzf "zlib-$ZLIB_VERSION.tar.gz"
    rm -f "zlib-$ZLIB_VERSION.tar.gz"
    cd "zlib-$ZLIB_VERSION" && emconfigure ./configure --static --prefix="$WASM_HOME"
    EMCC_CFLAGS=-Wno-deprecated-non-prototype emmake make && emmake make install )
else
  echo "  libz.a present — skip"
fi
# --- 2. build swipl-web from OUR pinned commit ------------------------------
log "2. swipl-web from $SWIPL_SRC @ $SWIPL_COMMIT (packages: $WASM_PACKAGE_LIST)"
# Cold runner / standalone CI (Batch 2): when $SWIPL_SRC is absent, clone the
# pin + init the WASM submodules (never mutates a pre-existing tree).
bootstrap_swipl_src "$SWIPL_SRC" "$SWIPL_COMMIT"
# Supply-chain guards (all read-only unless SWIPL_ALLOW_CHECKOUT=1), each in
# verify-pin.sh so they unit-test in ~1s without the ninja build:
#  1. superproject at the pin (won't move a shared HEAD unless allowed);
#  2. clean tree FIRST, so a dirty superproject is caught before submodule checks
#     (SWIPL_ALLOW_DIRTY=1 overrides for intentional local hacking);
#  3. the WASM-compiled submodules (pins.sh $WASM_SUBMODULES) at the
#     superproject's recorded gitlinks.
verify_superproject_pin "$SWIPL_SRC" "$SWIPL_COMMIT"
verify_clean_tree "$SWIPL_SRC"
verify_submodules "$SWIPL_SRC"
mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR"
# Direct emcmake instead of upstream's `../scripts/configure -y wasm`: the
# flags below ARE that preset (Release, find-root, USE_GMP=OFF, no docs/src;
# emcmake supplies the Emscripten toolchain file; the wasm binary
# self-bootstraps boot.prc + library qlf under node — no native friend), plus
# -DSWIPL_PACKAGE_LIST, which the preset cannot pass. Re-verify the flag set
# against scripts/configure's wasm case on every SWIPL_COMMIT bump.
# cmake honours an explicit SWIPL_PACKAGE_LIST cache var over its package-set
# assembly (cmake/PackageSelection.cmake) and auto-adds hard package deps
# (sgml, for clib+http) — no swipl-devel source is modified.
# NB: cmake 4.2.0–4.3.3 warns "does not support emscripten shared libraries" on
# every probe — harmless, the kernel links statically.
emcmake cmake -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_FIND_ROOT_PATH="$WASM_HOME" -DUSE_GMP=OFF \
  -DINSTALL_DOCUMENTATION=OFF -DINSTALL_PROLOG_SRC=OFF \
  -DSWIPL_PACKAGE_LIST="$WASM_PACKAGE_LIST" -G Ninja ..
ninja swipl-web

# --- 2.6 crosswordsmith-web profile image (payload plan Phases 2+3) ----------
# SWI's wasm build packs its ENTIRE library into swipl-web.data (~1 MB raw for
# ~250 KB actually reachable from the browser qlf) and statically links every
# plugin the configured packages compile. Rather than patching the pinned
# tree, replay the swipl-web link (extracted from ninja itself — not
# hand-maintained, so it always carries exactly the flags/objects of the full
# link it shadows; the grep asserts the expected --preload-file token so a
# cmake reshape fails loudly instead of silently shipping the wrong image)
# with three substitutions:
#   a. preload: the tracked keep-list (preload-profile.txt) assembled into a
#      CLEAN sibling dir (pruning the real wasm-preload/ in place is NOT an
#      option: install_in_wasm_preload() files are ninja outputs and any later
#      ninja run restores them mid-graph)              [Phase 2, .data]
#   b. static extensions: a static_packages.h generated from pins.sh
#      $WASM_STATIC_EXTENSIONS compiles a variant pl-load.c.o listed BEFORE
#      src/libswipl.a, so the linker never pulls the archive's full-registry
#      member; every plugin archive not in the keep-list is dropped from the
#      link inputs (their only referent was that registry)  [Phase 3, .wasm]
#   c. output: -o into the profile dir (fixed name: the only path-sensitive
#      bytes are the PACKAGE_NAME literal, so a constant dir keeps replays
#      byte-identical — verified by double-link)
# The full-registry node swipl.js is untouched: library qlfmake and the value
# goldens still get every plugin there.
log "2.6 profile image: crosswordsmith-web (preload keep-list + static extensions: $WASM_STATIC_EXTENSIONS)"
PROFILE_FILE="$REPO_ROOT/wasm/build/preload-profile.txt"
PRELOAD_SRC="$BUILD_DIR/src/wasm-preload"
PRELOAD_CW="$BUILD_DIR/src/wasm-preload-crosswordsmith-web"
PROFILE_OUT="$BUILD_DIR/src/profile-crosswordsmith-web"
PROFILE_INC="$BUILD_DIR/src/profile-crosswordsmith-web-inc"
rm -rf "$PRELOAD_CW" "$PROFILE_OUT" "$PROFILE_INC"
mkdir -p "$PRELOAD_CW" "$PROFILE_OUT" "$PROFILE_INC"
while IFS= read -r f; do
  case "$f" in ''|\#*) continue ;; esac
  if [ ! -f "$PRELOAD_SRC/$f" ]; then
    echo "preload profile FAILED: '$f' is in preload-profile.txt but not in $PRELOAD_SRC" >&2
    echo "  (stale profile after an SWI pin bump, or the staging step did not run)" >&2
    exit 1
  fi
  mkdir -p "$PRELOAD_CW/$(dirname "$f")"
  cp "$PRELOAD_SRC/$f" "$PRELOAD_CW/$f"
done < "$PROFILE_FILE"

# b1. the trimmed static-extension registry (same shape cmake's
# write_static_extensions emits into src/static_packages.h)
{ for e in $WASM_STATIC_EXTENSIONS; do echo "extern void install_$e(void);"; done
  echo
  echo "static_extension static_extensions[] = {"
  for e in $WASM_STATIC_EXTENSIONS; do echo "  { \"$e\", install_$e },"; done
  echo "  { NULL, NULL }"
  echo "};"
} > "$PROFILE_INC/static_packages.h"

# b2. compile the variant pl-load.c.o against it: take the real compile
# command from ninja, front-run the include path, retarget the output.
CC_CMD="$(cd "$BUILD_DIR" && ninja -t commands src/CMakeFiles/swiplobjs.dir/pl-load.c.o | tail -n1)"
if [ -z "$CC_CMD" ] || [[ "$CC_CMD" != *" -I"* ]]; then
  echo "profile image FAILED: could not extract the pl-load.c.o compile command" >&2
  exit 1
fi
CC_CMD="${CC_CMD/ -I/ -I$PROFILE_INC -I}"
CC_CMD="${CC_CMD//src\/CMakeFiles\/swiplobjs.dir\/pl-load.c.o/$PROFILE_INC/pl-load.c.o}"
( cd "$BUILD_DIR" && bash -c "$CC_CMD" )
[ -f "$PROFILE_INC/pl-load.c.o" ] || { echo "profile image FAILED: no variant pl-load.c.o" >&2; exit 1; }

LINK_CMD="$(cd "$BUILD_DIR" && ninja -t commands swipl-web \
            | grep -F -- "--preload-file $PRELOAD_SRC@swipl" | tail -n1)"
if [ -z "$LINK_CMD" ]; then
  echo "profile image FAILED: could not extract the swipl-web link command" >&2
  echo "  (expected a ninja command containing '--preload-file $PRELOAD_SRC@swipl')" >&2
  exit 1
fi
LINK_CMD="${LINK_CMD/--preload-file $PRELOAD_SRC@swipl/--preload-file $PRELOAD_CW@swipl}"
LINK_CMD="${LINK_CMD/-o src\/swipl-web.js/-o $PROFILE_OUT/swipl-web.js}"
LINK_CMD="${LINK_CMD/ src\/libswipl.a/ $PROFILE_INC/pl-load.c.o src/libswipl.a}"

# b3. keep only the plugin archives in $WASM_STATIC_EXTENSIONS; a plugin whose
# registry entry is gone has no referent and would only bloat the wasm. Token
# filter (the ninja command has no quoted arguments), and count what we drop
# so a cmake reshape that renames archives fails loudly.
TRIMMED_CMD=""
dropped=0
for tok in $LINK_CMD; do
  case "$tok" in
    packages/*.a)
      base="$(basename "$tok" .a)"
      keep=0
      for e in $WASM_STATIC_EXTENSIONS; do [ "$base" = "$e" ] && keep=1; done
      if [ "$keep" = 1 ]; then TRIMMED_CMD+=" $tok"; else dropped=$((dropped+1)); fi
      ;;
    *) TRIMMED_CMD+=" $tok" ;;
  esac
done
LINK_CMD="${TRIMMED_CMD# }"
if [ "$dropped" -eq 0 ]; then
  echo "profile image FAILED: no plugin archive was dropped — link command shape changed?" >&2
  exit 1
fi
( cd "$BUILD_DIR" && bash -c "$LINK_CMD" )
for a in swipl-web.js swipl-web.wasm swipl-web.data; do
  [ -f "$PROFILE_OUT/$a" ] || { echo "profile image FAILED: link produced no $a" >&2; exit 1; }
done
echo "  profile image: $(wc -c < "$PROFILE_OUT/swipl-web.wasm" | tr -d ' ') bytes wasm, $(wc -c < "$PROFILE_OUT/swipl-web.data" | tr -d ' ') bytes data ($dropped plugin archives dropped)"

# --- 3. copy the three web artifacts into wasm/client -----------------------
# The PROFILE image is what ships; the full-image swipl-web.{js,wasm,data}
# stay in $BUILD_DIR/src/ as the SWI-default build (diagnostics, comparison).
log "3. artifacts (crosswordsmith-web profile) -> $CLIENT_DIR"
cp "$PROFILE_OUT/swipl-web.js"   "$CLIENT_DIR/"
cp "$PROFILE_OUT/swipl-web.wasm" "$CLIENT_DIR/"
cp "$PROFILE_OUT/swipl-web.data" "$CLIENT_DIR/"

# --- 4. qcompile the app to crosswordsmith.qlf (WASM word size) -------------
# CRITICAL: produce the qlf with the wasm build's own swipl (run under node), not
# native swipl — a native x86-64 qlf encodes 64-bit VM variants/offsets and would
# load-and-misbehave under wasm32. include(user) folds all seven project modules
# into one self-contained file (no source tree needed at runtime). NODERAWFS lets
# node read the sources by absolute path.
# NB: qcompile writes the .qlf ALONGSIDE ITS SOURCE — i.e. $CLIENT_DIR/solve_browser.qlf
# (basename from the source, independent of CWD) — so rename it in place. (Earlier
# recipes cp'd from $BUILD_DIR/src, which never existed; caught by a clean rebuild.)
log "4. qcompile -> $CLIENT_DIR/crosswordsmith.qlf"
node "$BUILD_DIR/src/swipl.js" \
  -g "qcompile('$CLIENT_DIR/solve_browser.pl', [include(user)]), halt" \
  -t 'halt(1)'
mv "$CLIENT_DIR/solve_browser.qlf" "$CLIENT_DIR/crosswordsmith.qlf"

# --- 5. provenance: build-manifest.json + content-hashed names ---------------
# Ties the four outputs together (buildId, per-artifact sha256, swipl commit +
# submodule SHAs, toolchain pins) and emits <stem>.<sha256:12><ext> copies the
# worker prefers when the manifest is present — a CDN redeploy can no longer
# pair a new js with a long-cached stale wasm/qlf.
log "5. provenance -> $CLIENT_DIR/build-manifest.json"
CLIENT_DIR="$CLIENT_DIR" SWIPL_SRC="$SWIPL_SRC" PRELOAD_PROFILE=crosswordsmith-web \
  "$SCRIPT_DIR/stamp-manifest.sh"

log "done — serve wasm/client/ (see wasm/README.md) and open the page."
