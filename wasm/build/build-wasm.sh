#!/usr/bin/env bash
# build-wasm.sh — reproducible build of the crosswordsmith browser bundle.
#
# Turns the prose recipe in docs/plans/wasm-browser-deployment.md §2 into one
# runnable script. It (0) activates emsdk, (1) stages zlib+pcre2 into $WASM_HOME,
# (2) builds swipl-web from our pinned swipl-devel commit, (3) copies the three
# web artifacts into wasm/client/, and (4) qcompiles the app into
# wasm/client/crosswordsmith.qlf using the *wasm* swipl (correct pointer size).
#
# Everything it produces under wasm/client/ is gitignored (see .gitignore) — the
# repo tracks this recipe, not the 3.9 MB of binaries.
#
# First run took a couple of minutes for configure (hundreds of cross-compile
# probes) + ~600 compile steps. Re-runs skip already-staged deps.
#
# Usage:
#   wasm/build/build-wasm.sh            # full build (verifies the swipl-devel pin)
#   SWIPL_ALLOW_CHECKOUT=1 …            # allow the script to `git checkout` the pin
#   WASM_HOME=… SWIPL_SRC=… …          # override locations
set -euo pipefail

# --- configuration (override via env) ---------------------------------------
: "${WASM_HOME:=$HOME/wasm}"                 # staging prefix for cross-compiled deps
: "${SWIPL_SRC:=$HOME/src/swipl-devel}"      # the swipl-devel checkout to build from
SWIPL_COMMIT="aa6289399"                     # our pin (V10.1.10-17-g…, 2026-07-01)
EMSDK_VERSION="6.0.1"                         # npm-swipl-wasm's pin; swipl-devel has no wasm CI
ZLIB_VERSION="1.3.2"
PCRE2_VERSION="10.47"

# repo root = two levels up from this script (wasm/build/ -> repo root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CLIENT_DIR="$REPO_ROOT/wasm/client"
BUILD_DIR="$SWIPL_SRC/build.wasm"

log() { printf '\n\033[1;36m== %s\033[0m\n' "$*"; }

# --- 0. toolchain: emsdk -----------------------------------------------------
log "0. emsdk $EMSDK_VERSION"
if [ ! -d "$WASM_HOME/emsdk" ]; then
  git clone https://github.com/emscripten-core/emsdk "$WASM_HOME/emsdk"
fi
( cd "$WASM_HOME/emsdk" && ./emsdk install "$EMSDK_VERSION" && ./emsdk activate "$EMSDK_VERSION" )
# shellcheck disable=SC1091
source "$WASM_HOME/emsdk/emsdk_env.sh"

# --- 1. stage wasm-side deps into $WASM_HOME --------------------------------
log "1. deps: zlib $ZLIB_VERSION + pcre2 $PCRE2_VERSION -> $WASM_HOME"
if [ ! -f "$WASM_HOME/lib/libz.a" ]; then
  ( cd "$WASM_HOME" && curl -sL "https://zlib.net/zlib-$ZLIB_VERSION.tar.gz" | tar xz
    cd "zlib-$ZLIB_VERSION" && emconfigure ./configure --static --prefix="$WASM_HOME"
    EMCC_CFLAGS=-Wno-deprecated-non-prototype emmake make && emmake make install )
else
  echo "  libz.a present — skip"
fi
if [ ! -f "$WASM_HOME/lib/libpcre2-8.a" ]; then
  ( cd "$WASM_HOME" && rm -rf pcre2
    git clone --depth 1 --branch "pcre2-$PCRE2_VERSION" https://github.com/PCRE2Project/pcre2
    cd pcre2 && mkdir -p build && cd build
    emcmake cmake -DCMAKE_INSTALL_PREFIX="$WASM_HOME" -DPCRE2GREP_SUPPORT_JIT=OFF -G Ninja ..
    ninja && ninja install )
else
  echo "  libpcre2-8.a present — skip"
fi

# --- 2. build swipl-web from OUR pinned commit ------------------------------
log "2. swipl-web from $SWIPL_SRC @ $SWIPL_COMMIT"
# QLF/boot are word-size-specific and this tree is SHARED with other work, so we
# do NOT silently move HEAD. Verify the pin; only check out if explicitly allowed.
HEAD_SHA="$(git -C "$SWIPL_SRC" rev-parse HEAD)"
if ! git -C "$SWIPL_SRC" merge-base --is-ancestor "$SWIPL_COMMIT" HEAD 2>/dev/null \
   || [ "$(git -C "$SWIPL_SRC" rev-parse "$SWIPL_COMMIT")" != "$HEAD_SHA" ]; then
  if [ "${SWIPL_ALLOW_CHECKOUT:-0}" = "1" ]; then
    echo "  checking out $SWIPL_COMMIT (SWIPL_ALLOW_CHECKOUT=1)"
    git -C "$SWIPL_SRC" checkout "$SWIPL_COMMIT"
  else
    echo "ERROR: $SWIPL_SRC is not at $SWIPL_COMMIT (HEAD=$HEAD_SHA)." >&2
    echo "  This tree may be shared with other work — refusing to move HEAD." >&2
    echo "  Fix it yourself:  git -C $SWIPL_SRC checkout $SWIPL_COMMIT" >&2
    echo "  or re-run with:   SWIPL_ALLOW_CHECKOUT=1 $0" >&2
    exit 1
  fi
fi
mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR"
# -y MUST be the first arg (without a tty, configure otherwise loops forever on a
# "Run Cmake?" prompt). The wasm profile sets the toolchain, USE_GMP=OFF,
# CMAKE_FIND_ROOT_PATH=$WASM_HOME, INSTALL_PROLOG_SRC=OFF, Release, and NO native
# friend (the wasm binary self-bootstraps boot.prc + library qlf under node).
# NB: cmake 4.2.0–4.3.3 warns "does not support emscripten shared libraries" on
# every probe — harmless, the kernel links statically.
../scripts/configure -y wasm
ninja swipl-web

# --- 3. copy the three web artifacts into wasm/client -----------------------
log "3. artifacts -> $CLIENT_DIR"
cp "$BUILD_DIR/src/swipl-web.js"   "$CLIENT_DIR/"
cp "$BUILD_DIR/src/swipl-web.wasm" "$CLIENT_DIR/"
cp "$BUILD_DIR/src/swipl-web.data" "$CLIENT_DIR/"

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

log "done — serve wasm/client/ (see wasm/README.md) and open the page."
