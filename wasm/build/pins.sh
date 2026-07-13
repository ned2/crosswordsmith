# pins.sh — the single source of truth for every supply-chain pin and for the
# crosswordsmith-web build profile, sourced by build-wasm.sh (which enforces
# them) and stamp-manifest.sh (which records them into
# wasm/client/build-manifest.json). Change a pin HERE, nowhere else.
#
# Bumping a pin re-runs its resolution procedure first (hardening plan §3.1/§4):
# confirm the version/tag resolves, THEN record the hash/commit of what resolved.

SWIPL_COMMIT="aa6289399"                     # our pin (V10.1.10-17-g…, 2026-07-01)
EMSDK_VERSION="6.0.1"                         # npm-swipl-wasm's pin; swipl-devel has no wasm CI
ZLIB_VERSION="1.3.2"
ZLIB_SHA256="bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16"

# The crosswordsmith-web package profile (payload plan Phase 3). Passed to
# cmake as -DSWIPL_PACKAGE_LIST, replacing SWI's 13-package WASM default.
# Rationale for each entry:
#   clib   library(uri) — required by library(wasm); also sha/filesex for the
#          NODE-side tooling (qcompile, value goldens run the full app there)
#   json   library(json) — the worker protocol (atom_json_dict/3)
#   http   NOT used at runtime; required so `library_qlf` can compile the core
#          library/dom.pl, which imports library(http/html_decl). Its Prolog
#          files never enter the preload image (preload-profile.txt gates that).
# cmake's dependency walker auto-adds sgml (required by clib and http), so the
# CONFIGURED set is clib+http+json+sgml — recorded as such in the manifest.
WASM_PACKAGE_LIST="clib;json;http"

# The foreign (C) extensions actually linked into swipl-web — the browser
# runtime's true native surface (build-wasm.sh step 2.6 relinks with a
# static_packages.h generated from THIS list and drops every other plugin
# archive). Everything else that the packages above compile (sgml2pl,
# http_stream, sha4pl, hashstream, md54pl, crypt, memfile, files,
# prolog_stream) stays in the build tree for node-side swipl.js but never
# reaches the browser. Keep in sync with the foreign-library needs of
# preload-profile.txt: uri.pl -> foreign(uri), json.pl -> foreign(json),
# readutil.pl -> foreign(readutil) (soft; it falls back to Prolog, but the
# accelerator is 3KB).
WASM_STATIC_EXTENSIONS="uri readutil json"

# The swipl-devel submodules the WASM profile build compiles: the explicit
# WASM_PACKAGE_LIST above plus sgml (cmake dependency of clib and http).
# Cross-check after a build: ls $SWIPL_SRC/build.wasm/packages (4 + swipl.home).
# Consumed by verify-pin.sh (bootstrap + drift guards) and stamp-manifest.sh
# (provenance). Env-overridable (unset-only, so an intentionally EMPTY override
# sticks) so guard logic can be unit-tested against a tiny local fake remote.
: "${WASM_SUBMODULES=packages/clib packages/http packages/json packages/sgml}"
