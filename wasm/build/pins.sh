# pins.sh — the single source of truth for every supply-chain pin, sourced by
# build-wasm.sh (which enforces them) and stamp-manifest.sh (which records them
# into wasm/client/build-manifest.json). Change a pin HERE, nowhere else.
#
# Bumping a pin re-runs its resolution procedure first (hardening plan §3.1/§4):
# confirm the version/tag resolves, THEN record the hash/commit of what resolved.

SWIPL_COMMIT="aa6289399"                     # our pin (V10.1.10-17-g…, 2026-07-01)
EMSDK_VERSION="6.0.1"                         # npm-swipl-wasm's pin; swipl-devel has no wasm CI
ZLIB_VERSION="1.3.2"
ZLIB_SHA256="bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16"
PCRE2_VERSION="10.47"
PCRE2_COMMIT="f454e231fe5006dd7ff8f4693fd2b8eb94333429"

# The thirteen swipl-devel submodules the WASM build compiles: the eight in
# cmake/PackageSelection.cmake's SWIPL_PACKAGE_LIST_WASM (chr clib clpqr http
# pcre plunit semweb utf8proc) plus the dependency closure cmake enables with
# them — json + sgml (http, clib) and RDF + nlp + zlib (semweb). packages/zlib
# here is the SWI library(zlib) *binding*; the C zlib it wraps is staged into
# $WASM_HOME from the ZLIB_VERSION tarball above — both are pinned, separately.
# Cross-check on a pin bump: ls $SWIPL_SRC/build.wasm/packages (13 + swipl.home).
# Consumed by verify-pin.sh (bootstrap + drift guards) and stamp-manifest.sh
# (provenance). Env-overridable (unset-only, so an intentionally EMPTY override
# sticks) so guard logic can be unit-tested against a tiny local fake remote.
: "${WASM_SUBMODULES=packages/RDF packages/chr packages/clib packages/clpqr \
packages/http packages/json packages/nlp packages/pcre packages/plunit \
packages/semweb packages/sgml packages/utf8proc packages/zlib}"
