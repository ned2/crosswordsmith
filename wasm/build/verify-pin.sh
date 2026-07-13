#!/usr/bin/env bash
# verify-pin.sh — sourced supply-chain guards for build-wasm.sh.
#
# Every guard here runs BEFORE the ~600-step `ninja swipl-web` compile and is
# pure shell + git + sha256, so each can be driven standalone in ~1s without a
# wasm rebuild (see docs/plans/wasm-supply-chain-hardening.md §7). build-wasm.sh
# sources this file once and calls each function at the point its guard belongs.
#
# Each function returns non-zero on a failed check; build-wasm.sh runs under
# `set -e`, so a failed guard aborts the build.

# $WASM_SUBMODULES — the thirteen submodules the WASM build compiles — lives in
# pins.sh with the other supply-chain pins (single source of truth, shared with
# stamp-manifest.sh). Source it here too so these guards keep working when this
# file is sourced standalone for the ~1s guard tests (idempotent re-source when
# build-wasm.sh already did). The set deliberately excludes packages/ssl and
# the 30+ submodules uninitialised on a normal checkout (asserting on them
# false-fails on '-').
# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/pins.sh"

# bootstrap_swipl_src SRC COMMIT — the cold-runner / standalone-CI path (Batch 2,
# docs/plans/wasm-supply-chain-batch2.md §2). When SRC does not exist: clone
# swipl-devel (override the remote with SWIPL_REPO_URL), check out the pin, and
# init the WASM-compiled submodules. NOT shallow — a shallow tip need not
# contain the pinned commit or the recorded gitlinks. A pre-existing tree is
# never mutated here (the read-only guards below own it, with their
# SWIPL_ALLOW_CHECKOUT-gated refusal to move a shared HEAD); a non-git SRC is a
# loud error rather than a surprise clone-over.
bootstrap_swipl_src() {
  local src="$1" want="$2"
  if [ -e "$src" ]; then
    if [ ! -d "$src/.git" ]; then
      echo "ERROR: $src exists but is not a git checkout — refusing to bootstrap over it." >&2
      echo "  Point SWIPL_SRC elsewhere, or remove it to let this script clone the pin." >&2
      return 1
    fi
    return 0
  fi
  echo "  $src absent — cloning swipl-devel @ $want (cold runner)"
  git clone "${SWIPL_REPO_URL:-https://github.com/SWI-Prolog/swipl-devel.git}" "$src"
  git -C "$src" checkout "$want"
  # shellcheck disable=SC2086
  git -C "$src" submodule update --init -- $WASM_SUBMODULES
}

# verify_superproject_pin SRC COMMIT — assert SRC's HEAD is the pinned commit.
# QLF/boot are word-size-specific and this tree is SHARED with other work, so we
# do NOT silently move HEAD: verify the pin, and only check it out (plus sync the
# compiled submodules the checkout drags along) when SWIPL_ALLOW_CHECKOUT=1.
verify_superproject_pin() {
  local src="$1" want="$2" head_sha
  head_sha="$(git -C "$src" rev-parse HEAD)"
  if ! git -C "$src" merge-base --is-ancestor "$want" HEAD 2>/dev/null \
     || [ "$(git -C "$src" rev-parse "$want")" != "$head_sha" ]; then
    if [ "${SWIPL_ALLOW_CHECKOUT:-0}" = "1" ]; then
      echo "  checking out $want (SWIPL_ALLOW_CHECKOUT=1)"
      git -C "$src" checkout "$want"
      # the checkout moved the recorded gitlinks — sync the compiled submodules
      # too so verify_submodules sees them in-sync (mutates the shared tree).
      git -C "$src" submodule update --init -- $WASM_SUBMODULES
    else
      echo "ERROR: $src is not at $want (HEAD=$head_sha)." >&2
      echo "  This tree may be shared with other work — refusing to move HEAD." >&2
      echo "  Fix it yourself:  git -C $src checkout $want" >&2
      echo "  or re-run with:   SWIPL_ALLOW_CHECKOUT=1 $0" >&2
      return 1
    fi
  fi
}

# verify_clean_tree SRC — refuse to build a dirty tree (right sha, modified
# sources, "success"). --ignore-submodules=none also surfaces submodule commit
# drift (belt-and-suspenders with verify_submodules). SWIPL_ALLOW_DIRTY=1 is the
# intentional escape hatch for local hacking (mirrors SWIPL_ALLOW_CHECKOUT).
verify_clean_tree() {
  local src="$1"
  if [ "${SWIPL_ALLOW_DIRTY:-0}" = "1" ]; then return 0; fi
  if [ -n "$(git -C "$src" status --porcelain=v1 --ignore-submodules=none)" ]; then
    echo "ERROR: $src has uncommitted changes — refusing to build modified sources." >&2
    git -C "$src" status --porcelain=v1 --ignore-submodules=none >&2
    echo "  Override intentionally with:  SWIPL_ALLOW_DIRTY=1 $0" >&2
    return 1
  fi
}

# verify_submodules SRC — assert every WASM-compiled submodule is at the gitlink
# the (already pin-verified) superproject records. `git submodule status`
# prefixes drift with '+' (checked-out != recorded), '-' (uninitialised) or 'U'
# (conflict); an in-sync submodule shows a leading space. So: fail on any line
# that is not leading-space. Auto-tracks a SWIPL_COMMIT bump (no SHA list here).
#
# SWIPL_ALLOW_CHECKOUT=1 makes this guard SYNC a drift rather than just fail,
# mirroring verify_superproject_pin's opt-in mutation of the shared tree. This
# also covers the already-at-pin case: verify_superproject_pin only syncs the
# submodules it drags along a superproject checkout, so when HEAD is already the
# pin but a submodule has drifted, THIS is the function that honours the flag.
verify_submodules() {
  local src="$1" drift
  drift="$(git -C "$src" submodule status $WASM_SUBMODULES | grep -vE '^ ' || true)"
  if [ -n "$drift" ]; then
    if [ "${SWIPL_ALLOW_CHECKOUT:-0}" = "1" ]; then
      echo "  syncing drifted WASM submodule(s) (SWIPL_ALLOW_CHECKOUT=1)"
      git -C "$src" submodule update --init -- $WASM_SUBMODULES
      drift="$(git -C "$src" submodule status $WASM_SUBMODULES | grep -vE '^ ' || true)"
      if [ -n "$drift" ]; then
        echo "ERROR: WASM submodule(s) still drifted after sync:" >&2
        echo "$drift" >&2
        return 1
      fi
    else
      echo "ERROR: swipl-devel WASM submodule(s) drifted or uninitialized:" >&2
      echo "$drift" >&2
      echo "  Sync them:  git -C $src submodule update --init -- $WASM_SUBMODULES" >&2
      echo "  (or re-run with SWIPL_ALLOW_CHECKOUT=1 to let this script sync them)" >&2
      return 1
    fi
  fi
}

# verify_pcre2_commit DIR WANT — assert the cloned pcre2 tag resolved to the
# pinned commit. A git tag is safer than a rotating tarball but is movable and
# carries no commit assertion; a retagged upstream pcre2-<ver> would otherwise
# compile silently different source.
verify_pcre2_commit() {
  local dir="$1" want="$2" got
  got="$(git -C "$dir" rev-parse HEAD)"
  if [ "$got" != "$want" ]; then
    echo "ERROR: pcre2 clone at $got != pinned $want (tag may have moved)." >&2
    return 1
  fi
}

# verify_emcc_version WANT — assert the ACTIVE emcc is the pinned toolchain. This
# is the load-bearing emsdk guard: `emsdk install/activate <version>` resolves
# the toolchain from the version number regardless of the emsdk repo HEAD, so the
# version assert — not a repo-tag pin — is what catches drift. Match the X.Y.Z
# substring (not full-line equality) to survive version-string format drift.
verify_emcc_version() {
  local want="$1" got
  got="$(emcc --version | head -n1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1)"
  if [ "$got" != "$want" ]; then
    echo "ERROR: active emcc $got != pinned $want." >&2
    return 1
  fi
}
