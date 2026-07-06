# Plan: WASM build supply-chain hardening

Status: **Batch 1 landed** (guards in `wasm/build/verify-pin.sh`, zlib
fossils-URL + sha256 — on `main`); **Batch 2 landed** (provenance manifest,
cold-runner bootstrap, committed lockfile — see
[`wasm-supply-chain-batch2.md`](./wasm-supply-chain-batch2.md)) · Drafted
2026-07-06.

Scope: close the still-OPEN build-reproducibility / supply-chain debt tracked in
[`wasm-browser-deployment.md`](./wasm-browser-deployment.md) §10.3. Every change
lands in the single build script
[`wasm/build/build-wasm.sh`](../../wasm/build/build-wasm.sh) (plus a new sourced
guard helper and, for the publish-gated batch, the test toolchain). This
document is the executable record; **it makes no code changes itself.**

Line references are against `wasm/build/build-wasm.sh` as of this draft (116
lines). Recorded facts (versions, SHAs, hashes) were verified live on
2026-07-06 and are cited inline.

---

## 0. Verified facts (checked 2026-07-06, before pinning anything)

| Fact | Value | How verified |
|---|---|---|
| swipl-devel pin | `aa6289399…` = `V10.1.10-17-gaa6289399` | `git -C ~/src/swipl-devel rev-parse HEAD` |
| zlib version | `1.3.2` **exists and resolves** | `curl -fsSLI https://zlib.net/zlib-1.3.2.tar.gz` → 200; fossils path → 200 |
| zlib sha256 | `bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16` | `curl -fsSL <url> \| sha256sum`; **current path and `/fossils/` path are byte-identical** |
| pcre2 version | `pcre2-10.47` **tag exists** | live clone at `~/wasm/pcre2` is at tag `pcre2-10.47` |
| pcre2 pinned commit | `f454e231fe5006dd7ff8f4693fd2b8eb94333429` | `git -C ~/wasm/pcre2 rev-parse HEAD` |
| emsdk activated | tool `6.0.1` (script pin) — but repo clone HEAD is at tag `6.0.2` | `git -C ~/wasm/emsdk describe --tags` → `6.0.2`; script activates `6.0.1` — live drift the `emcc` assert catches |

Config constants live inline at the top of the script:
`SWIPL_COMMIT="aa6289399"` (`build-wasm.sh:25`), `EMSDK_VERSION="6.0.1"`
(`:26`), `ZLIB_VERSION="1.3.2"` (`:27`), `PCRE2_VERSION="10.47"` (`:28`). New
pins (`ZLIB_SHA256`, `PCRE2_COMMIT`) are recorded **inline alongside these** —
consistent with the existing style; a single-file `checksums.txt` is not
justified for two deps.

The pin guard is inline at `build-wasm.sh:70-83`; there is no separate
verify/guard script today. `wasm/test/run_all.sh` re-qcompiles the qlf but never
builds the engine or touches pins, so **all §10.3 changes land in
`build-wasm.sh`** (plus a new sourced helper).

---

## 1. Deliverable: extract guards into `wasm/build/verify-pin.sh`

Create a **new sourced helper** `wasm/build/verify-pin.sh` holding the guard
logic as functions, and have `build-wasm.sh` `source` it and call them. This is
not cosmetic — it is what makes every guard **unit-testable in ~1s without the
~600-step ninja build** (§7). Proposed functions:

- `verify_superproject_pin` — the existing `build-wasm.sh:70-83` logic, moved verbatim.
- `verify_submodules` — item (a).
- `verify_clean_tree` — the dirty-tree guard.
- `verify_pcre2_commit` — item pcre2.
- `verify_emcc_version` — the emcc assertion (load-bearing emsdk guard).

`build-wasm.sh` sources `verify-pin.sh` once near the top (after the config
block, ~`:28`) and calls each function at the point where its guard belongs
(mapped per-item below).

---

## 2. Item (a) — HIGH: pin swipl-devel submodules

**Problem** (`wasm-browser-deployment.md:618-623`): the pin guard checks only the
superproject SHA (`build-wasm.sh:70-72`). The WASM build compiles
`packages/{pcre,http,json,clib,utf8proc,plunit,chr,clpqr,semweb}`, which are git
**submodules**; `git checkout <pin>` does not update them and the guard never
inspects them. Two runs can compile different sources with zero signal.

**What to pin against — no new lockfile needed.** Once the superproject SHA is
asserted `== SWIPL_COMMIT` (already done at `:70-72`), the superproject's
recorded **gitlinks are the manifest**. `git submodule status` compares each
submodule's checked-out commit to that recorded gitlink and prefixes drift with
`+` (checked-out ≠ recorded), `-` (uninitialized), or `U` (conflict); an
in-sync submodule shows a **leading space**. So the assertion is: every
WASM-compiled submodule must show the leading-space state. This auto-tracks a
`SWIPL_COMMIT` bump (no manual SHA list to maintain).

**Scope — the nine compiled submodules** (WASM package set per
`wasm-browser-deployment.md:205` = `clpqr plunit chr clib http semweb pcre
utf8proc` + `json`, which `http` pulls in). Recorded SHAs at the pin (for
audit / optional lockfile only — the assertion does **not** hardcode them):

| submodule | recorded SHA |
|---|---|
| packages/chr | `08f42a653c0d2d13e6d7841fa8889113384f904c` |
| packages/clib | `121c9375ba4323839ba27fbce3a9264e94528732` |
| packages/clpqr | `150ab4d51ee28949b2c7c1e4e97044b81b177d1c` |
| packages/http | `d111640558ea3c0434d227f76e9d1ff4d5a1a007` |
| packages/json | `8da0890d7a1b4c4956056349384f3059816ce885` |
| packages/pcre | `7d561e5eb9ef58dac0b277081de93e991cbe7f8f` |
| packages/plunit | `4de27bb655cb2cb314a6c6434f42aaffa492e5bd` |
| packages/semweb | `cda444fbb8031c24e9829f3a6b7fdd425ef595b1` |
| packages/utf8proc | `99e724163532d2b07ea3d6d250e2d4f2ac73ff7f` |

Deliberately **out of scope**: `packages/zlib`, `packages/ssl`, and the other
30+ submodules. `packages/zlib` is the SWI `library(zlib)` binding, not the C
zlib the build links (that comes from `$WASM_HOME`, §3); the rest are
uninitialized on a normal checkout and asserting on them would false-fail on `-`.

**Where it lands:** call `verify_submodules` immediately after the superproject
pin guard closes (after `build-wasm.sh:83`), before `mkdir -p "$BUILD_DIR"`
(`:84`).

**Assertion (read-only; always runs):**

```bash
# verify_submodules — in verify-pin.sh; called after the superproject pin guard.
verify_submodules() {
  local src="$1"
  local subs="packages/chr packages/clib packages/clpqr packages/http \
packages/json packages/pcre packages/plunit packages/semweb packages/utf8proc"
  local drift
  drift="$(git -C "$src" submodule status $subs | grep -vE '^ ' || true)"
  if [ -n "$drift" ]; then
    echo "ERROR: swipl-devel WASM submodule(s) drifted or uninitialized:" >&2
    echo "$drift" >&2
    echo "  Sync them:  git -C $src submodule update --init -- $subs" >&2
    echo "  (or re-run with SWIPL_ALLOW_CHECKOUT=1 to let this script sync them)" >&2
    return 1
  fi
}
```

**Mutation gated by `SWIPL_ALLOW_CHECKOUT`** (mirrors the superproject
auto-checkout at `:73-75`, which mutates the shared tree). When
`SWIPL_ALLOW_CHECKOUT=1`, run the fix **before** the assert, right after the
superproject `git checkout` at `:75`:

```bash
git -C "$SWIPL_SRC" submodule update --init -- $WASM_SUBMODULES
```

The read-only assert always runs; only the auto-sync is gated. This matches the
"refuse to move a shared HEAD" ethos at `:77-81`.

**Optional lockfile (defense-in-depth, deferred):** a checked-in
`wasm/build/swipl-submodules.lock` of the nine `path SHA` pairs above, asserted
exactly, would additionally catch a *wrong gitlink* (superproject tampered) and
give a human-auditable diff. Cost: regenerate on every `SWIPL_COMMIT` bump.
Recommendation: **skip until Batch 2 / standalone-CI**, where a locally
checked-out superproject can't be trusted; the gitlink-status assert is
sufficient and maintenance-free until then.

---

## 3. Item (b) — HIGH: zlib integrity + pinned URL

**Problem** (`wasm-browser-deployment.md:624-627`): `curl -sL zlib.net/… | tar`
at `build-wasm.sh:50` — no `--fail`, no checksum, and `zlib.net` rotates
non-current releases into `/fossils/`, so this 404s once 1.3.2 is superseded.
Behind an `if [ ! -f "$WASM_HOME/lib/libz.a" ]` cache-skip (`:49`), so the
breakage stays invisible until a from-clean build.

### 3.1 Version-resolution FIRST (before pinning any hash)

The pinned `ZLIB_VERSION="1.3.2"` sits behind the cache-skip, so a stale/wrong
version could already 404 unnoticed. **First action:** confirm the pinned
version resolves, then record the hash of *that* tarball. If the version were
wrong, correcting it is part of this hardening — do not hardcode a hash for an
unverified version.

Resolution procedure (run once at pin time; re-run when bumping the version):

```bash
V=1.3.2
# 1. confirm the version exists at the durable (fossils) path
curl -fsSLI "https://zlib.net/fossils/zlib-$V.tar.gz"        # expect HTTP 200
# 2. record the sha256 OF THE RESOLVED TARBALL
curl -fsSL "https://zlib.net/fossils/zlib-$V.tar.gz" | sha256sum
# 3. (optional) cross-check the current path is byte-identical to fossils
```

**Result of running this on 2026-07-06:** `1.3.2` resolves at both the current
and `/fossils/` paths, and the two are **byte-identical**. Recorded hash:

```
bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16
```

So the version is correct; pin it, the fossils URL, and this hash.

### 3.2 The three fixes

1. **Pin the fossils URL** — `https://zlib.net/fossils/zlib-$ZLIB_VERSION.tar.gz`
   (durable archive path; the bare `/zlib-…` path serves only the current release).
2. **`curl -fsSL`** — `-f` fails on 4xx/5xx (no HTML error body saved to disk),
   `-S` surfaces errors under `-s`, `-L` follows redirects.
3. **Download → `sha256sum -c` → extract** — you cannot checksum a stream already
   piped into `tar`, so the `:50` pipe splits into three steps. Also `rm -rf` the
   extract dir first, closing the lower §10.3 note "zlib branch re-extracts over a
   possibly-dirty source dir (pcre2 branch `rm -rf`s first)"
   (`wasm-browser-deployment.md:638-639`).

**Where the hash is recorded:** inline constant at `build-wasm.sh:27`, next to
`ZLIB_VERSION`:

```bash
ZLIB_VERSION="1.3.2"
ZLIB_SHA256="bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16"
```

**Replacement for `build-wasm.sh:49-55`** (the zlib block):

```bash
if [ ! -f "$WASM_HOME/lib/libz.a" ]; then
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
```

---

## 4. Item pcre2 — MEDIUM: assert the clone's commit SHA

**Problem** (was omitted from the first pass): `build-wasm.sh:56-64` does
`git clone --depth 1 --branch "pcre2-$PCRE2_VERSION"` (`:58`). A git tag is
safer than a rotating tarball but is **movable** and there is no commit-SHA
assertion — a retagged upstream `pcre2-10.47` would compile silently different
source. Same class of gap as zlib, lower priority (tags rarely move).

**Two additions:**

1. **Verify the tag exists** (version-resolution, analogous to zlib §3.1) — at
   pin time: `git ls-remote --tags https://github.com/PCRE2Project/pcre2 pcre2-10.47`
   must return a ref. **Confirmed 2026-07-06:** `pcre2-10.47` exists →
   commit `f454e231fe5006dd7ff8f4693fd2b8eb94333429`.
2. **Assert the resolved commit** after clone. Record inline at `build-wasm.sh:28`:

   ```bash
   PCRE2_VERSION="10.47"
   PCRE2_COMMIT="f454e231fe5006dd7ff8f4693fd2b8eb94333429"
   ```

   Add a `verify_pcre2_commit` call inside the pcre2 block, immediately after the
   clone at `build-wasm.sh:58` (before `mkdir -p build`):

   ```bash
   # verify_pcre2_commit — in verify-pin.sh
   verify_pcre2_commit() {
     local dir="$1" want="$2" got
     got="$(git -C "$dir" rev-parse HEAD)"
     if [ "$got" != "$want" ]; then
       echo "ERROR: pcre2 clone at $got != pinned $want (tag may have moved)." >&2
       return 1
     fi
   }
   ```
   called as `verify_pcre2_commit "$WASM_HOME/pcre2" "$PCRE2_COMMIT"`.

Note the pcre2 block already `rm -rf pcre2` before cloning (`:57`), so no
dirty-dir concern here.

---

## 5. Lower items

### 5.1 Dirty working tree (`verify_clean_tree`)

**Problem** (`wasm-browser-deployment.md:628-630`): a right-SHA-but-dirty tree
builds modified sources and reports success. Add alongside `verify_submodules`
(after `build-wasm.sh:83`), running **first** so a dirty superproject is caught
before submodule checks:

```bash
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
```

`--ignore-submodules=none` makes `--porcelain` also surface submodule commit
drift (belt-and-suspenders with `verify_submodules`, which additionally
distinguishes uninitialized vs. wrong-commit). **New env var `SWIPL_ALLOW_DIRTY`**
gives an intentional escape hatch (mirrors `SWIPL_ALLOW_CHECKOUT`); it MUST be
documented in the same commit (§8).

### 5.2 `emcc --version` assertion — the load-bearing emsdk guard (`verify_emcc_version`)

**Problem** (`wasm-browser-deployment.md:638`): no assertion that the activated
compiler is the pinned one. This is the **primary** emsdk guard: `./emsdk install
6.0.1 && ./emsdk activate 6.0.1` (`build-wasm.sh:43`) resolves the toolchain
deterministically from the version number regardless of the emsdk *repo* HEAD, so
the version assert — not a repo pin — is what catches drift (and the live clone is
at repo tag `6.0.2` while activating `6.0.1`, exactly the case it catches).

Add after `source …/emsdk_env.sh` (`build-wasm.sh:45`):

```bash
verify_emcc_version() {
  local want="$1" got
  got="$(emcc --version | head -n1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1)"
  if [ "$got" != "$want" ]; then
    echo "ERROR: active emcc $got != pinned $want." >&2
    return 1
  fi
}
```
called as `verify_emcc_version "$EMSDK_VERSION"`.

### 5.3 emsdk repo-tag pin — **optional**, demoted

**Problem** (`wasm-browser-deployment.md:638`): `git clone …/emsdk` (`:41`)
floats HEAD. Because `emsdk install/activate <version>` is deterministic (§5.2),
this is cosmetic once `verify_emcc_version` is in place. **Optional** hardening if
desired: `git clone --branch "$EMSDK_VERSION" --depth 1 …` at `:41` (emsdk tags
per SDK version). Not part of the Batch-1 critical set — the `emcc` assert is.

---

## 6. From-clean validation (Batch-1 acceptance gate)

Both deps sit behind cache-skips (`libz.a` present → skip, `:49`;
`libpcre2-8.a` present → skip, `:56`), so a `sha256sum -c` you *add* may **never
execute** on a warm `$WASM_HOME`. The guards must be proven against an **empty
`$WASM_HOME`**:

```bash
WASM_HOME="$(mktemp -d)/wasm-clean" SWIPL_ALLOW_CHECKOUT=0 wasm/build/build-wasm.sh
```

This forces the real download → `curl -fsSL` → `sha256sum -c` → build path for
both deps and proves the pinned URLs + hashes + tags resolve and compile from
scratch. It is also the first concrete step toward the Batch-2
"standalone-CI-runnable" item (§7 of the deployment plan's lower list) —
de-risking it early. Run this once when the guards land; keep the command in the
PR description as the acceptance evidence. (Full run includes the multi-minute
ninja build; that is expected here — this gate is deliberately the *only* place a
full clean build is required. Guard-logic edits themselves are tested cheaply per
§7.)

---

## 7. Testing guard changes WITHOUT a full wasm rebuild

Every Batch-1 guard runs **before** `ninja swipl-web` (`build-wasm.sh:92`) and is
pure shell + git + sha256. Extracting them into `verify-pin.sh` (§1) lets each be
driven standalone in ~1s. Recipes (all avoid the ~600-step compile):

**Submodule guard fires on drift / passes when synced.** Use a throwaway worktree
so the shared tree is untouched:

```bash
git -C ~/src/swipl-devel worktree add /tmp/swtest aa6289399
source wasm/build/verify-pin.sh
git -C /tmp/swtest/packages/http checkout HEAD~1        # induce '+' drift
verify_submodules /tmp/swtest && echo BAD || echo "OK: fired"      # expect fired (exit 1)
git -C /tmp/swtest submodule update --checkout packages/http
verify_submodules /tmp/swtest && echo "OK: passed" || echo BAD     # expect passed (exit 0)
git -C ~/src/swipl-devel worktree remove --force /tmp/swtest
```
(For full isolation from shared submodule state, use a fresh `git clone` instead
of a worktree.)

**Dirty-tree guard.** In the same worktree: `touch /tmp/swtest/src/pl-main.c` →
`verify_clean_tree /tmp/swtest` exits 1; `git -C /tmp/swtest checkout -- .` →
exits 0; `SWIPL_ALLOW_DIRTY=1 verify_clean_tree /tmp/swtest` → exits 0 (escape
hatch).

**zlib checksum guard.** No download or build needed:

```bash
printf tampered > z.tgz
echo "bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16  z.tgz" | sha256sum -c -   # exit 1
```
And the URL fix by a HEAD probe only:
`curl -fsSLI https://zlib.net/fossils/zlib-1.3.2.tar.gz` → 200.

**pcre2 commit guard.** `verify_pcre2_commit ~/wasm/pcre2 <wrong-sha>` → exit 1;
with `f454e231…` → exit 0. (Reuses the existing clone; no re-clone.)

**emcc guard.** Put a stub `emcc` printing a wrong version on `PATH` →
`verify_emcc_version 6.0.1` exits 1; real emcc → exit 0.

**Integration (once, the expensive one):** the §6 from-clean build — the single
place a full compile is warranted.

---

## 8. Doc discipline (same commit as the code)

This repo moves tracker status with the work. In the commit that lands Batch 1:

1. **Tick the closed §10.3 rows** in `wasm-browser-deployment.md:616-639`:
   - submodules unpinned (HIGH) → closed by §2
   - zlib integrity + rotating URL (HIGH) → closed by §3
   - pin guard ignores dirtiness → closed by §5.1
   - no `emcc --version` assertion → closed by §5.2
   - zlib re-extracts over a possibly-dirty dir → closed by §3.2 (`rm -rf`)
   - (pcre2 SHA assertion, §4, is new hardening beyond the listed rows — note it)
   Leave OPEN (Batch 2): test toolchain reproducible (lockfile), artifact
   provenance / cache-busting, standalone-CI-runnable, optional emsdk repo pin.
2. **Document `SWIPL_ALLOW_DIRTY`** next to `SWIPL_ALLOW_CHECKOUT` in the script
   header usage block (`build-wasm.sh:16-19`) and in `wasm/README.md` (the env
   list near line 83).

---

## 9. Sequencing + OQ-5 coupling

**Batch 1 — quick wins + the two HIGHs (land now; independent of publish).** All
in `build-wasm.sh` + `verify-pin.sh`; read-only guards or a URL/hash swap; no
dependency on npm/CDN:

| # | Item | Where |
|---|---|---|
| 1 | zlib fossils URL + `curl -fsSL` + `sha256sum -c` + `rm -rf` (b, HIGH) | `:27`, `:49-55` |
| 2 | submodule gitlink-status assert (a, HIGH) | after `:83` |
| 3 | dirty-tree assert + `SWIPL_ALLOW_DIRTY` | after `:83` |
| 4 | `emcc --version` assert | after `:45` |
| 5 | pcre2 commit-SHA assert | after `:58` |
| — | extract guards into `verify-pin.sh` + guard tests (§7) | new file |
| — | from-clean validation (§6) + doc ticks (§8) | acceptance |

**Batch 2 — publish-blockers (coupled to OQ-5; do before any npm/CDN deploy, not
before).** Per `wasm-sdk-strategy.md:486-488` (OQ-1) npm packaging is deferred and
"**OQ-5/OQ-6 gate the eventual publish**"; per `wasm-sdk-strategy.md:495-496`
OQ-5 *is* "CI qlf rebuild reproducibility (deployment §10.3 'not
standalone-CI-runnable') + the asset-copy-into-tarball / package-size step":

| Item | Notes |
|---|---|
| ✅ **DONE** standalone-CI-runnable | `bootstrap_swipl_src` (verify-pin.sh): clone `swipl-devel` @ pin + `submodule update --init` the nine WASM submodules when `SWIPL_SRC` absent; `SWIPL_REPO_URL` overridable; never mutates a pre-existing tree. §6's from-clean run remains owed at first real CI/CDN use. **This is OQ-5.** |
| ✅ **DONE** artifact provenance manifest + content-hash filenames | `wasm/build/stamp-manifest.sh` (pins single-sourced in `wasm/build/pins.sh`): `build-manifest.json` = buildId + per-artifact sha256 + swipl commit/submodule SHAs + toolchain pins; `<stem>.<sha256:12><ext>` copies; `worker.js` prefers hashed names when the manifest is present (unhashed fallback). Called by `build-wasm.sh` step 5 + `run_all.sh` step 2.5. **OQ-5 asset step.** |
| ✅ **DONE** committed `package-lock.json` + `npm ci` | lockfile un-ignored + committed; Playwright pinned exact (`1.61.1`); `run_all.sh` preflight + wasm/README say `npm ci`. |
| ✅ **DONE** (optional) emsdk repo `--branch` pin | fresh clones use `--branch "$EMSDK_VERSION" --depth 1`; cosmetic given the emcc assert. |
| ✅ **SEEDED** OQ-6 licensing manifest | `wasm/THIRD_PARTY_NOTICES.md` (SWI BSD-2-Clause, zlib, PCRE2, Emscripten glue; recorded future UKACD18 verbatim obligation); verbatim texts inlined at publish. |

There is **no dependency from Batch 1 into Batch 2** — Batch 1 hardens today's
local build and can ship immediately; Batch 2 is only worth building when
packaging is scheduled.

---

## 10. Risks

- **False-positive submodule failures blocking legit builds.** Over-broad scope
  (all 40+ submodules) false-fails on uninitialized `-` lines. Mitigation: scope to
  the nine compiled submodules (§2); `packages/zlib`/`ssl` excluded by design.
- **Lockfile staleness** (if the optional lockfile is adopted): hardcoded SHAs go
  stale on a `SWIPL_COMMIT` bump. Mitigation: prefer the gitlink-status assert
  (auto-tracks bumps); adopt the lockfile only at standalone-CI.
- **Dirty-tree guard blocks intentional local hacking.** Mitigation:
  `SWIPL_ALLOW_DIRTY=1` escape hatch (§5.1), documented (§8).
- **zlib hash/version wrong blocks every clean build.** Mitigation: the §3.1
  version-resolution-FIRST procedure records the hash of the *resolved* tarball;
  verified 2026-07-06. Re-run §3.1 on any version bump.
- **Fossils path not yet populated for a still-current release.** If a bumped
  version isn't in `/fossils/` yet, `-f` 404s. Mitigation: for 1.3.2 both paths
  serve byte-identical content today; the integrity check (the load-bearing HIGH
  fix) is independent of which path is used — keep the current path only if a
  future version isn't in fossils yet, but always keep the sha256.
- **pcre2 tag moved / `git ls-remote` needed at bump.** Mitigation: the
  commit-SHA assert (§4) catches a moved tag; re-record `PCRE2_COMMIT` on bump.
- **emcc version-string format drift across SDK lines.** Mitigation: match the
  pinned `X.Y.Z` substring, not full-line equality (§5.2).
- **Cache-skip hides a broken guard.** A warm `$WASM_HOME` skips the download, so
  an added checksum never runs. Mitigation: the §6 from-clean validation is a
  mandatory Batch-1 acceptance step.

---

## 11. Effort estimate (per item)

| Item | Where | Effort |
|---|---|---|
| (b) zlib fossils URL + `curl -fsSL` + `sha256sum -c` + `rm -rf` | `:27`, `:49-55` | **S** (~30 min; hash already recorded) |
| (a) submodule gitlink-status assert | after `:83` | **S** (~30-45 min incl. worktree test) |
| dirty-tree assert + `SWIPL_ALLOW_DIRTY` | after `:83` | **XS** (~15 min) |
| `emcc --version` assert | after `:45` | **XS** (~15 min) |
| pcre2 commit-SHA assert | after `:58` | **XS** (~15 min; SHA recorded) |
| extract into `verify-pin.sh` + guard tests | new file | **S** (~45 min) |
| from-clean validation (§6) | run once | **S** wall (~few min build) / **XS** author time |
| doc ticks + `SWIPL_ALLOW_DIRTY` docs (§8) | 2 docs | **XS** (~15 min) |
| **Batch 1 total** | — | **~0.5-0.75 day** |
| standalone-CI-runnable (OQ-5) | `:24`, `:70` + clone path | **M** (~0.5-1 day) |
| artifact provenance manifest + content-hash (OQ-5) | new step + loader wiring | **M-L** (~1-1.5 days) |
| committed lockfile + `npm ci` | `.gitignore`, `package.json`, `run_all.sh`, README | **S** (~30 min) |
| **Batch 2 total** | — | **~2-3 days** (schedule with packaging) |
