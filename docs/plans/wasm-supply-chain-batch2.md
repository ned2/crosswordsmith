# Plan: WASM supply-chain hardening — Batch 2 (publish blockers)

Status: **executed on branch `hardening/wasm-batch2`** · Drafted 2026-07-06.

Scope: the Batch 2 items defined in
[`wasm-supply-chain-hardening.md`](./wasm-supply-chain-hardening.md) §9 — the
OQ-5 publish-blockers ("do before any npm/CDN deploy, not before"), plus the
optional emsdk repo pin and an OQ-6 licensing-manifest seed. Batch 1 (the
guards in [`wasm/build/verify-pin.sh`](../../wasm/build/verify-pin.sh), zlib
fossils URL + sha256) is on `main` and is **not** revisited here.

Tracker rows this closes: `wasm-browser-deployment.md` §10.3, the three
`⏳ OPEN (Batch 2)` bullets. Strategy context:
[`wasm-sdk-strategy.md`](./wasm-sdk-strategy.md) OQ-5 (CI qlf rebuild
reproducibility + asset step) and OQ-6 (redistribution licensing manifest).

---

## 0. Verified facts (2026-07-06, this branch's baseline)

| Fact | Value | How verified |
|---|---|---|
| `wasm/test/package-lock.json` | exists on the dev machine, **untracked** (ignored by `wasm/test/.gitignore:2`) | `git check-ignore -v` |
| lockfile resolutions | playwright `1.61.1` (spec `^1.55.0`), playwright-core `1.61.1`, typescript `5.9.3`, fsevents `2.3.2` (darwin-only optional) | read the lockfile |
| artifact consumers | `wasm/client/worker.js` is the ONLY loader of `swipl-web.{js,wasm,data}` + `crosswordsmith.qlf` (importScripts + `locateFile` + `consult`); `probe.html` references `./swipl-web.js` directly (diagnostic page); the SDK facade fetches nothing itself | grep `swipl-web\|qlf` across `wasm/` |
| re-stamp points | `build-wasm.sh` steps 3–4 emit the artifacts; `wasm/test/run_all.sh` steps 1–2 re-stage + re-qcompile (so any manifest must be re-stampable there too) | read both scripts |
| pins location | `SWIPL_COMMIT`/`EMSDK_VERSION`/`ZLIB_*`/`PCRE2_*` inline at `build-wasm.sh` config block | read the script |

---

## 1. Deliverable (a) — artifact provenance manifest + content-hashed names

**Problem** (`wasm-browser-deployment.md` §10.3): the four build outputs
(`swipl-web.{js,wasm,data}`, `crosswordsmith.qlf`) carry no build id tying them
together; a redeploy that changes one while a CDN long-caches the rest
"loads-and-misbehaves" (the qlf word-size warning's own failure mode, §6).

**Design — manifest as the only mutable name.**

1. New script **`wasm/build/stamp-manifest.sh`** (standalone, ~1s, no build):
   - sha256s each of the four artifacts in `$CLIENT_DIR`;
   - copies each to a content-hashed sibling `<stem>.<sha256:12><ext>`
     (e.g. `swipl-web.1a2b3c4d5e6f.wasm`), removing stale hashed copies first;
   - writes **`wasm/client/build-manifest.json`**:

     ```jsonc
     { "v": 1,
       "buildId": "<sha256 over the 4 artifact hashes>",
       "swipl":   { "commit": "<HEAD of $SWIPL_SRC>",
                    "submodules": { "packages/chr": "<sha>", … } },  // the 9 compiled ones
       "toolchain": { "emsdk": "6.0.1", "zlib": "1.3.2",
                      "zlibSha256": "…", "pcre2": "10.47", "pcre2Commit": "…" },
       "artifacts": { "swipl-web.js": { "hashed": "swipl-web.<h12>.js",
                                        "sha256": "…", "bytes": N }, … } }
     ```

   - **No timestamps**: identical inputs ⇒ byte-identical manifest (buildId is
     content-derived, deterministic).
   - Provenance degrades explicitly: if `$SWIPL_SRC` is not a git tree the
     `swipl` fields are `null` (visible, never fabricated).
2. Pins move from `build-wasm.sh`'s config block into a sourced
   **`wasm/build/pins.sh`** so `stamp-manifest.sh` and `build-wasm.sh` share
   one source of truth (same extraction ethos as `verify-pin.sh`).
3. **Wiring:** `build-wasm.sh` gains step 5 (stamp after qcompile);
   `run_all.sh` re-stamps after its stage/qcompile steps so the battery tests
   the manifest-driven load path.
4. **Loader (`wasm/client/worker.js`):** reads `./build-manifest.json` at
   worker start via **synchronous XHR** (legal and non-janky off the main
   thread; `importScripts` is synchronous so an async fetch can't precede it
   at top level). When present, every artifact resolves through
   `manifest.artifacts[name].hashed` — importScripts, `locateFile`, and the
   qlf consult; when absent (hand-staged dev dir) it falls back to the
   unhashed names, so nothing new is required for local hacking.
   `probe.html` (diagnostic) keeps using unhashed names — the originals stay
   in place next to the hashed copies.
5. `wasm/client/.gitignore` grows the hashed patterns + `build-manifest.json`
   (all build outputs, same policy as the artifacts themselves).

**Deploy story this buys:** a CDN deploy ships the hashed set + `worker.js` +
`build-manifest.json`; only the latter two need short/no cache. A partial
deploy can no longer pair a new js with a stale wasm — the names differ per
build and the manifest records the expected sha256 of each.

Out of scope (noted, not silently dropped): `worker.js` itself is a tracked
source file, versioned by the app deploy, not a build output — it is the
manifest's *consumer*, so it cannot itself be named by the manifest. The OQ-7
qlf-hash-in-capabilities echo stays at packaging time.

## 2. Deliverable (b) — standalone-CI-runnable (cold-runner bootstrap)

**Problem:** `build-wasm.sh` assumes `~/src/swipl-devel` exists at the pin; a
cold CI runner has nothing there. **This is OQ-5's "CI qlf rebuild" leg.**

**Design:** new `bootstrap_swipl_src SRC COMMIT` in `verify-pin.sh` (kept with
the other supply-chain functions so it unit-tests the same way), called at the
top of build-wasm.sh step 2, *before* the read-only guards:

- `$SWIPL_SRC` **absent** → `git clone $SWIPL_REPO_URL` (new env override,
  default `https://github.com/SWI-Prolog/swipl-devel.git`), `git checkout
  $SWIPL_COMMIT`, `git submodule update --init -- $WASM_SUBMODULES` (the nine
  compiled submodules; **not** shallow — a shallow tip need not contain the
  recorded gitlinks).
- exists but not a git checkout → loud error.
- exists as a git tree → no-op; the existing guards own it (including their
  `SWIPL_ALLOW_CHECKOUT`-gated refusal to move a shared HEAD — the bootstrap
  never mutates a pre-existing tree).
- `WASM_SUBMODULES` becomes env-overridable so the function is testable
  against a tiny local fake remote (no network, no real clone).

## 3. Deliverable (c) — committed npm lockfile + `npm ci`

- Delete `package-lock.json` from `wasm/test/.gitignore`; commit the lockfile
  (carried over from the dev machine where it was generated — resolutions in
  §0).
- Pin Playwright **exact** in `package.json` (`1.61.1`, the lock's resolved
  version); update the lock's root spec to match so `npm ci` stays in sync.
  (`typescript` stays caret — the committed lock pins its resolution at 5.9.3
  regardless; Playwright is the one that drives the browser protocol and was
  the row's named risk.)
- `run_all.sh` preflight fix-it text + `wasm/README.md` prereq: `npm install`
  → `npm ci` (reproducible, lockfile-exact, fails on drift).

## 4. Folded in

- **emsdk repo `--branch` pin** (hardening plan §5.3, "optional/cosmetic"):
  `git clone --branch "$EMSDK_VERSION" --depth 1` — emsdk tags per SDK
  version; `verify_emcc_version` remains the load-bearing guard. Only affects
  fresh `$WASM_HOME`s.
- **OQ-6 seed:** new [`wasm/THIRD_PARTY_NOTICES.md`](../../wasm/THIRD_PARTY_NOTICES.md)
  — the redistribution manifest for what `wasm/client/` actually ships
  (SWI-Prolog BSD-2-Clause, zlib, PCRE2 BSD-3-Clause, Emscripten MIT glue) plus
  the recorded **future** obligation: UKACD18's verbatim freeware licence must
  ship the day `fill` browserifies. Marked draft-until-publish; OQ-6 closes at
  npm packaging by inlining verbatim texts. The manifest gains a `licenses`
  pointer to it.

## 5. Verification (what runs, what deliberately does not)

| Check | Cost | Gate |
|---|---|---|
| `stamp-manifest.sh` against a synthetic `$CLIENT_DIR` (4 fake artifacts) — manifest fields, hashed copies, stale-hash cleanup, null-provenance degrade | ~1s | (a) |
| `bootstrap_swipl_src` against a local fake remote (clone+checkout path, not-a-git-repo error path, existing-tree no-op) | ~1s | (b) |
| every Batch-1 guard still green (`verify_superproject_pin`/`clean_tree`/`submodules`/`pcre2`/`emcc` driven standalone per hardening plan §7) | ~1s each | regression |
| `bash -n` + shellcheck-by-eye on the three shell files | ~1s | syntax |
| `npm ci` in `wasm/test` (network/cache permitting) | ~1min | (c) |
| `make test` (native suite) | ~min | no regression outside wasm/ |
| **NOT run:** full `ninja swipl-web` build, `make test-wasm` (needs the swipl-devel wasm build tree + Chrome), a real swipl-devel clone from GitHub | hours/network | recorded in the report; the §6-style from-clean acceptance run remains owed **when a real CDN deploy or CI runner is first stood up** |

## 6. Doc discipline (same commit as the code)

- `wasm-browser-deployment.md` §10.3: tick the three Batch-2 rows (+ status
  header updated to "Batch 2 closed").
- `wasm-supply-chain-hardening.md`: Batch-2 table rows marked done with
  landing points; header status line updated.
- `wasm/README.md`: `npm ci`, the manifest + bootstrap story, notices file.
- `docs/STATUS.md`: no wasm row exists — nothing to move (checked).
