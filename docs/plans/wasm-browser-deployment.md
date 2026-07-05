# Plan: crosswordsmith in the browser (SWI-Prolog WASM)

Status: **strategy / pre-implementation** · Drafted 2026-07-05.

Goal: run the crosswordsmith solver in the browser under SWI-Prolog's WASM
runtime, at the **same SWI version we build natively**, without regressing the
CLI. This doc is the strategy; the runnable proof-of-shape is the spike under
[`spikes/wasm-browser/`](../../spikes/wasm-browser/).

Grounded in: our own `swipl-devel` checkout at `~/src/swipl-devel`
(`git describe` = `V10.1.10-17-gaa6289399`, HEAD 2026-07-01), the version-matched
WASM manual under [`docs/reference/swi-manual/`](../reference/swi-manual/)
(`wasm-loading.md`, `wasm-calling.md`, `wasm-js-call.md`), the research note
[`docs/research/swi-vm-wasm-performance.md`](../research/swi-vm-wasm-performance.md),
a 2026-07 web sweep of the SWI build docs / Discourse / `npm-swipl-wasm`, and a
native spike that exercises the file-free entry point end to end. Source URLs
are collected at the foot.

---

## 1. Bottom line

This is a **low-risk port**. Three facts make it easy:

1. **Version parity is free.** The full Emscripten toolchain lives *inside our
   own source tree* (`cmake/port/Emscripten.cmake`, `cmake/EmscriptenTargets.cmake`,
   `src/wasm/`, `library/wasm.pl`). We compile the browser target from the exact
   checkout we already build natively — there is no separate version to mirror.
2. **The pinned commit is already fast.** `aa6289399` (2026-07-01) post-dates
   the Sept-2025 fix that restored WASM performance (`7b4b138ba`, "Move
   setjmp() out of the VM main function"); `git merge-base --is-ancestor`
   confirms it. Building from our tree yields the ~2.7–3×-native path with **no
   special flags** — see §2 and the correction in §8.
3. **The dependency surface is fully covered and the porting seam is tiny.**
   Every library the solver uses is available under WASM (§4); the only code
   work is swapping two file-I/O predicates for in-memory I/O (§5).

**Decision summary:** self-build `swipl-web` from our tree (npm can't match our
commit); ship the app as one `.qlf`; run the solver in a **Web Worker** with
**input via query binding, output via JSON**; reuse one instance; cancel by
`worker.terminate()`.

## 2. Build — exact-version parity from our tree

`$WASM_HOME` is exported to `$HOME/wasm` (the staging prefix, see below).
`emsdk` is **not yet installed** on this box (`$WASM_HOME/emsdk` absent, no
`emcc`); `ninja` is present. Installing emsdk + staging deps is the one real
prerequisite.

```bash
# 0. toolchain — emsdk 6.0.1 (SWI's own CI pin; anything ≥4.0.15 works post-fix,
#    3.1.37 in the stale official doc is obsolete, emsdk 4.0.0 briefly broke CI)
git clone https://github.com/emscripten-core/emsdk $WASM_HOME/emsdk
cd $WASM_HOME/emsdk && ./emsdk install 6.0.1 && ./emsdk activate 6.0.1 && source ./emsdk_env.sh

# 1. stage wasm-side deps into $WASM_HOME — see "Dependency staging" below

# 2. configure + build from OUR source tree
cd ~/src/swipl-devel && mkdir -p build.wasm && cd build.wasm
emcmake cmake -DCMAKE_TOOLCHAIN_FILE=$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake \
  -DCMAKE_FIND_ROOT_PATH=$WASM_HOME -DCMAKE_BUILD_TYPE=Release \
  -DUSE_GMP=OFF -DINSTALL_QLF=ON -DINSTALL_PROLOG_SRC=OFF -DINSTALL_DOCUMENTATION=OFF \
  -DSWIPL_NATIVE_FRIEND=$HOME/src/swipl-devel/build -G Ninja ..
ninja swipl-web        # or: ninja swipl-web swipl-bundle
```

Prefer the in-tree `scripts/configure wasm` profile (it encodes these flags —
and defaults `WASM_HOME` to `$HOME/wasm`, writing it into `.envrc`) over the
**official `WebAssembly.md`, which is ~3 years stale** (still names Emscripten
3.1.37, a non-`emcmake` cmake line).

### Dependency staging (`$WASM_HOME`)

`$WASM_HOME` is just the **staging prefix** for emscripten-cross-compiled deps:
the `wasm` profile passes `-DCMAKE_FIND_ROOT_PATH=$WASM_HOME` so CMake's
`find_package` resolves them there. The path is arbitrary (`export WASM_HOME=…`
to move it). Using a *separate* prefix rather than the emsdk sysroot is the
**idiomatic, deliberate** choice — it survives emsdk upgrades and `emcc
--clear-cache`, and stays reproducible. The `npm-swipl-wasm` reference build
does exactly this into `/wasm`, so we are matching the canonical convention.

Two facts fix the approach:

- **zlib is mandatory** — SWI core does `find_package(ZLIB REQUIRED)`, resolved
  via `find_package`. The emscripten built-in *port* (`--use-port=zlib`) is
  awkward here: ports build lazily at link time and aren't seen by configure-
  time `find_package` unless you `embuilder build zlib` first. So a manual build
  into `$WASM_HOME` is the smooth path.
- **pcre2 has no emscripten port at all** — manual build (or a package manager)
  regardless.

Least-friction path = mirror npm/SWI, building both into `$WASM_HOME` (versions
per `npm-swipl-wasm/build-config.json`; confirm exact tags at build time):

```bash
# zlib 1.3.2 — autotools
cd $WASM_HOME && curl -sL https://zlib.net/zlib-1.3.2.tar.gz | tar xz
cd zlib-1.3.2 && emconfigure ./configure --static --prefix=$WASM_HOME
EMCC_CFLAGS=-Wno-deprecated-non-prototype emmake make && emmake make install

# pcre2 10.47 — cmake; JIT off (no JIT under wasm)
cd $WASM_HOME && git clone --depth 1 --branch pcre2-10.47 https://github.com/PCRE2Project/pcre2
cd pcre2 && mkdir build && cd build
emcmake cmake -DCMAKE_INSTALL_PREFIX=$WASM_HOME -DPCRE2GREP_SUPPORT_JIT=OFF -G Ninja ..
ninja && ninja install
```

pcre2 is only pulled by the `pcre` package; it *could* be dropped by overriding
the WASM package list, but toggles are per-set (no clean per-package flag), so
building it is simpler than trimming for a first build. Other WASM-set packages
(`utf8proc`, `clib`, `http`, `json`, …) need no external deps.

Landscape, for the record: **ports** (zero-install but zlib-only-here +
`find_package` friction), **manual prefix** (this; sysroot-vs-separate — we use
separate), and **package managers** (emscripten-forge ships both zlib & pcre2;
vcpkg's `wasm32-emscripten` is community/flaky; Conan viable) — the last is
overkill for two deps.

**Flag decisions:**

| Flag | Setting | Why |
|---|---|---|
| `USE_GMP` | **OFF** | LibBF bignums; the WASM norm. arrange/fill are integer-only, so the weaker LibBF rational arithmetic is irrelevant. |
| `INSTALL_QLF` / `INSTALL_PROLOG_SRC` | **ON / OFF** | qlf-only image → smaller `.data`, no source tracing. |
| `WASM_EXCEPTIONS` | **OFF (default)** | *Alternative* to the setjmp fix, ~30% slower; superseded by the fix already in our pin (§8). Not additive. |
| `SWIPL_NATIVE_FRIEND` | our native build dir | Skips the slow node-emulator bootstrap (the build is inherently two-stage: a node/native swipl generates `boot.prc`+`.qlf`, then those are preloaded into the browser target). Friend must be the **same commit**. |
| stack/mem link flags | in-tree defaults | `STACK_SIZE=1M`, `STACK_OVERFLOW_CHECK`, `ALLOW_MEMORY_GROWTH`, `WASM_BIGINT`, `LZ4` are already set in `EmscriptenTargets.cmake`. No manual tuning. |
| `VMI_FUNCTIONS` | **default (OFF for Emscripten)** | The tree forces it off for WASM; the old note's "ON, >10% faster" is stale (§8). |

**Targets & artifacts** (sizes from the qlf-only npm build of 10.1.10; ours
lands within a few %):

| Target | Files | Size | Use |
|---|---|---|---|
| `swipl-web` | `.js` 188K + `.wasm` 2.1M + `.data` 1.6M | ~3.9M | **production** (3 files, cacheable) |
| `swipl-bundle` | single `.js` | 6.0M | one request, trivial embed (base64 overhead) |
| `swipl-bundle-no-data` | single `.js` | 2.5M | code only, supply data yourself |
| `swipl` (node) | `.js`/`.wasm` | — | dev/test + build bootstrap only |

Serve `.wasm` as `application/wasm`, gzip/brotli + long-cache. **Single-threaded
build ⇒ no COOP/COEP headers, no SharedArrayBuffer.**

## 3. Distribution — self-build, not npm

The `swipl-wasm` npm package (latest `v8.0.3`, 2026-06-26) pins the bare
`V10.1.10` tag — GitHub's compare API confirms it is **exactly the 17 commits
behind** our `aa6289399`. It ships a prebuilt `dist/` and cannot be retargeted,
so it can never match our commit. Self-build is therefore required — and easy.
Two routes:

- **(a) Fork their build harness** *(recommended)* — clone `npm-swipl-wasm`,
  set `build-config.json → swipl.commit = aa6289399`, `npm run build` (Docker).
  Reuses their exact toolchain + zlib/pcre2 recipe and emits a drop-in
  `swipl-wasm`-shaped `dist/` (TS types, ESM/CJS wrapper, all four targets).
- **(b) Build the targets directly** from our tree (§2). Leaner, no Docker; you
  wire the JS loader yourself.

The npm package builds SWI *from source in CI via Emscripten* (it does not
vendor a foreign binary), so its `docker/Dockerfile` is a ready reference for
either route. Its CDN (GitHub Pages / jsDelivr) only serves released versions —
moot for us. Note the npm build's own flags (`USE_GMP=OFF`, `STACK_SIZE=1M`,
package subset) match §2, so behaviour parity holds.

## 4. Dependency surface — fully covered

The solver's libraries are `http/json`, `apply`, `aggregate`, `ordsets`,
`lists`, `assoc`, `pairs`, `random`, plus **tabling** (`answer_letters/2`,
`pair_crossings/3`). The Emscripten package set is
`clpqr plunit chr clib http semweb pcre utf8proc`; core autoload libs are always
present; `library(http/json)` resolves via a 2025 compat shim to the standalone
`json` package. **Tabling is a core engine feature, not a package — available.**
Nothing the solver needs is on the "missing under WASM" list (threads, socket,
ssl, crypto, process, …). The CLI layer (`optparse`, `argv`, `halt/1`) is simply
not used in the browser.

## 5. Code porting — one small seam

The engine is plain Prolog. File I/O is isolated to two predicates in
`core.pl` / the CLI: `load_clues/2` (JSON/PL in via `open`+`json_read_dict`) and
`with_output/2` (results out). Replace both with in-memory I/O. Crucially, use a
**deliberate asymmetry** (validated in the spike):

- **Input → query binding.** Pass the JS payload as `Prolog.query(Goal, {In:…})`;
  the JS object becomes a Prolog dict. `doc_to_words/2` already maps that dict to
  the internal `Words` list — reuse it verbatim. Clean and lossless.
- **Output → JSON, not a bound dict.** The layout schema is full of JSON `null`
  (every empty cell) and `across`/`down` values. Per the WASM **reverse**-
  translation table, Prolog's `null`/`true`/`false` atoms come back as the JS
  *strings* `"null"`/`"true"`/`"false"` — so a bound dict would corrupt every
  empty cell. Instead serialise with `json_write_dict` and `JSON.parse` on the
  JS side: this preserves null/bool fidelity **and** is byte-identical to the
  golden CLI output. Return the JSON as a Prolog **atom** (→ a clean JS string;
  a Prolog string would arrive wrapped as a `Prolog.String` instance).

The spike's `solve_browser.pl` implements both (`solve_browser/2` dict→dict for
illustration; `solve_browser_json/2` dict→JSON-atom as the one actually used).
Production should promote this to a proper **exported** `solve_browser/2` in a
small `browser.pl` module and reuse the CLI's size/mode resolvers. This change is
**additive** — the CLI/file paths stay for native use.

## 6. Runtime architecture

- **Package** the app as one `crosswordsmith.qlf` via
  `qcompile('…solve_browser.pl', [include(user)])`, loaded with
  `Prolog.consult(url)`. Validated natively: a 46 KB qlf folds all seven modules
  and runs standalone (loading raw `.pl` from URLs is explicitly slow).
- **Off the main thread: Web Worker.** SWI's blessed design is main-thread
  `forEach`+`heartbeat` yielding, but a genuinely CPU-bound multi-second DFS
  still steals main-thread time between paints. The Worker gives a jank-free UI
  and a trivial hard "Stop". Inside it, run `forEach(goal, in, onAnswer,
  {engine:true, heartbeat:N})`.
- **One instance, reused.** There is no `destroy()`; instantiate `SWIPL()` once.
  `{engine:true}` gives each solve a throwaway engine so per-run state doesn't
  accumulate. Set `stack_limit` generously at init (Prolog stacks are heap-backed
  and grow with `ALLOW_MEMORY_GROWTH`, separate from the 1 MB C stack). Keep a
  warm spare worker so respawn-after-Stop is instant.
- **Cancellation:** `worker.terminate()` (hard) or `Query.close()` at a yield
  point (cooperative). `Prolog.Promise.abort()` only cancels goals blocked on a
  JS promise — **useless for a CPU DFS**; don't rely on it.

## 7. Risks & the validation gate

Do these before trusting the browser build (or the ratchet as a WASM proxy):

1. **Auto-yield starvation — top risk.** Yield fires at the *exit port*, so a
   tight DFS calling only built-ins can run to completion without yielding. The
   `heartbeat` knob (inference-counted, default 10k) is the primary mitigation
   and composes with arrange's inference-based budget. **Verify** the arrange/
   fill search actually yields under a low heartbeat; tune the spike's `50000`.
2. **Inference-count parity is not certified** across native↔WASM. Run
   `benchmarks/workloads.pl` once under `node src/swipl.js` and diff against
   `benchmarks/baseline.json`.
3. **Stack.** 1 MB C stack; Prolog recursion is heap-backed so usually fine, but
   validate the deepest ladder rung. Raise `STACK_SIZE` only if it overflows.
4. **Memory.** Budget conservatively on mobile (~300 MB observed, unconfirmed as
   a hard limit); reuse the instance to avoid cross-run leak accumulation.
5. **qlf under WASM.** Confirm `Prolog.consult(qlf)` folds the modules under WASM
   as it does natively (it should).
6. **Output fidelity.** Confirm the `null`/bool-atom → JS-string behaviour from
   §5 against the real build; it's the reason we chose the JSON output path.

## 8. Corrections owed to `docs/research/swi-vm-wasm-performance.md`

Applied in that file alongside this plan:

- ❌ *"Our pinned 10.1.10 predates that fix — build from a post-fix checkout"* →
  **the pin includes it** (2026-07-01 > 2025-09-29; `merge-base` confirms).
  Building from our own tree is the fast path automatically.
- ❌ *"`-DVMI_FUNCTIONS=ON` (>10% faster)"* → the tree **forces it OFF for
  Emscripten** (`DEFAULT_VMI_FUNCTIONS OFF` under `if(EMSCRIPTEN)`); drop it.
- ➕ The "recommended build flags" are already the in-tree defaults; note
  `WASM_EXCEPTIONS` should stay OFF (alternative, ~30% slower, superseded).
- ➕ Tabling upgraded from "unconfirmed" to **available (core feature)**.

## 9. Phased plan

1. Install emsdk 6.0.1 + stage zlib/pcre2 into `$WASM_HOME` (§2 "Dependency
   staging"); `ninja swipl-web` from our tree; smoke `node src/swipl.js`.
2. Validation gate #2: ladder under node, diff inference counts.
3. Promote `solve_browser.pl` to an exported `browser.pl`; wire the CLI size/mode
   resolvers; qcompile to `crosswordsmith.qlf`.
4. Wire the Worker harness (spike): one solve round-trip, tune `heartbeat`,
   `terminate()` cancel; run gates #1, #5, #6.
5. Only then: the real UI.

---

## Sources

- Build doc (stale, use with care): https://www.swi-prolog.org/build/WebAssembly.md
- WASM manual (version-matched, local): `docs/reference/swi-manual/wasm-{loading,calling,js-call,version}.md`;
  online https://www.swi-prolog.org/pldoc/man?section=wasm
- setjmp perf fix / WASM_EXCEPTIONS tradeoff: https://swi-prolog.discourse.group/t/wasm-performance/9320
- npm-swipl-wasm: https://github.com/SWI-Prolog/npm-swipl-wasm (build-config.json, docker/Dockerfile) ·
  https://www.npmjs.com/package/swipl-wasm · leak issue https://github.com/SWI-Prolog/npm-swipl-wasm/issues/23
- emsdk v4 CI break: https://swi-prolog.discourse.group/t/help-wanted-swipl-wasm-build-not-compatible-with-emsdk-v4/8735
- Blocking/yield vs Worker: https://swi-prolog.discourse.group/t/blocking-operations-in-the-wasm-version-javascript-expertise-needed/5655
- In-browser overview + qlf loading: https://swi-prolog.discourse.group/t/swi-prolog-in-the-browser-using-wasm/5650 ·
  chat80 demo https://wasm.swi-prolog.org/wasm/chat80.html
- Memory/leak guidance: https://swi-prolog.discourse.group/t/firebase-functions-swi-prolog-wasm-node-memory-build-up-issue/6435
- In-tree grounding: `~/src/swipl-devel/{CMakeLists.txt,cmake/port/Emscripten.cmake,cmake/EmscriptenTargets.cmake,cmake/PackageSelection.cmake,scripts/configure,src/wasm/README.md}`
