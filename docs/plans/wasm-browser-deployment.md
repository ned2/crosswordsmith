# Plan: crosswordsmith in the browser (SWI-Prolog WASM)

Status: **strategy / pre-implementation** · Drafted 2026-07-05.

Goal: run the crosswordsmith solver in the browser under SWI-Prolog's WASM
runtime, at the **same SWI version we build natively**, without regressing the
CLI. This doc is the strategy; the runnable client (promoted from a spike, and
validated end-to-end in headless Chrome) lives under
[`wasm/`](../../wasm/) — client, `build/build-wasm.sh`, and headless tests.

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
**JSON in and JSON out**; reuse one instance; cancel by `worker.terminate()`.
**The spike works end-to-end in headless Chrome** (renders a solved grid); the
full toolchain, build, and browser round-trip are validated (2026-07-05, §6/§9),
and **both validation gates are certified**: gate #2 inference parity (§7) and
gate #1 large-search responsiveness + cancel (§7 — the Worker isolates the UI, and
`terminate()` cancels; the heartbeat turned out not to be load-bearing).

**The one non-obvious subtlety is build *provenance*, not the app.** WASM is
32-bit, so anything the build "compiles" — the `boot.prc`/library `.qlf` during
bootstrap, and our app `.qlf` — must be produced by a *same-pointer-size* Prolog.
That means: no native "friend" build (let the wasm binary self-bootstrap under
node), and produce `crosswordsmith.qlf` with the wasm build's `node src/swipl.js`,
not native `swipl`. This is the one place the "version parity is free" framing
can mislead — parity of *source* is free; word-size compatibility is a separate
constraint. (Adversarial review, 2026-07-05, caught this; §2/§6.)

## 2. Build — exact-version parity from our tree

`$WASM_HOME` is exported to `$HOME/wasm` (the staging prefix, see below). **As of
2026-07-05 this is all done:** emsdk 6.0.1 is installed at `$WASM_HOME/emsdk`,
zlib/pcre2 are staged, and `swipl-web` is built + node-smoke-tested at
`~/src/swipl-devel/build.wasm` (see the validation note under the recipe). The
recipe below stays as the reproducible record.

```bash
# 0. toolchain — emsdk 6.0.1 (the version npm-swipl-wasm pins in build-config.json;
#    swipl-devel itself has NO wasm CI job, so this is not attested in-tree — smoke-
#    test it. ≥4.0.15 works post-fix; 3.1.37 in the stale official doc is obsolete;
#    emsdk 4.0.0 briefly broke swipl-wasm CI.)
git clone https://github.com/emscripten-core/emsdk $WASM_HOME/emsdk
cd $WASM_HOME/emsdk && ./emsdk install 6.0.1 && ./emsdk activate 6.0.1 && source ./emsdk_env.sh

# 1. stage wasm-side deps into $WASM_HOME — see "Dependency staging" below

# 2. pin the source to our exact commit, then configure via the in-tree profile
cd ~/src/swipl-devel && git checkout aa6289399      # detached — else a later pull drifts HEAD
mkdir -p build.wasm && cd build.wasm
../scripts/configure -y wasm # -y MUST be the first arg (else configure hangs on an
                             # interactive "Run Cmake?" prompt with no tty). The profile
                             # sets toolchain, CMAKE_FIND_ROOT_PATH=$WASM_HOME, USE_GMP=OFF,
                             # INSTALL_PROLOG_SRC=OFF, Release, and NO native friend.
ninja swipl-web              # or: ninja swipl-web swipl-bundle
```

> **Validated end-to-end on 2026-07-05** with emsdk 6.0.1 + cmake 4.2.3: builds
> clean, `node src/swipl.js --version` → `10.1.10 for wasm-emscripten`,
> `library(http/json)` loads and `json_write_dict` round-trips JSON `null`. Two
> operational gotchas surfaced (folded into the recipe above):
> - **`configure` needs `-y` first** — without a tty it otherwise loops forever
>   on its confirm prompt.
> - **cmake 4.2.3 is in emscripten's unsupported-shared-lib window** (it warns
>   "does not support emscripten shared libraries; use cmake <4.2.0 or >4.3.3" on
>   every compile probe). **Harmless here** — the wasm kernel links statically, so
>   `swipl-web` builds fine — but the warnings are loud; don't be alarmed. Only a
>   concern if you later need actual emscripten shared modules.

Use the in-tree `scripts/configure wasm` profile — do **not** hand-roll the
`emcmake cmake …` line, and do **not** pass `-DSWIPL_NATIVE_FRIEND`. The build is
two-stage: a Prolog of the **same pointer size** must generate `boot.prc` + the
library `.qlf`, and the profile's no-friend path runs the freshly-built wasm
binary under the node emulator to do exactly that (`cmake/Ports.cmake`,
`cmake/EmscriptenTargets.cmake`). Our native `~/src/swipl-devel/build` is x86-64
(8-byte pointers) — pointing a friend at it would make a 64-bit swipl emit the
wasm32 target's boot/qlf: a silent word-size mismatch. The node bootstrap is
slower but correct. (Prefer all this over the **official `WebAssembly.md`, which
is ~3 years stale** — Emscripten 3.1.37, a non-`emcmake` line.)

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
| `USE_GMP` | **OFF** (profile sets it) | LibBF bignums; the WASM norm. arrange/fill are integer-only, so the weaker LibBF rational arithmetic is irrelevant. |
| `INSTALL_QLF` / `INSTALL_PROLOG_SRC` | **ON / OFF** | qlf-only image → smaller `.data`, no source tracing. `INSTALL_QLF=ON` is the global default (don't pass it); the profile sets `INSTALL_PROLOG_SRC=OFF`. |
| `WASM_EXCEPTIONS` | **OFF (default)** | *Alternative* to the setjmp fix, ~30% slower; superseded by the fix already in our pin (§8). Not additive. |
| `SWIPL_NATIVE_FRIEND` | **do NOT set** | A friend must match the *target* pointer size (wasm32 = 4 bytes); our native build is x86-64 = 8 bytes, so a friend would emit a mismatched `boot.prc`/`.qlf`. Leave unset → the profile bootstraps boot/qlf by running the wasm binary under node (`cmake/Ports.cmake` "compatible ⇒ same pointer size"). |
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
`pair_crossings/3`). The Emscripten package **set** is
`clpqr plunit chr clib http semweb pcre utf8proc`; core autoload libs are always
present. **`json` is not in that literal list, but it is still built:**
`json` is now its own package (`packages/json/`, with the 2025
`library(http/json)` compat shim reexporting `library(json)`), and
`SWIPL_PKG_DEPS_http` includes `json` under `EMSCRIPTEN`
(`cmake/PackageSelection.cmake`). The build's dependency-resolution fixpoint
(`check_package_dependencies`, run in a `while(new)` loop) pulls any missing dep
into the built set, so `http` drags `json` in — and npm-swipl-wasm's own package
list ships `json` for the same reason. **Tabling is a core engine feature, not a
package — available.** Nothing the solver needs is on the "missing under WASM"
list (threads, socket, ssl, crypto, process, …). The CLI layer (`optparse`,
`argv`, `halt/1`) is simply not used in the browser.
> **Confirmed on the 2026-07-05 build** — with one wasm nuance the browser run
> exposed. Configure logged `-- Added packages json: Required by package http`, so
> `json` **is** built. But the two runtimes differ on the *file alias*:
> - Under `node src/swipl.js` (reads the full library tree from disk),
>   `use_module(library(http/json))` loads and `json_write_dict` round-trips
>   `{"grid":null,…}` with `null` preserved.
> - In the **web image** (`swipl-web.data`), `library(json)` resolves but the
>   compat-shim *alias* `library(http/json)` does **not** — so every app module's
>   explicit `:- use_module(library(http/json))` **fails with a warning** at qlf
>   load. It's non-fatal: `json_write_dict`/`json_read_dict` **autoload from
>   `library(json)`**, so output still works (verified in-browser — the grid
>   renders). This is the kernel of truth in the review's CRITICAL #2: the
>   *package* is present and the *predicates* work, but the `http/json` alias has
>   a web-image resolution gap. Production could silence the noise by importing
>   `library(json)` directly, or registering the `ext/json` dir on the library
>   search path before load.

## 5. Code porting — one small seam

The engine is plain Prolog. File I/O is isolated to two predicates in
`core.pl` / the CLI: `load_clues/2` (JSON/PL in via `open`+`json_read_dict`) and
`with_output/2` (results out). Replace both with in-memory I/O. **JSON in, JSON
out** — the symmetric choice, *corrected from an earlier "input via query
binding" plan after the browser run disproved it* (see below):

- **Input → JSON string, parsed in Prolog.** The tempting path is a query binding
  (`Prolog.query(Goal, {In:…})`, JS object → Prolog dict). It works on the main
  thread for simple cases but has a **silent trap the browser exposed**: the
  JS→Prolog binding delivers JS strings as Prolog **atoms**, so a clue's
  `answer` fails `doc_to_words/2`'s `string(RawAnswer)` check
  (`every entry needs a string "answer"`). Instead pass `JSON.stringify(input)`
  and `json_read_dict/3` it in Prolog: every value lands with the right type
  (string/number/null/bool). This is the input-side twin of the output quirk.
- **Output → JSON, not a bound dict.** The layout schema is full of JSON `null`
  (every empty cell) and `across`/`down` values. Per the WASM **reverse**-
  translation table, Prolog's `null`/`true`/`false` atoms come back as the JS
  *strings* `"null"`/`"true"`/`"false"` — so a bound dict would corrupt every
  empty cell. Instead serialise with `json_write_dict` and `JSON.parse` on the
  JS side: this preserves null/bool fidelity **and** is byte-identical to the
  golden CLI output. Return the JSON as a Prolog **atom** (→ a clean JS string;
  a Prolog string would arrive wrapped as a `Prolog.String` instance).

The spike's `solve_browser.pl` implements the tolerant path: `solve_browser_str/2`
(JSON string → JSON string, the one the Worker calls) wraps `solve_browser_json/2`
(dict → JSON). Production should promote this to a proper **exported**
`solve_browser/2` in a small `browser.pl` module and reuse the CLI's size/mode
resolvers. This change is **additive** — the CLI/file paths stay for native use.

## 6. Runtime architecture

- **Package** the app as one `crosswordsmith.qlf` via
  `qcompile('…solve_browser.pl', [include(user)])`, loaded with
  `Prolog.consult(url)` (loading raw `.pl` from URLs is explicitly slow). A 46 KB
  qlf folds all seven modules and runs standalone — but that was validated in
  *native x86* swipl, which only proves 64→64. **QLF is word-size-specific** (it
  encodes word-sized instruction variants and byte offsets, and the load-time
  check is version+VM-signature only — both word-size-independent, so a
  mismatched qlf loads-and-misbehaves rather than failing cleanly). So the app
  qlf **must be produced by a same-pointer-size Prolog** — i.e.
  `node build.wasm/src/swipl.js -g "qcompile('…/solve_browser.pl',[include(user)]),halt"`
  (the wasm build's own node binary, with NODERAWFS reaching the app sources by
  absolute path), **not** native `swipl`.
- **Off the main thread: Web Worker.** SWI's blessed design is main-thread
  `forEach`+`heartbeat` yielding, but a genuinely CPU-bound multi-second DFS
  still steals main-thread time between paints. The Worker gives a jank-free UI
  and a trivial hard "Stop". Inside it, run `forEach(goal, in, onAnswer,
  {engine:true, heartbeat:N})`.
- **One instance, reused.** There is no `destroy()`; instantiate `SWIPL()` once.
  `{engine:true}` gives each solve a throwaway engine so per-run state doesn't
  accumulate. Set `stack_limit` **under the memory budget** (~200–256 MB, *not*
  1 GB): Prolog stacks are heap-backed and grow with `ALLOW_MEMORY_GROWTH`, so an
  over-large limit lets a runaway search grow the wasm heap until the browser
  `abort()`s the tab (uncatchable) instead of throwing a recoverable
  `resource_error(stack)` first. This logical limit is separate from the 1 MB C
  `STACK_SIZE`. Keep a warm spare worker so respawn-after-Stop is instant.
- **Cancellation:** `worker.terminate()` is the hard Stop — a thread kill that is
  **independent of yielding** (it does not need a heartbeat or a cooperative
  check), and **gate #1 measured it prompt** (0 ms to call, search killed, no late
  result). The cooperative `Query.close()` alternative is *not reachable* in this
  single-worker design: gate #1 showed a running solve does **not** service the
  worker's incoming messages at any heartbeat (the yield cooperates among Prolog
  engines, not the JS message queue), so there is no way to *tell* the worker to
  close the query mid-search. `terminate()` is therefore the only practical cancel.
  `Prolog.Promise.abort()` only cancels goals blocked on a JS promise — **useless
  for a CPU DFS**; don't rely on it.

### Worker gotchas — validated end-to-end in headless Chrome (2026-07-05)

The spike **runs**: headless Chrome loads the page, the Worker instantiates
`swipl-web`, consults `crosswordsmith.qlf`, and renders a solved 5×5 grid. Three
Worker-specific fixes were needed — none of which node testing can surface,
because they are all about the *browser Worker* environment:

- **`self.window = self` before `importScripts("./swipl-web.js")`.** SWI's wasm
  URL helpers (e.g. `prolog.url_properties`, used by `Prolog.consult(URL)`) reach
  for a browser `window`, which does **not** exist in a Worker → `ReferenceError:
  window is not defined` and the consult fails. Aliasing the Worker global as
  `window` (its `location` then resolves to the worker URL) fixes it.
- **Absolute URL for `Prolog.consult`.** With no `window`, SWI can't derive a base
  URL, so a relative `"./crosswordsmith.qlf"` 404s. Pass
  `new URL("./crosswordsmith.qlf", self.location.href).href`.
- **JSON in, not a query binding** (see §5) — the JS-string→atom trap.

Non-fatal noise to expect: each app module's `:- use_module(library(http/json))`
prints `source_sink library(http/json) does not exist` **during qlf load**. In
the web image the compat-shim *file alias* `library(http/json)` doesn't resolve
(only `library(json)` does), so the explicit directive fails — but
`json_write_dict`/`json_read_dict` **autoload from `library(json)`**, so output
works anyway (the grid renders). Cosmetic; §4 has the detail. The
`crosswordsmith.qlf` request also shows one `net::ERR_ABORTED` (SWI's consult
probes then fetches); harmless — the qlf loads and the solve completes.

Reproduce headlessly with Playwright + system Chrome — `node wasm/test/headless.mjs`
(see `wasm/README.md`).

## 7. Risks & the validation gate

Do these before trusting the browser build (or the ratchet as a WASM proxy):

1. ✅ **Auto-yield starvation — RESOLVED, and the framing was wrong (2026-07-05).**
   Measured a real ~38.3M-inference search (`ladder_15x15_36w`) in headless Chrome
   via `wasm/test/yield_probe.mjs`, sweeping `heartbeat` 10k→2e9. Findings:
   - **The search completes correctly in-browser** (36 words placed, ~5 s, 256 MB
     stack held) — the yield path scales far past the 5-word toy.
   - **The page's main thread stays perfectly responsive throughout** — main-thread
     `setInterval` drift `0 ms` at *every* heartbeat, *including the no-yield
     control*. Responsiveness comes from running the search in a **Web Worker (off
     the UI thread)**, **not** from the heartbeat. The heartbeat is *not*
     load-bearing here; the no-yield control even ran ~3 % faster.
   - **The heartbeat does not let the worker service its own messages mid-search**
     (0 mid-search pongs at every heartbeat, 10k included): the yield cooperates
     among *Prolog engines*, it does not drain the *worker's JS message queue*. So
     you cannot gracefully cancel a running solve by posting it a message.
   - **Cancellation is `worker.terminate()`** — measured prompt (0 ms to call, the
     search is killed, no late result), which is exactly what `main.js` already
     does. `Prolog.abort()` would only interrupt `await/2` points, not a pure DFS.

   Net: the "top risk" is fully mitigated by the **Worker + `terminate()`**
   architecture already in place. `heartbeat: 50000` is kept as a mild, harmless
   default (a few % throughput tax, no UX cost) and left caller-overridable.
2. ✅ **Inference-count parity — CERTIFIED (2026-07-05).** Ran the full arrange
   ladder (all 12 rungs, 9×9/15×15/21×21, ~28k–38.5M inferences) under both native
   `swipl` and the wasm `node src/swipl.js`, measuring the same search-layer count
   the ratchet records (`call_time` around `arrange_best_layout/6`, no `once/1`).
   **All 12 rungs byte-identical native↔wasm**; native reproduced
   `benchmarks/baseline.json` to within 0.0035% (one rung off by a single
   inference — a first-warm-iteration wobble, and wasm shows the *same* +1, so it's
   an engine artifact, not a VM difference). Wall was ~2.7× slower under wasm on
   the heavy set — i.e. *only the wall-per-inference constant changed*, exactly the
   ratchet's premise. Reproduce: `wasm/test/inference_parity.pl` (see wasm/README).
3. **Stack.** 1 MB C stack; Prolog recursion is heap-backed so usually fine, but
   validate the deepest ladder rung. Raise `STACK_SIZE` only if it overflows.
4. **Memory.** Budget conservatively on mobile (~300 MB observed, unconfirmed as
   a hard limit); set `stack_limit` under it (§6); reuse the instance to avoid
   cross-run leak accumulation.
5. ✅ **qlf must be wasm-produced — DONE.** Built `crosswordsmith.qlf` with
   `node build.wasm/src/swipl.js`; it loads via `Prolog.consult()` under WASM and
   is **self-contained** (ran the selftest from a dir with no source reachable,
   and rendered a grid in-browser). A native-x86 qlf would not have (§6).
6. ✅ **Output fidelity — DONE (in browser).** `json_write_dict` → `JSON.parse`
   preserved JSON `null` for empty cells; the headless-Chrome grid rendered
   correctly. Input fidelity too: JSON-in (§5) fixed the JS-string→atom trap.

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

1. ✅ **DONE 2026-07-05.** Installed emsdk 6.0.1 + staged zlib 1.3.2 / pcre2 10.47
   into `$WASM_HOME`; `ninja swipl-web` from our tree (three artifacts +
   `src/swipl.js` node bootstrap); smoked under `node src/swipl.js` — version
   `10.1.10 for wasm-emscripten`, `http/json` + `json_write_dict` verified. Build
   dir `~/src/swipl-devel/build.wasm` (gitignored; source tree left at the pin).
2. ✅ **DONE 2026-07-05.** Validation gate #2: ran the arrange ladder under both
   VMs — **all 12 rungs byte-identical native↔wasm** (~28k–38.5M inferences),
   native within 0.0035% of `baseline.json`. Inference-as-wasm-proxy certified;
   harness kept as `wasm/test/inference_parity.pl`.
3. ✅ **DONE 2026-07-05 (spike level).** `qcompile`d the app to
   `crosswordsmith.qlf` **with `node build.wasm/src/swipl.js`**; verified
   self-contained (loads with no source reachable). *Still to productionise:*
   promote `solve_browser.pl` → an exported `browser.pl`, reuse the CLI size/mode
   resolvers (the spike takes size/mode straight from the payload). **Scoped
   predicate-by-predicate in §9.1.**
4. ✅ **DONE 2026-07-05.** Worker harness runs end-to-end in **headless Chrome** —
   one solve round-trip renders a grid; gates #5, #6 passed. Fixes folded in:
   `window` shim, absolute-URL consult, JSON-in (§6). **Gate #1 also cleared**
   (`wasm/test/yield_probe.mjs`): a real ~38.3M-inference search completes
   in-browser with the UI fully responsive (main-thread drift 0 ms — the Worker
   boundary, not the heartbeat, is what isolates the UI), and `worker.terminate()`
   is a prompt cancel. Heartbeat found *not* load-bearing; kept at 50000 (see §7).
5. Only then: the real UI.

### 9.1 Sub-plan — promote the browser entry to an exported `browser.pl` module

*Stress-tested 2026-07-05: a read-only red-team pass re-verified every claim against
source and surfaced three correctness gaps (per-request state reset, error channel,
in-process golden capture) plus missing build/loader wiring — all folded in below.*

Scopes phase-3's productionization tail (item 3). Two concerns in one deliverable:
**(A) API hygiene** — replace the spike's `Module:Pred` reaches with a declared
public surface; **(B) resolver reuse** — drive the browser through the *same*
option-resolution as the CLI so the two can't drift. Nothing here fixes a live bug:
the worker's JSON path already flows through the **exported** `arrange_solve/4` and
is byte-identical to golden CLI output for the features it exposes (fixed/max size,
strict/best-effort). Parity + hygiene, not a fix.

**Architecture (corrected by the review).** The *logic* moves into a new **project
module `prolog/crosswordsmith/browser.pl`** (`:- module(crosswordsmith_browser,
[browser_solve_json/2])`), loaded via `load.pl` and folded into the qlf by
`include(user)`. The client file **`wasm/client/solve_browser.pl` STAYS** as the
thin worker entry: the worker calls `solve_browser_str/2` **unqualified**
(worker.js:88) and `build-wasm.sh` qcompiles `solve_browser.pl` by name into
`crosswordsmith.qlf` (build-wasm.sh:108/110) — so that file keeps
`solve_browser_str/2`/`_json/2` in `user`, now as **thin adapters delegating to
`crosswordsmith_browser:browser_solve_json/2`**. Worker and build stay untouched.
(Renaming the client entry instead would force edits to worker.js, the build script,
and ~6 README refs — avoid.)

**‼️ Mandatory new invariant — per-request state reset (top gap the review found).**
Seed and check-target live in **instance-global mutable state** the reused WASM
instance carries across solves, and `{engine:true}` does **not** clear it:
`set_search_seed/1` reseeds the global RNG (`set_random(seed(N))`, core.pl) and
asserts `search_seed/1`; `set_check_target/1` asserts `check_target_override/1`
(arrange.pl). The CLI is immune only because it's a fresh process **and**
re-establishes a clean baseline every run (`set_check_target(-1)`,
`apply_seed_mode(none)→set_search_seed(-1)`). So `browser_solve_json/2` **must reset
unconditionally at the top of every request** — `set_search_seed(-1)`,
`set_check_target(-1)`, `set_verbose(false)` — before applying request keys, or one
user's `seed:42`/`checkTarget:N` leaks into the next seedless solve. The spike has no
such step; here it is required.

**Error channel — resolvers can't be lifted verbatim.** The CLI resolvers report
errors via `cli_error/1`, which is **script-local in `crosswordsmith`** and does
`format(user_error,…), fail` — it neither throws nor halts. So (a) a resolver clause
moved into a module still calls `cli_error` → **undefined predicate** at runtime;
(b) even if reachable, a *failing* resolver yields no `forEach` binding, which the
worker reports as the generic `"no layout (search failed)"` (worker.js:93) with the
real reason lost to a wasm console. The lift must therefore **refactor error branches
to `throw(error(browser_…, _))`** (caught by the worker's `try/catch`, worker.js:98,
→ an error-JSON reply), with the CLI adapter catching-and-`cli_error`-ing to preserve
today's CLI UX. This refactor is a **prerequisite of Part B**, not an afterthought.

**Part A — the module + public API.**

| predicate | file:line | today | action |
|---|---|---|---|
| `arrange_solve/4` | arrange.pl:21 | ✅ exported | reuse as-is (output path already does) |
| `doc_to_words/2` | core.pl:204 | ❌ internal | **export from core.** It's the JSON-only post-parse mapping (schema→words), *not* the full `load_clues/2` (which also dispatches `.pl` terms) — correct for a JSON-only browser. |
| `arrange_best_layout/5`, `arrange_diag_layout_dict/5` | arrange.pl | ❌ internal | only the *illustrative* dict path (`solve_browser/2`) uses them — **drop that path** (JSON-in/out is all the browser needs). Nothing else depends on them (arrange uses them internally at 290/392/890; `inference_parity.pl` uses the `/6` arity). |

Dropping the dict path means **`browser_selftest` must be updated** — it calls
`solve_browser/2` first (solve_browser.pl:114), so it is *not* "re-run unchanged".

**Part B — lift the CLI resolvers into a shared module.** All live **script-local in
`crosswordsmith`**; "reuse" = extract the logic into a module (e.g.
`crosswordsmith_request`, or `core`/`arrange`), leaving the CLI a thin `Opts`→values
adapter. Every entry below is **logic-pure but errors via `cli_error`** — i.e. each
needs the error-channel refactor above before it can be shared:

| resolver | `crosswordsmith` line | nature | split action |
|---|---|---|---|
| `resolve_size_flags/4` | 234–237 | logic-pure; err→cli_error | move + export; browser calls directly |
| `validate_size_flag/2`, `resolve_size/2` | 226–228, 241–243 | logic-pure; err→cli_error | move + export (positive-int / default-15 guards) |
| `drop_contract/2` | 165–171 | reads `Opts` | split: pure `drop_of(Strict,Best,Drop)` + CLI adapter |
| `seed_mode/3` | 194–197 | logic-pure; err→cli_error | move + export |
| `validate_seed/1`, `guard_seed_combos/4` | 187–189, 203–210 | logic-pure; err→cli_error | move + export |
| `apply_seed_mode/1` | 214–219 | state-setting **+ side effects** | delegates to core's exported `set_search_seed/1`·`set_shuffle_seed/1`, but its `shuffle` clause also `verbose_report`s to `user_error` and the fresh seed needs surfacing to JS, not stderr — prefer calling the two core setters directly and returning the seed in the result over moving the wrapper |
| `set_check_target/1`, `validate_check_target/1` | arrange (exported) / 178–180 | mixed | `set_check_target/1` exported already; move `validate_check_target/1` |
| `load_fragment/3` | arrange.pl:538 (exported) | **file I/O** | pure core genuinely is `fragment_dict_words/3` (open+read wraps it); add in-memory `load_fragment_dict/3` (JSON text→dict→`fragment_dict_words/3`) — the exact `doc_to_words`-vs-`load_clues` split. Export `fragment_dict_words/3` if needed. |
| `reconcile_fragment_size/3` | arrange (exported) | pure (throws on mismatch) | reuse as-is |
| `arrange_action/8` | 247–268 | **not pure** — file I/O (`with_output`, `load_fragment`) + `cli_error` | mirror as `browser_action/…`: swap `with_output`→`with_output_to(atom(...))`, `load_fragment`→`load_fragment_dict`, `cli_error`→throw; routes plain/fragment/candidates/enumerate |

**Part C — input schema growth.** `browser_solve_json/2`'s dict gains optional keys
mapping 1:1 to CLI flags, each defaulting to today's behaviour when absent (**after**
the mandatory reset above): `maxSize` (↔ `--max-size`, exclusive with `size`),
`fragment` (a fragment **dict**, not a path), `seed`/`shuffle`, `checkTarget`,
`candidates`, `enumerate`. Absent ⇒ current fixed-size / strict path. (Fragment
sub-dicts arrive already-parsed — watch `default_tag` consistency vs the CLI's
file-read path; `fragment_dict_words/3` is tag-agnostic, but assert it in Part D.)

**Part D — verification (parity gate).** Prove `browser_solve_json/2` output === the
golden arrange bytes for the *same* request. **Capture golden IN-PROCESS**, not by
shelling to the CLI: `wasm/test/inference_parity.pl` deliberately avoids
`library(process)`, `with_output/2` has **no `-`→stdout case** (`--out -` writes a
file literally named `-`), and loading the `crosswordsmith` script fires its
`:- initialization(main)`→`halt`. So the harness (load via `load.pl`, no `process`)
diffs the browser JSON atom against a sibling `with_output_to(string(Golden),
arrange_solve(Words,Grid,Drop,Mode))` (+ the fragment/candidates/enumerate variants).
Matrix: {fixed, max} × {strict, best-effort} × {no-seed, seed N, fragment} on a
couple of ladder fixtures — **skip the `best-effort × seed` cells**
(`guard_seed_combos/4` rejects them by design; assert they error instead).
Byte-identical across the valid matrix = done. Also re-run the (updated)
`browser_selftest` + `wasm/test/headless.mjs`.

**Build & loader wiring (don't forget).** Adding the module means **`load.pl`** gains
`:- use_module(crosswordsmith(browser)).` (after arrange, load.pl:26+) so
`browser_solve_json/2` is defined at qcompile time and folded into the qlf.
**`build-wasm.sh` and `worker.js` stay as-is** *because* the client entry keeps its
`solve_browser.pl` name and `solve_browser_str/2`/`_json/2` in `user`. A later rename
of the client entry would require updating build-wasm.sh:108/110, worker.js:88, and
the `wasm/README.md` refs together.

**Cut lines (incremental landing).** *Minimal* = Part A only, **keeping** the local
size/mode/drop helpers (moved into `browser.pl`; the "delete, superseded by B" step
waits for B): a new module + one `core` export + one `load.pl` line — **low
conflict** (no shared *logic*, and the `crosswordsmith` script untouched), so it can
land on the first post-rebase pass. *Full* = A+B+C+D: CLI feature parity, incl. the
error-channel refactor and per-request reset. **Only Part B is genuinely
conflict-heavy** (surgery on the `crosswordsmith` script + `cli_error`); schedule Full
behind the real UI (item 5), which is what actually needs fragments/seeds.

*(This whole sub-plan is doc-only regardless: the branch is parked by decision until
the fill campaign quiets and this branch rebases — see §1/the branch status. The
Minimal/Full split above governs the order of work* once *it unparks.)*

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
