# WASM browser spike

Throwaway spike for the browser milestone — see
[`docs/plans/wasm-browser-deployment.md`](../../docs/plans/wasm-browser-deployment.md)
for the full strategy. This proves the end-to-end shape: a file-free `arrange`
entry point, called from a Web Worker, rendering a grid in the page.

## What's here

| File | Role |
|---|---|
| `solve_browser.pl` | file-free entry points: `solve_browser/2` (dict→dict) and `solve_browser_json/2` (dict→JSON atom, the one the worker uses). Loads the project via `load.pl`. |
| `worker.js` | the SWI-Prolog WASM instance, off the main thread. Input via query binding, output via `JSON.parse` of a JSON atom. |
| `main.js` | page controller; "Stop" = `worker.terminate()` + warm-spare respawn. |
| `index.html` | minimal UI (words in, grid out). |

The Prolog side is testable **today** under native swipl (no build needed):

```bash
swipl -q -g browser_selftest -t halt spikes/wasm-browser/solve_browser.pl
# dict path OK: gridLength=5, 5 placed words
# json path OK: 1706 chars
```

## Building the WASM artifacts (needed for the browser half)

Prereq: **emsdk 6.0.1** (SWI's current CI pin) — not yet installed on this box.

```bash
# 1. toolchain (WASM_HOME=$HOME/wasm is exported as the staging prefix)
git clone https://github.com/emscripten-core/emsdk $WASM_HOME/emsdk
cd $WASM_HOME/emsdk && ./emsdk install 6.0.1 && ./emsdk activate 6.0.1 && source ./emsdk_env.sh
# stage wasm-side deps (zlib mandatory, pcre2) into $WASM_HOME — exact commands
# in docs/plans/wasm-browser-deployment.md § "Dependency staging".

# 2. build swipl-web from OUR pinned source tree (exact-version parity)
cd ~/src/swipl-devel && mkdir -p build.wasm && cd build.wasm
emcmake cmake -DCMAKE_TOOLCHAIN_FILE=$EMSDK/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake \
  -DCMAKE_FIND_ROOT_PATH=$WASM_HOME -DCMAKE_BUILD_TYPE=Release \
  -DUSE_GMP=OFF -DINSTALL_QLF=ON -DINSTALL_PROLOG_SRC=OFF -DINSTALL_DOCUMENTATION=OFF \
  -DSWIPL_NATIVE_FRIEND=$HOME/src/swipl-devel/build -G Ninja ..
ninja swipl-web

# 3. drop the three artifacts next to these files
cp src/swipl-web.js src/swipl-web.wasm src/swipl-web.data <this-dir>/
```

## Building the app `.qlf`

Fold the app into one Quick-Load File served alongside the artifacts:

```bash
swipl -g "qcompile('spikes/wasm-browser/solve_browser.pl', [include(user)]), halt" -t 'halt(1)'
cp spikes/wasm-browser/solve_browser.qlf spikes/wasm-browser/crosswordsmith.qlf
```

> ✅ **Validated natively:** `include(user)` *does* fold all seven modules —
> a 46 KB `solve_browser.qlf` loaded on its own into a fresh native swipl runs
> `browser_selftest` cleanly (both paths). Residual gate: confirm
> `Prolog.consult('./crosswordsmith.qlf')` behaves identically under WASM (it
> should). Fallback if not: serve the source tree and `consult('./load.pl')`.

## Running

Serve this directory over HTTP (the worker + `.wasm`/`.data` need real URLs,
and `.wasm` must be served as `application/wasm`):

```bash
python3 -m http.server 8080   # then open http://localhost:8080/
```

Single-threaded build ⇒ **no COOP/COEP headers required.**

## Known spike shortcuts (not production)

- Reaches module internals (`Module:Pred`); production adds a proper exported
  `solve_browser/2` in a small `browser.pl`.
- Grid size/mode read straight from the input dict; the CLI's `--max-size`
  crop / fragment / seed machinery is out of scope.
- `heartbeat: 50000` in `worker.js` is a guess — tune against a real search
  (the auto-yield-starvation risk in the plan).
