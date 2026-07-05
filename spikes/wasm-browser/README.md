# WASM browser spike

Throwaway spike for the browser milestone — see
[`docs/plans/wasm-browser-deployment.md`](../../docs/plans/wasm-browser-deployment.md)
for the full strategy. This proves the end-to-end shape: a file-free `arrange`
entry point, called from a Web Worker, rendering a grid in the page.

## What's here

> ✅ **Validated end-to-end in headless Chrome (2026-07-05):** loading the page,
> clicking Solve, and rendering a solved 5×5 grid all work — the solver runs in a
> Web Worker off the main thread. See "Browser gotchas" below for the three
> Worker-specific fixes that were needed.

| File | Role |
|---|---|
| `solve_browser.pl` | file-free entry points: `solve_browser_str/2` (**JSON string → JSON string**, what the worker calls), `solve_browser_json/2` (dict → JSON atom), `solve_browser/2` (dict→dict, illustrative). Loads the project via `load.pl`. |
| `worker.js` | the SWI-Prolog WASM instance, off the main thread. **JSON in, JSON out.** |
| `main.js` | page controller; "Stop" = `worker.terminate()` + warm-spare respawn. |
| `index.html` | minimal UI (words in, grid out). |
| `probe.html` | throwaway diagnostic page (library resolution / json autoload); handy if the browser path breaks again. |

The Prolog side is testable **today** under native swipl (no build needed):

```bash
swipl -q -g browser_selftest -t halt spikes/wasm-browser/solve_browser.pl
# dict path OK: gridLength=5, 5 placed words
# json path OK: 1706 chars
# str  path OK: 1706 chars
```

### Browser gotchas (learned the hard way; all fixed in `worker.js`/`solve_browser.pl`)

None of these show up under `node` — they are specific to the browser **Worker**:

1. **`self.window = self`** before `importScripts("./swipl-web.js")`. SWI's wasm
   URL helpers reach for `window` (absent in a Worker) → `window is not defined`
   and `consult` fails.
2. **Absolute URL for `Prolog.consult`** (`new URL("./crosswordsmith.qlf",
   self.location.href).href`) — no `window` ⇒ no base URL ⇒ relative path 404s.
3. **JSON in, not a query binding.** The JS→Prolog binding turns JS strings into
   Prolog **atoms**, so a clue's `answer` fails core's `string` check. The worker
   sends `JSON.stringify(input)`; `solve_browser_str/2` parses it with
   `json_read_dict`. Symmetric with the JSON-out path.

Expected non-fatal noise: `source_sink library(http/json) does not exist` warnings
during qlf load — in the web image the `http/json` *alias* doesn't resolve, but
`json_write_dict`/`json_read_dict` autoload from `library(json)`, so output works.

### Reproduce the headless browser test

Needs a running server (below) + Playwright driving system Chrome (no browser
download):

```bash
mkdir -p /tmp/pwtest && cd /tmp/pwtest && npm init -y >/dev/null
PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 npm i playwright
# then a ~30-line script: chromium.launch({channel:'chrome',headless:true}),
# goto http://127.0.0.1:8080/, click #solve, wait for #status to match /placed/,
# read #grid. Expect: "placed 5 words on 5×5" and the RAT/ARC/CAT/TAR/CAR grid.
```

## Building the WASM artifacts (needed for the browser half)

Prereq: **emsdk 6.0.1** (the version npm-swipl-wasm pins; swipl-devel has no wasm
CI job so it isn't attested in-tree). **Done once on 2026-07-05** — the full
recipe below has been run and validated; it's kept here as the reproducible
record.

```bash
# 1. toolchain (WASM_HOME=$HOME/wasm is exported as the staging prefix)
git clone https://github.com/emscripten-core/emsdk $WASM_HOME/emsdk
cd $WASM_HOME/emsdk && ./emsdk install 6.0.1 && ./emsdk activate 6.0.1 && source ./emsdk_env.sh
# stage wasm-side deps (zlib mandatory, pcre2) into $WASM_HOME — exact commands
# in docs/plans/wasm-browser-deployment.md § "Dependency staging". (zlib 1.3.2 +
# pcre2 10.47 used; both resolved by find_package from $WASM_HOME.)

# 2. build swipl-web from OUR pinned commit via the in-tree profile.
#    Do NOT pass -DSWIPL_NATIVE_FRIEND: the friend must match the TARGET pointer
#    size (wasm32 = 4 bytes); our native build is x86-64 = 8 bytes, so a friend
#    would emit a mismatched boot.prc/.qlf. No friend => the wasm binary self-
#    bootstraps under node (slower, but the only correct path). The profile sets
#    the toolchain, CMAKE_FIND_ROOT_PATH=$WASM_HOME, USE_GMP=OFF, Release, etc.
cd ~/src/swipl-devel && git checkout aa6289399     # detached — pin the exact commit
mkdir -p build.wasm && cd build.wasm
../scripts/configure -y wasm     # -y MUST be first — without a tty the profile
                                 # otherwise loops forever on a "Run Cmake?" prompt
ninja swipl-web

# 3. drop the three artifacts next to these files
cp src/swipl-web.js src/swipl-web.wasm src/swipl-web.data <this-dir>/
```

> **cmake note:** cmake 4.2.0–4.3.3 sits in emscripten's unsupported-shared-lib
> window and prints `does not support emscripten shared libraries` on *every*
> compile probe (seen loudly with cmake 4.2.3). It's **harmless** — the wasm
> kernel links statically and `swipl-web` builds fine — just don't mistake the
> warnings for errors. Configure runs hundreds of cross-compile probes (each
> spawns emcc+node), so it takes a couple of minutes; the compile is ~600 steps.

## Building the app `.qlf`

Fold the app into one Quick-Load File served alongside the artifacts. **QLF is
word-size-specific**, so the file the browser loads must be produced by the
*wasm* build's own swipl (run under node), **not** native `swipl` — a native
x86-64 qlf encodes 64-bit instruction variants/offsets and would load-and-
misbehave under wasm32 (the load check is version+VM-signature only, both
word-size-independent). Build it from the wasm tree:

```bash
# run the wasm build's node binary; NODERAWFS lets it read the app sources by
# absolute path. Point it at THIS repo's solve_browser.pl.
cd ~/src/swipl-devel/build.wasm
node src/swipl.js -g "qcompile('$OLDPWD/spikes/wasm-browser/solve_browser.pl', [include(user)]), halt" -t 'halt(1)'
cp src/solve_browser.qlf "$OLDPWD/spikes/wasm-browser/crosswordsmith.qlf"
```

> ✅ **Validated natively (64→64 only):** `include(user)` *does* fold all seven
> modules — a 46 KB `solve_browser.qlf` loaded on its own into a fresh native
> swipl runs `browser_selftest` cleanly (both paths). That proves the *packaging*
> works but **not** the word-size crossing — the residual gate is to produce the
> qlf with the wasm/node swipl (above) and confirm `Prolog.consult('./crosswordsmith.qlf')`
> loads it under WASM. Fallback if not: serve the source tree and `consult('./load.pl')`.

## Running

Serve this directory over HTTP (the worker + `.wasm`/`.data` need real URLs,
and `.wasm` should be served as `application/wasm`):

```bash
python3 -m http.server 8080   # then open http://localhost:8080/
```

Single-threaded build ⇒ **no COOP/COEP headers required.** Caveat: Python's
`http.server` may not send `Content-Type: application/wasm` for `.wasm` on
older stdlib versions — SWI's loader uses streaming instantiate, so if the
browser complains about the MIME type, serve with a dev server that sets it
(e.g. `npx serve`) or add the mapping.

## Known spike shortcuts (not production)

- Reaches module internals (`Module:Pred`); production adds a proper exported
  `solve_browser/2` in a small `browser.pl`.
- Grid size/mode read straight from the input dict; the CLI's `--max-size`
  crop / fragment / seed machinery is out of scope.
- `heartbeat: 50000` in `worker.js` is a guess — tune against a real search
  (the auto-yield-starvation risk in the plan).
