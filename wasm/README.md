# WASM browser client

Runs crosswordsmith's `arrange` solver in the browser under SWI-Prolog's WASM
runtime — no server, the solve happens client-side in a Web Worker. See
[`docs/plans/wasm-browser-deployment.md`](../docs/plans/wasm-browser-deployment.md)
for the full strategy and [`docs/research/swi-vm-wasm-performance.md`](../docs/research/swi-vm-wasm-performance.md)
for the VM cost model.

> ✅ **Validated end-to-end in headless Chrome (2026-07-05):** load the page,
> click Solve, a solved 5×5 grid renders — the solver runs off the main thread
> and calls the **same `arrange_solve/4`** the CLI's `arrange` verb calls (only
> the output sink differs). Promoted from a `spikes/` spike to this maintained
> layout the same day; the client code is still spike-grade pending the
> productionization items in plan §9 (an exported `browser.pl`, CLI resolver
> reuse, heartbeat tuning).

## Layout

| Path | What |
|---|---|
| `client/` | the browser app — served as static files |
| `build/build-wasm.sh` | reproducible build: emsdk + deps → `swipl-web` → app `.qlf` |
| `test/` | headless-Chrome regression tests (Playwright + system Chrome) |

### `client/`
| File | Role |
|---|---|
| `solve_browser.pl` | file-free entry points: `solve_browser_str/2` (**JSON string → JSON string**, what the worker calls), `solve_browser_json/2` (dict → JSON atom), `solve_browser/2` (dict→dict, illustrative). Loads the project via `../../load.pl`. |
| `worker.js` | the SWI-Prolog WASM instance, off the main thread. **JSON in, JSON out.** |
| `main.js` | page controller; "Stop" = `worker.terminate()` + warm-spare respawn. |
| `index.html` | minimal UI (words in, grid out). |
| `probe.html` | diagnostic page (library resolution / json autoload); use if the browser path breaks. |
| `swipl-web.{js,wasm,data}`, `crosswordsmith.qlf` | **build outputs — gitignored**, produced by `build/build-wasm.sh`. |

## Build

One script does the whole thing (emsdk activate → stage zlib/pcre2 → build
`swipl-web` from our pinned commit → copy artifacts → qcompile the app `.qlf`):

```bash
wasm/build/build-wasm.sh          # verifies the swipl-devel pin; won't move a shared HEAD
```

It's idempotent-ish (skips already-staged deps) and drops all outputs into
`client/` (gitignored). Overridable via env: `WASM_HOME`, `SWIPL_SRC`,
`SWIPL_ALLOW_CHECKOUT=1`. The prose rationale (emsdk 6.0.1 pin, no native friend,
qlf word-size, cmake-4.2.x noise) lives in plan §2 — the script is the executable
record of it.

> **Why the app `.qlf` must be wasm-produced:** QLF is word-size-specific. The
> file the browser loads must come from the *wasm* build's own swipl (run under
> node), **not** native `swipl` — a native x86-64 qlf encodes 64-bit VM
> variants/offsets and would load-and-misbehave under wasm32 (the load check is
> version+VM-signature only, both word-size-independent). `build-wasm.sh` step 4
> does this correctly; `include(user)` folds all seven project modules into one
> self-contained file (verified: it loads with no source tree reachable).

## Run

Serve `client/` over HTTP (the worker + `.wasm`/`.data` need real URLs):

```bash
( cd wasm/client && python3 -m http.server 8080 )   # open http://localhost:8080/
```

Single-threaded build ⇒ **no COOP/COEP headers required.** Caveat: some stdlib
`http.server` versions don't send `Content-Type: application/wasm`; SWI's loader
uses streaming instantiate, so if the browser complains, serve with something
that sets it (e.g. `npx serve`).

## Test

The Prolog side is testable under native swipl (no build needed):

```bash
swipl -q -g browser_selftest -t halt wasm/client/solve_browser.pl
# dict path OK: gridLength=5, 5 placed words
# json path OK: 1706 chars
# str  path OK: 1706 chars
```

The browser path is tested headlessly (Playwright drives **system Chrome** — no
browser download). Needs the server above running:

```bash
( cd wasm/test && PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 npm install )
node wasm/test/headless.mjs     # exits 0 on "placed …"; expect the RAT/ARC/CAT/TAR/CAR grid
node wasm/test/probe.mjs        # low-level load/json diagnostic (run first if things break)
```

## Browser gotchas (fixed in `client/worker.js` + `client/solve_browser.pl`)

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

## Known shortcuts (client is still spike-grade — plan §9 productionizes)

- Reaches module internals (`Module:Pred`); production adds a proper exported
  `solve_browser/2` in a small `browser.pl`.
- Grid size/mode read straight from the input dict; the CLI's `--max-size`
  crop / fragment / seed machinery is not wired in.
- `heartbeat: 50000` in `worker.js` is a guess — tune against a real (large)
  search (the auto-yield-starvation risk in plan §7).
