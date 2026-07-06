# WASM browser client + SDK

Runs crosswordsmith's `arrange` solver in the browser under SWI-Prolog's WASM
runtime — no server, the solve happens client-side in a Web Worker, consumed
through a typed JS SDK facade. Strategy and wire contract:
[`docs/plans/wasm-sdk-strategy.md`](../docs/plans/wasm-sdk-strategy.md); the
build substrate: [`docs/plans/wasm-browser-deployment.md`](../docs/plans/wasm-browser-deployment.md);
VM cost model: [`docs/research/swi-vm-wasm-performance.md`](../docs/research/swi-vm-wasm-performance.md).

> ✅ **Validated end-to-end in headless Chrome (2026-07-05, re-validated on the
> SDK 2026-07-06):** the solver runs off the main thread and drives the **same
> arrange engine** the CLI's `arrange` verb calls. Gates cleared: **#2**
> inference parity (native↔wasm) and **#1** a real ~38M-inference search
> completing in-browser with the UI responsive and a prompt `terminate()`
> cancel. The spike's `Module:Pred` reaches are gone: the Prolog side is the
> exported `browser_dispatch/3` spine (`prolog/crosswordsmith/browser.pl`).

## Layout

| Path | What |
|---|---|
| `sdk/` | the JS/TS SDK: `crosswordsmith.mjs` (`createCrosswordsmith()` facade) + hand-written `crosswordsmith.d.ts` |
| `client/` | the engine assets — served as static files; `harness.html` is the minimal SDK consumer (smoke/diagnostic page) |
| `build/build-wasm.sh` | reproducible build: emsdk + deps → `swipl-web` → app `.qlf` |
| `test/` | native + headless-Chrome regression tests (Playwright + system Chrome); `run_all.sh` is the one-command battery (`make test-wasm`) |

### `sdk/`
| File | Role |
|---|---|
| `crosswordsmith.mjs` | `createCrosswordsmith()` — spawns + warms the engine Worker (plus a standby spare), id-correlated single-flight queue, `AbortSignal` cancel = `terminate()` + spare promote with generation fencing, engine-sourced `capabilities()`. `randomSeed()` draws the app-side regenerate seed in `[0, 1e9]`. |
| `crosswordsmith.d.ts` | the envelope + `ArrangeInput`/`Layout` types (DEC-5); locked against the goldens by `test/golden_type_check.mjs` |

### `client/`
| File | Role |
|---|---|
| `worker.js` | the SWI-Prolog WASM instance, off the main thread. Speaks the v1 envelope RPC: `{type:"request", request:{v,id,verb,params}}` → `{type:"response", id, envelope}`; **JSON in, JSON out**, verb/payload as **bound** forEach variables. |
| `solve_browser.pl` | the qcompile load root (thin adapter): loads `load.pl`, which exports `browser_dispatch/3` into `user` for the Worker's goal. |
| `harness.html` | minimal SDK consumer + test host (`?noauto` for the raw-worker probes) |
| `probe.html` | diagnostic page (library resolution / json autoload); use if the browser path breaks |
| `swipl-web.{js,wasm,data}`, `crosswordsmith.qlf` | **build outputs — gitignored**, produced by `build/build-wasm.sh` |

## Wire contract (v1)

One generic Prolog entry (`browser_dispatch/3`) dispatches every verb; replies
are a discriminated envelope (strategy §4):

```jsonc
{ "v": 1, "id": "r-1", "verb": "arrange", "status": "success", "result": { /* layout */ } }
{ "v": 1, "id": "r-1", "verb": "arrange", "status": "failure",
  "detail": { "reason": "no_interlock" | "grid_too_small" | "unplaceable_words", "words": [] } }
{ "v": 1, "id": "r-1", "verb": "arrange", "status": "error",
  "error": { "type": "validation" | "budget_exceeded" | "resource_exhausted"
           | "unknown_verb" | "unsupported_version" | "internal", "message": "…" } }
```

Verbs on the spine: **arrange**, **lint**, **export**, **capabilities**
(`fill` is deliberately not browserified yet — strategy §6 phase 4).

- **arrange** params (DEC-6): `size` (default 15), `mode` (`"fixed"` | `"max"`
  tight crop), `bestEffort`, `seed` (integer ≥ 0; not combinable with
  `bestEffort`).
- **lint** params: `layout` (a canonical arrange result), `profile`
  (`toc` | `blocked-uk` | `american` | `barred-ximenean`, required),
  `allowAsymmetry` (default false). The report is the success result — a FAIL
  verdict is an answer, not an error (the CLI's verdict-as-exit-code is a
  shell convention).
- **export** params: `layout` + `to` (`"ipuz"` → the ipuz v2 JSON document
  itself; `"exolve"` → `{format:"text", body}`).

"Cancelled" is a JS-side `CancelledError` rejection, never an envelope — a
terminated worker posts nothing.

## Build

One script does the whole thing (emsdk activate → stage zlib/pcre2 → build
`swipl-web` from our pinned commit → copy artifacts → qcompile the app `.qlf`):

```bash
wasm/build/build-wasm.sh          # verifies the swipl-devel pin; won't move a shared HEAD
```

It's idempotent-ish (skips already-staged deps) and drops all outputs into
`client/` (gitignored). Overridable via env: `WASM_HOME`, `SWIPL_SRC`,
`SWIPL_ALLOW_CHECKOUT=1`. The prose rationale (emsdk 6.0.1 pin, no native friend,
qlf word-size, cmake-4.2.x noise) lives in the deployment plan §2 — the script
is the executable record of it.

> **Why the app `.qlf` must be wasm-produced:** QLF is word-size-specific. The
> file the browser loads must come from the *wasm* build's own swipl (run under
> node), **not** native `swipl` — a native x86-64 qlf encodes 64-bit VM
> variants/offsets and would load-and-misbehave under wasm32 (the load check is
> version+VM-signature only, both word-size-independent). `build-wasm.sh` step 4
> does this correctly; `include(user)` folds all project modules into one
> self-contained file (verified: it loads with no source tree reachable).

## Run

Serve the `wasm/` tree over HTTP (`sdk/` and `client/` are siblings, and the
worker + `.wasm`/`.data` need real URLs):

```bash
( cd wasm && python3 -m http.server 8080 )
# open http://localhost:8080/client/harness.html
```

Single-threaded build ⇒ **no COOP/COEP headers required.** Caveat: some stdlib
`http.server` versions don't send `Content-Type: application/wasm`; SWI's loader
uses streaming instantiate, so if the browser complains, serve with something
that sets it (e.g. `npx serve`).

## Test

**One command** — from a checkout with the swipl-devel wasm build tree present
(built per `wasm/build/build-wasm.sh`; override `WASM_BUILD`/`WASM_SWIPL`):

```bash
make test-wasm                            # stage artifacts + run the full wasm battery
make test-wasm WASM_TEST_ARGS=--yield     # + the ~60s gate #1 responsiveness/cancel probe
```

`wasm/test/run_all.sh` copies `swipl-web.{js,wasm,data}` into `client/` and
re-qcompiles `crosswordsmith.qlf` **with the WASM node swipl** (see "Why the
app `.qlf` must be wasm-produced" above) — both skipped when nothing changed —
then runs the value golden, the type lock, and the headless-Chrome probes
behind a self-managed static server on a free port, aggregating exit codes
(any failing step ⇒ non-zero). One-time prereq — the tests drive **system
Chrome** (Playwright `channel:'chrome'`), so skip the browser download:

```bash
( cd wasm/test && PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 npm install )
```

The same battery layer by layer, cheapest first — useful when bisecting a
failure:

**Native (no build, no browser)** — the dispatch spine itself is plunit-locked
in `tests/browser.plt` (envelope fuzz, determinism/reset, edge matrix, value
parity vs the committed CLI goldens, CLI resolver parity): part of `make test`.
Plus a quick self-test of the client load root:

```bash
swipl -g browser_selftest -t halt wasm/client/solve_browser.pl
# dispatch OK: status success, gridLength=5, 5 placed words
```

**WASM value golden + same-instance determinism (node, no browser)** — the
same requests through `browser_dispatch/3` under the REAL wasm VM must
deep-equal the CLI's output (DEC-8: value-equality after parse, not bytes;
covers seedless / `seed:42` / the `max` crop / strict / best-effort / lint /
both export formats, all through ONE reused instance so the per-request reset
is proven under wasm):

```bash
wasm/test/value_golden.sh          # WASM_SWIPL=… to point at the build tree
```

Finding (2026-07-06): **the VM RNG is not portable — so the engine owns its
PRNG.** `set_random(seed(N))` seeds GMP's RNG natively but SWI's builtin RNG
under the `USE_GMP=OFF` wasm build — the same seed drew a different sequence,
making seeded layouts diverge CLI vs browser. Fixed by replacing the seeded
path's RNG with a module-owned portable splitmix64 (`core.pl`; known-answer
locked in `tests/arrange.plt`): the same seed now reproduces the same layout
on every build, and the harness deep-equals wasm seed:42 against the CLI
alongside determinism, perturbation, and provenance checks.

**Type drift lock** — the committed goldens must structurally satisfy the
hand-written `Layout` type (DEC-5):

```bash
node wasm/test/golden_type_check.mjs
```

**Headless browser (SDK end-to-end)** — drives `client/harness.html` through
the real facade: smoke render, envelope echo + typed validation, seed
determinism, overlapping-call concurrency (id routing), queued + in-flight
cancel (CancelledError + warm-spare recovery), the arrange → lint → export
composition, engine-sourced capabilities. Needs the server above running:

```bash
node wasm/test/headless.mjs
node wasm/test/probe.mjs        # low-level load/json diagnostic (run first if things break)
```

**Worker-level unhappy paths** — envelope taxonomy under the real wasm build,
fresh worker per case + same-worker liveness after:

```bash
node wasm/test/error_probe.mjs
# PASS validation:bad-seed         — typed validation envelope; worker survives
# PASS failure:infeasible-strict   — {status:failure, reason:unplaceable_words}
# PASS heavy-default-completes     — ~38M search under the 256MB cap → 36 words
# PASS heavy-capped-recoverable    — 300KB cap → {error.type:resource_exhausted}
#      (the cap is derived from a measured ~750KB search peak — see the
#       HEAVY_STACK_CAP comment in error_probe.mjs)
```

**Inference-count parity (gate #2)** — proves the ratchet is a valid wasm proxy:
the arrange search yields the *same* inference count native and under wasm. Run
the ladder under both VMs and diff — certified byte-identical on all 12 rungs
(2026-07-05):

```bash
swipl -q                                wasm/test/inference_parity.pl -- --heavy > /tmp/nat.csv  2>/dev/null
node ~/src/swipl-devel/build.wasm/src/swipl.js -q wasm/test/inference_parity.pl -- --heavy > /tmp/wasm.csv 2>/dev/null
diff /tmp/nat.csv /tmp/wasm.csv && echo "parity OK"     # (drop --heavy for just the 5 fast core rungs)
```

**Large-search responsiveness + cancel (gate #1)** — a real ~38.3M-inference
search in headless Chrome: completes, page stays responsive, `terminate()`
cancels promptly. Needs the server running:

```bash
node wasm/test/yield_probe.mjs
```

Finding: the **Worker boundary** (not the heartbeat) is what keeps the UI
responsive, and `worker.terminate()` is the only working cancel — so the
heartbeat is not load-bearing (kept at 50000, a few-% throughput tax, no UX
cost). Full write-up in the deployment plan §7 gate #1.

Finding (folded into `worker.js`): the stack cap **must be set inside the
`forEach` goal**. A `{engine:true}` search engine does *not* inherit the
parent's `stack_limit` (defaults to 1GB) — so an init-time cap never constrains
the search. Prepended to the query, an over-budget search throws a recoverable
`resource_error(stack)` → a `resource_exhausted` envelope (worker keeps
working) instead of a heap-growth tab abort.

## Browser gotchas (fixed in `client/worker.js` + the Prolog spine)

None of these show up under `node` — they are specific to the browser **Worker**:

1. **`self.window = self`** before `importScripts("./swipl-web.js")`. SWI's wasm
   URL helpers reach for `window` (absent in a Worker) → `window is not defined`
   and `consult` fails.
2. **Absolute URL for `Prolog.consult`** (`new URL("./crosswordsmith.qlf",
   self.location.href).href`) — no `window` ⇒ no base URL ⇒ relative path 404s.
3. **JSON in, not a query binding.** The JS→Prolog binding turns JS strings into
   Prolog **atoms**, so a clue's `answer` fails core's `string` check; and Prolog
   `null`/`true`/`false` return as the JS *strings* `"null"`/`"true"`/`"false"`.
   The worker sends `JSON.stringify(request)`; `browser_dispatch/3` parses it
   with `json_read_dict` and returns the envelope as one JSON **atom**.

Expected non-fatal noise: `source_sink library(http/json) does not exist` warnings
during qlf load — in the web image the `http/json` *alias* doesn't resolve, but
`json_write_dict`/`json_read_dict` autoload from `library(json)`, so it all works.
The warnings cite source locations *baked into the qlf* (provenance for error
messages) — not live fetches; the qlf is self-contained (`include(user)`). You may
also see a `crosswordsmith.qlf … net::ERR_ABORTED` line: a cancelled duplicate
request from SWI's URL consult — the qlf still loads, and results are correct.

## Scope / open items (strategy §8)

- `fill` is the one unbrowserified verb (strategy §6 phase 4 — hard: dict
  delivery, memory ceiling, fastrw-under-wasm). `lint` and `export` landed
  2026-07-06 as additive dispatch clauses, proving the spine generalises.
- npm packaging (asset-copy step, dual exports map, licensing manifest) is
  deferred: the SDK is consumed in-repo (OQ-1 decision 2026-07-06; OQ-5/OQ-6
  gate the publish).
- The JS↔qlf version assert / honoured-params echo (OQ-7) lands with
  version-stamped asset paths at packaging time.
