# Plan: crosswordsmith client-side SDK (Prolog/WASM → JS/TS)

Status: **strategy / pre-implementation** · Drafted 2026-07-06 · **Red-teamed
2026-07-06** (three adversarial agents — source-grounding, architecture,
completeness; findings folded in below and ledgered in §10, each marked [RT]).

Goal: turn the working WASM browser spike (arrange-only, runs in a Web Worker —
see [`wasm-browser-deployment.md`](./wasm-browser-deployment.md)) into a
forward-looking, importable **JavaScript/TypeScript SDK** whose first shipped
capability is `arrange` (list of words → aesthetically-pleasing grid, optionally
randomised, JSON out for an HTML/CSS renderer), but whose architecture lets
`lint` / `export` / `fill` land later as **purely additive** methods — never a
re-plumb.

Grounded in: the existing spike (`wasm/client/{worker.js,main.js,solve_browser.pl}`),
the version-matched WASM manual (`docs/reference/swi-manual/wasm-{loading,calling,js-call}.md`,
SWI 10.1.10), the module source under `prolog/crosswordsmith/`, and a **four-agent
research sweep + a three-agent red-team (both 2026-07-06)**. External sources at
the foot.

---

## 1. Bottom line

The research converged on a **four-layer spine**; the red-team then corrected two
load-bearing implementation claims (the entry seam and byte-identity — §10). It is
not exotic — it is the shape duckdb-wasm and trealla-js already ship, adapted to
our **single-threaded** SWI reality (no threads, no `SharedArrayBuffer`, no
COOP/COEP).

```
┌─ Consumer app (yours) ─────────────────────────────────────────────┐
│   const cw = await createCrosswordsmith()                           │
│   const res = await cw.arrange(input, { signal })                   │
└─────────────────────────────────────────────────────────────────────┘
      │  typed method-per-verb:  cw.arrange / .lint / .export / .fill
┌─ SDK facade (new) ─ Worker baked in · single-flight queue · abort→terminate ─┐
│   private RPC:  postMessage({ v, id, verb, params, meta })          │
└─────────────────────────────────────────────────────────────────────┘
      │  ONE generic dispatch · JSON string both ways · bound forEach vars
┌─ Prolog spine (new: prolog/crosswordsmith/browser.pl) ─────────────┐
│   browser_dispatch(Verb, PayloadAtom, JsonEnvelope)                 │
│     ├ reset 3 globals · json_read_dict · dispatch · classify        │
│     └ dispatch(arrange, …) → (Outcome, LayoutDict)  [not *_solve/N] │
└─────────────────────────────────────────────────────────────────────┘
      │  (Outcome atom + Dict) OR thrown term → discriminated envelope
   { v, id, verb, status: "success" | "failure" | "error", result? | detail? | error? }
```

**Decision summary.** Build `arrange` on this generic spine + discriminated
envelope *from day one* — it makes every later verb additive. Public API is
**typed method-per-verb over a private generic RPC** (never a raw `query(goal)` —
it leaks Prolog and can't be typed). Cancellation stays `worker.terminate()` +
warm-spare, dressed as an `AbortSignal`. Ship our own hand-written `.d.ts`. Bake
the Worker *into* the SDK — our one real differentiator (swipl-wasm and trealla
both punt it to the consumer; the spike already solved it, 0 ms main-thread drift).

**The single highest-value change:** a **discriminated result envelope**. Today
the worker collapses three distinct outcomes into one generic `"no layout (search
failed)"`. They split into `success` / `failure` (legitimately no layout) /
`error` (a typed exception) — mapping onto arrange's existing outcome atoms.

**[RT] Two corrections the red-team made load-bearing (details §10):**
1. **The envelope seam is NOT the exported `arrange_solve/4`.** That predicate
   emits JSON to stdout and *fails* on any non-`placed` outcome — it hands a
   caller nothing to wrap. The envelope must drive the **lower, currently
   unexported** layer: `arrange_best_layout/5` (returns `placed|not_proven|
   infeasible`) + `arrange_diag_layout_dict/5` (builds the dict *without*
   emitting). The slice must export these (removing the spike's `Module:Pred`
   leak), not just `doc_to_words/2`.
2. **`result` is NOT byte-identical to the golden CLI output when nested.**
   `json_write_dict` pretty-prints (`width(72), step(2)`); nesting under `result`
   shifts indentation and flips spaces→tabs past a depth threshold, and drops the
   trailing `\n`. So `result` is a live typed **dict** and correctness is tested
   by **value-equality after `JSON.parse`**, not a byte diff.

---

## 2. Goal & non-goals

**In scope (this plan governs):** the Prolog **dispatch spine** (`browser.pl`), the
**wire contract** (envelope + error taxonomy), the **JS/TS SDK facade** (packaging,
cancellation, typing), a **phased verb roadmap**, and the concrete
**arrange-first slice** (§7).

**Non-goals (deferred, tracked elsewhere):** real UI beyond a smoke-test consumer;
`fill` in the browser (deployment plan §11); build reproducibility / supply-chain
hardening (deployment §10.3 — but its asset-provenance item is the hook in §5.3
here); the arrange `fragment` / `candidates` / `enumerate` sub-modes (out of the
slice by decision; the spine takes them as additive params later).

---

## 3. The Prolog spine — `browser_dispatch/3`

One generic entry point, **not** per-predicate glue. Adding a verb = one dispatch
clause + one classifier row + one TS type; never new message wiring.

```prolog
:- module(crosswordsmith_browser, [ browser_dispatch/3 ]).

% browser_dispatch(+Verb:atom, +PayloadAtom, -JsonEnvelope:atom)
browser_dispatch(Verb, PayloadAtom, JsonEnvelope) :-
    reset_request_state,                                  % §3.2 — unconditional
    catch( ( parse_request(PayloadAtom, Req),
             dispatch(Verb, Req, Outcome) ),              % Outcome is a TAGGED term
           Err,
           Outcome = thrown(Err) ),                       % never let a throw escape
    outcome_envelope(Verb, Req, Outcome, Envelope),       % §4.2 — total classifier
    with_output_to(atom(JsonEnvelope),
                   json_write_dict(current_output, Envelope)).

% dispatch/3 MUST unify a tagged Outcome for EVERY input — never plain-fail.
% (The cross-cutting invariant §10 — a plain fail here becomes the generic
% "no layout" masquerade the envelope exists to kill.)
dispatch(arrange, Req, Outcome) :- arrange_browser(Req, Outcome).   % §7
dispatch(Verb, _, thrown(error(unknown_verb(Verb), _))).            % catch-all last
```

`arrange_browser/2` (in the slice) resolves options (validate-and-**throw** on bad
input), drives `arrange_best_layout/5` + `arrange_diag_layout_dict/5`, and unifies
one of `placed(Dict) | infeasible(Detail) | not_proven | best_effort(Dict)` — a
**tagged term, never a bare fail**. JSON-string in / JSON-string out is mandatory
(§4.3); verb + payload ride as **bound `forEach` variables**, never
string-concatenated. The client entry (`wasm/client/solve_browser.pl`) stays the
thin `user`-level adapter the Worker + `build-wasm.sh` name by string; it delegates
to `crosswordsmith_browser:browser_dispatch/3`. `load.pl` gains one
`:- use_module(crosswordsmith(browser)).` so the module folds into
`crosswordsmith.qlf` under `include(user)`.

### 3.1 Prerequisite A — the `cli_error`→throw refactor (partly pre-existing)

The CLI's **option/usage** resolvers report errors through `cli_error/1` →
`fail_with/2` = `format(user_error, …), fail` (`crosswordsmith:273–283`) — they
never throw, so a browser path that reuses them would swallow the reason as a bare
fail. **[RT]** But the **input-schema** layer already throws shaped `error/2`
terms the browser `catch` captures unchanged (`doc_to_words/2` throws
`json_no_clues_array` / `json_invalid_answer(Entry)` / `json_invalid_meta(_)`;
`check_unique_answers/1` throws `duplicate_answer(_)`; the fragment parser throws a
dozen shaped errors). So "validation → throw" is **partly free**; the real new work
is only the **size/seed/mode resolvers** (make them validate-and-throw, §4.2) and
reading the outcome atoms. Phase 0 for the whole SDK is refactoring the shared
option resolvers to throw with a CLI catch-and-`cli_error` adapter; the arrange
slice writes **local** throwing resolvers instead (§7, with a parity lock, §10).

### 3.2 Prerequisite B — per-request state reset (`reset_request_state/0`)

**[RT verified — this is the ONLY mutable global state.]** A repo-wide grep found
no `nb_setval`/`b_setval`/`set_flag`; exactly three `:- dynamic` module globals
survive a reused instance, and `{engine:true}` clears none:

| Dynamic (decl) | Setter | Note |
|---|---|---|
| `search_seed/1` (core.pl:479) | `set_search_seed/1` (core.pl:487) | also `set_random(seed(N))` — reseeds SWI's global RNG |
| `check_target_override/1` (arrange.pl:101) | `set_check_target/1` (arrange.pl:112) | the `--check-target` ceiling |
| `verbose_mode/0` (core.pl:142) | `set_verbose/1` (core.pl:144) | stderr chatter gate |

`reset_request_state :- set_search_seed(-1), set_check_target(-1),
set_verbose(false).`, run **unconditionally at the top of every dispatch** (all
three, defensively, even though the slice only exposes `seed`). **[RT] The reset is
provably sufficient:** the deterministic path is gated on the *presence of the
`search_seed/1` fact* (`seed_word_order/2`, `order_candidates/2`:
`( search_seed(_) -> shuffle ; identity )`), so a leftover RNG seed is never *read*
after the fact is retracted — no residue can leak.

### 3.3 Surface map — CLI stdout seam vs envelope seam

**[RT — recaptioned.]** Every verb's *exported* `*_solve/N` is a **stdout emitter
that `fail`s-to-stderr on the non-happy outcome** — unusable for the envelope,
which needs a nestable **dict** + a discriminated **outcome**. Those live in a
**lower, mostly-unexported** layer. The "envelope seam" column is the real target.

| Verb | CLI stdout seam (exported) | **Envelope seam** (mostly unexported) | Output | State reset | File-I/O to replace | Readiness |
|---|---|---|---|---|---|---|
| **arrange** | `arrange_solve/4` (emits+fails) | `arrange_best_layout/5` (Outcome) + `arrange_diag_layout_dict/5` (dict); best-effort: `arrange_best_effort/6` + `unplaceable_words/2` — **all unexported** | JSON dict | all three | `doc_to_words/2` (unexported); fragment `load_fragment/3`→pure `fragment_dict_words/3` | **spike reaches these as `Module:Pred`** — slice must export them |
| **lint** | `lint_solve/4` (emits) | `lint_dict_layout/3` (**unexported**) + `lint_run/5` (verdict in the report dict) | JSON | `verbose_mode` | the `open` in `lint_load/3` | **easy** — no engine/seed; export the dict-builder |
| **export** | `export_solve/2` (emits) | in-memory twin of `export_load/2` **incl. its shape validation** + `layout_to_ipuz/2` / `layout_to_exolve/2` | ipuz = JSON dict; **Exolve = plain text** | none | the `open` in `export_load/2` | **easy logic**, but load+validate+2 branches (not "one call") |
| **fill** | `fill_solve/4` (emits+fails) | fill lower layer (same emit/fail pattern) | JSON | `verbose_mode` | mask (`stockgrid_load/2` incl. **shape validation** `stockgrid_invalid`), dict (`load_dict/3`→`read_file_lines/2`), seeds, `fastrw` index, SHA-256 | **hard** — see §6 |

---

## 4. The wire contract — envelope + error taxonomy

### 4.1 Request / response envelope

Every verb routes through one request and one response shape, transported as a JSON
string across the Worker boundary.

```jsonc
// Request
{ "v": 1, "id": "r-42", "verb": "arrange",
  "params": { /* verb-specific — the CLI's per-verb JSON input */ },
  "meta":   { /* optional transport/engine concerns */ } }
```
```jsonc
// Response — discriminated on `status`, correlated by `id`
{ "v": 1, "id": "r-42", "verb": "arrange", "status": "success", "result": { /* live layout dict */ } }
{ "v": 1, "id": "r-42", "verb": "arrange", "status": "failure",
  "detail": { "reason": "no_interlock", "words": [] } }
{ "v": 1, "id": "r-42", "verb": "arrange", "status": "error",
  "error": { "type": "budget_exceeded", "message": "…", "detail": { } } }
```

**Field rationale.** `v` = envelope version, bumped **only** on a breaking
envelope-shape change. **[RT] `id`** is echoed on the response so the facade routes
each reply to its pending promise — mandatory for concurrent / rapid same-verb
calls, since `verb` alone can't disambiguate two in-flight `arrange`s. `verb` is
echoed for self-describing logs. `params` = the CLI's per-verb JSON input (reuses
`doc_to_words/2` etc.). `meta` carries transport/engine concerns (kept out of
`params` so the domain schema stays clean). `result` is a **live dict** (§1
correction); a text-producing verb (Exolve) uses `result: {"format":"text",
"body":"…"}` so the envelope is always valid JSON.

**Status ← arrange's outcomes** (strict is the slice's path; best-effort also
exposed):

| arrange outcome (atom / predicate) | envelope | notes |
|---|---|---|
| `placed(Dict)` — strict | `success` (`result` = Dict) | |
| `best_effort(Dict)` — best-effort with drops | `success` (`result` = Dict; drops in `result.diagnostics.arrange`) | **[RT]** drops object must be typed in the `.d.ts` |
| `infeasible(Detail)` — strict, no full interlock | `failure` (`detail.reason` ∈ `no_interlock`\|`grid_too_small`\|`unplaceable_words`; `detail.words` **only** when `unplaceable_words/2` is non-empty — [RT] the typical case has **no** word list) | not an error |
| `not_proven` — inference budget hit | `error` (`type: budget_exceeded`) | **[RT] synthesised from the outcome ATOM**, not a throw (see §4.2) |
| best-effort "nothing placeable" | `failure` (`reason: grid_too_small`) | **[RT]** was unmapped |

**Forward/backward compat + the stale-qlf trap [RT].** Additive-by-default (new
verb → new literal, old engine → `error/unknown_verb`; new result field → additive).
**But** the "old client + old *qlf* silently ignores a new optional param" case is a
*correctness* hazard: a new client sending `seed`/`checkTarget` to a long-cached old
qlf (deployment §5.3 warns qlfs get cached) gets a layout computed **without** it,
believing it applied. `capabilities()` returns verb names, not the per-verb param
surface, so it can't catch this. **Mitigation:** the engine **echoes the params it
actually honoured** (`result.meta.honoured` or similar) so a client can detect an
ignored key; and `createCrosswordsmith()` asserts JS-version == qlf-version at init
via the §5.3 version-stamped asset path. (OQ-7.)

### 4.2 Error taxonomy + the `outcome_envelope/3` classifier

**[RT — this section was materially wrong in the first draft and is rewritten.]**
Two facts corrected:

- **Budget is an ATOM, not a throw, and the engine already separates it.**
  `call_with_inference_limit/3` *binds-and-succeeds*; `construct_one/7` already
  branches `inference_limit_exceeded → budget`, and `arrange_best_layout/6` already
  maps that to the distinct atom `not_proven` vs `infeasible`. The collapse the
  envelope fixes happens **only** downstream in `arrange_solve/4` + the worker,
  which the spine bypasses. So there is **no "branch and throw in the search"** to
  do — the dispatcher reads the existing outcome atom. (The first draft's §4.2
  prescription would have patched the wrong layer.)
- **A bad seed FAILS, it does not throw.** `set_search_seed/1` is `integer(N),
  N>=0` — a negative/float/oversized seed makes it *fail*, and a failing resolver
  makes `dispatch` fail → the generic "no layout" masquerade returns. So the
  size/seed/mode resolvers **must validate-and-throw**, never rely on a builtin
  failing.

`outcome_envelope/3` is therefore a **total** function with two channels:

1. **Tagged outcome → status** (the atom map above).
2. **`thrown(Err)` → `error.type`** via an explicit, **total** table — unmatched
   ⇒ `internal` (never crash, never mislabel):

| Thrown term (real, existing) | `error.type` |
|---|---|
| `json_no_clues_array`, `json_invalid_answer(_)`, `json_invalid_meta(_)`, `duplicate_answer(_)`, `browser_missing_size`, local `bad_seed(_)`/`bad_size(_)`/`bad_mode(_)` | `validation` |
| `resource_error(stack)` (native, from the cap) | `resource_exhausted` |
| `unknown_verb(_)` | `unknown_verb` |
| `unsupported_version(_)` | `unsupported_version` |
| **anything else** | `internal` |

**`cancelled` is NOT a wire type. [RT]** After `worker.terminate()` the worker is
dead and posts nothing, so no `cancelled` envelope can ever arrive. It is a
**JS-side promise-rejection type only** (§5.2) — removed from the wire `CwError`
union.

**`resource_exhausted` is best-effort on mobile [RT].** The clean
`resource_error(stack)` only fires if the 256 MB cap trips *before* the device
memory ceiling (~300 MB); otherwise the tab may `abort()` uncatchably with no
envelope. Document it as "delivered when the cap trips first."

### 4.3 Why JSON-both-ways (the translation traps)

SWI's native binding is lossy both ways for our schemas (documented in
`solve_browser.pl:14–31`, re-confirmed against `wasm-calling.md`): a JS `String`
lands as a Prolog **atom** (breaking `doc_to_words/2`'s `string/1` guard), and any
Prolog atom returns as a JS `String` (so `null`/`true`/`false` cells corrupt to
`"null"`/`"true"`/`"false"`). `JSON.stringify` → `json_read_dict(…,
[default_tag(json)])` in, `json_write_dict` → `JSON.parse` out lands every value
with the right type. Return the envelope JSON as a Prolog **atom** (→ clean JS
`String`), not a Prolog string (→ a `Prolog.String` wrapper). Native binding is
reserved for a possible flat `version`/`capabilities` handshake only.

---

## 5. The JS/TS SDK facade

### 5.1 Shape — `createCrosswordsmith()`

```ts
const cw = await createCrosswordsmith({ assetBaseUrl? });   // spawns + warms the Worker
const res = await cw.arrange(input, { signal });            // typed envelope back
cw.version;            // { sdk, engine, qlf } — three distinct things (OQ-4)
await cw.capabilities();  // engine-sourced verb list, NOT a JS hardcode (OQ-4)
```

Each method is a thin wrapper over one private `postMessage({ v, id, verb, params,
meta })`. **[RT] The facade is self-serializing:** it keeps an `id → pending-promise`
map and a **single-flight FIFO queue** — overlapping `cw.arrange()` calls queue and
each reply routes by `id`; it does **not** "trust the page to serialize" or lean on
the empirical "a running solve ignores messages" behaviour. The warm-spare /
pre-warm logic moves *inside* the SDK. **Do not expose a generic `query(goal)`** —
it leaks Prolog and defeats typing.

### 5.2 Cancellation — `AbortSignal` over `terminate()`

`worker.terminate()` is the only reliable cancel for a yield-free CPU-bound DFS
(`Prolog.abort()` only unwinds a pending `await/2`; the heartbeat never drains the
Worker's message queue — measured). So `{ signal }` → `terminate()` + promote the
warm spare, and **reject every in-flight *and* queued promise on that instance** as
`cancelled` (terminate can't kill one request among several under the
disposable-instance model — state that semantics explicitly). **[RT] Backpressure &
fencing:** (a) a solve dispatched onto a not-yet-warm spare must *await* readiness,
not assume instant; (b) tag each instance with a generation id and **drop stray
`result`/`pong` messages from a terminated worker** (a late reply must not resolve a
`cancelled` promise); (c) coalesce cancel storms.

### 5.3 Packaging & assets

- **Lazy-fetch, don't inline.** Four sidecars (~3.9 MB) rule out base64 inlining.
  Resolve via `assetBaseUrl` defaulting to `new URL('./assets/', import.meta.url)`
  (Vite/webpack "just work"), CDN override, wired through SWI's `locateFile`.
- **[RT] The assets are built into a gitignored `build.wasm` tree** — packaging
  needs an explicit "copy built artifacts into the tarball" step and a ~4 MB
  package-size decision (OQ-5).
- **Cache-bust by version-stamped path** (`/assets/<pkgversion>/swipl-web.wasm`) —
  filenames fixed, path carries the version (Pyodide's `/vX.Y.Z/`). Also the fix
  for deployment §10.3 artifact-provenance, and the anchor for the §4.1
  JS↔qlf-version assert.
- **[RT] Redistribution licensing** (distinct from supply-chain hardening): the SWI
  runtime embedded in `swipl-web.{wasm,data}` needs BSD attribution; a future
  bundled dictionary for `fill` (UKACD18) ships a *verbatim* freeware licence
  (README.md). An npm-published SDK carries these obligations (OQ-6).
- **ESM-first + CJS fallback** via a dual `exports` map (tree-sitter's template).

### 5.4 TypeScript typing strategy

Discriminated union on `status`, `id`-correlated:

```ts
type Envelope<V extends string, R> =
  | { v: 1; id: string; verb: V; status: "success"; result: R }
  | { v: 1; id: string; verb: V; status: "failure"; detail?: FailureDetail }
  | { v: 1; id: string; verb: V; status: "error";   error: CwError };
type FailureDetail = { reason: "no_interlock" | "grid_too_small" | "unplaceable_words"; words?: string[] };
type CwError = { type: "validation" | "budget_exceeded" | "resource_exhausted"
               | "unknown_verb" | "unsupported_version" | "internal";   // NOT "cancelled"
               message: string; detail?: unknown };
// "cancelled" is a JS-side rejection type, never an envelope (§4.2/§5.2).
```

**Decision: hand-write the `.d.ts` now** (envelope + `ArrangeInput`/`Layout` incl.
the best-effort drops object), **adopt JSON-Schema-as-source-of-truth at
lint/export**. **[RT] But add a drift lock from day one:** a `tsd`/`expectType`
test asserting the committed golden fixtures (`tests/golden/arrange_*.json`)
structurally satisfy the hand-written `Layout` type (and/or seed the initial type
from a golden via quicktype). A few lines — it catches emitter↔type drift
immediately instead of letting the hand-written types free-run for arrange's whole
lifetime.

---

## 6. Phased roadmap to the full surface

| Phase | Deliverable | Difficulty | Notes |
|---|---|---|---|
| **0** | shared `cli_error`→throw refactor + `reset_request_state/0` | prerequisite | slice does a **local** throwing-resolver variant instead (§7) |
| **1** | **arrange** on the spine — §7 | low–med | envelope-seam exports + classifier; not the trivially-"reuse `arrange_solve/4`" the first draft implied |
| **2** | **lint** | **easy** | no engine/seed; export `lint_dict_layout/3`, dict-in adapter over `lint_run/5`, verdict already in the report dict |
| **3** | **export** | easy logic | in-memory twin of `export_load/2` **+ its shape validation**; ipuz → dict, Exolve → `{format:"text",body}` |
| **4** | **fill** | **hard** | **[RT] beyond the three named risks** (mask/dict/seeds file seams, `fastrw`-under-wasm → maybe `.qlf`, precomputed bignum masks): also (a) the *same* stdout-emit/fail outcome-channel lift as arrange (§4.2), (b) a real UKACD18 dictionary's memory footprint under the ~300 MB wasm ceiling **and** how the multi-MB dict is *delivered* (a second large fetched asset), (c) re-implementing `stockgrid_load/2`'s mask **shape validation** in-memory, not just `mask_white_cells/3`. |

After arrange, **lint is the cheapest second verb** — a good proof the spine
generalises before fill.

*(Phases 2–3 landed 2026-07-06: `lint` and `export` are dispatch clauses on the
spine — thin resolvers over `lint_dict_layout/3`+`lint_run/5` and
`export_layout_dict/2`+`layout_to_*/2` — golden-locked native and wasm against
the committed CLI goldens. `fill` is the only unbrowserified verb.)*

---

## 7. The arrange-first slice (landed 2026-07-06)

Product shape (decided with the user): **strict path is the randomisable one**
(`seed` = reproducible; the app draws a random seed in JS for "regenerate"), **both
`fixed` and `max` framing exposed** (`max` = tight square crop to side max(H,W) —
usually nicer for an HTML/CSS renderer than a fixed N×N with dead border rows;
default `fixed`), **skip** fragment/candidates/enumerate. `checkTarget` is
**deferred** (a power-user knob the words→grid→regenerate UX doesn't need; the reset
still covers `check_target_override` defensively). Exposed params: **`size`, `mode`
(`fixed`|`max`), `bestEffort`, `seed`**.

**Prolog**
- New `prolog/crosswordsmith/browser.pl` with `browser_dispatch/3` + the `arrange`
  clause; `reset_request_state/0`; **local size/mode/seed resolvers that
  validate-and-throw** (`bad_size`/`bad_mode`/`bad_seed`; `mode` ∈ {`fixed`,`max`} —
  un-block the spike's `browser_max_mode_unsupported`, routing `max` through
  `arrange_diag_layout_dict/5`'s `SizeMode` arg), guarding `seed`+`bestEffort` →
  throw (randomisation perturbs the strict search only — matches `guard_seed_combos`).
- **[RT] Export the envelope seam** (not just `doc_to_words/2`): `doc_to_words/2`
  (core), and from arrange `arrange_best_layout/5`, `arrange_diag_layout_dict/5`,
  `arrange_best_effort/6`, `unplaceable_words/2` — or an explicit new exported
  `arrange_outcome(+Words,+Size,+Mode,-Outcome,-Dict)` wrapper that returns the
  tagged term + dict (preferred: keeps the public surface one predicate and removes
  the spike's `Module:Pred` leak).
- `arrange_browser/2` maps `placed→success` / `infeasible→failure(reason,words?)` /
  `not_proven→error(budget_exceeded)` / best-effort success+drops (§4.2).
- `load.pl` +1 line; `solve_browser.pl` → thin `user`-level adapter (Worker + build
  untouched by name).

**JS**
- `worker.js`: generalise the `forEach` goal to `browser_dispatch(Verb, Payload, Json)`
  with `Verb`/`Payload` bound; keep the stack cap, warm, terminate.
- New `createCrosswordsmith()` facade (`cw.arrange`, `version`, `capabilities`) with
  the Worker + warm-spare + **`id`-correlated single-flight queue** (§5.1) inside,
  and `{signal}` → terminate cancel with generation-fencing (§5.2).
- **[RT] Draw the JS seed in `[0, 1e9]` integer** (matching the CLI's
  `random_between(0, 1_000_000_000, N)`); a float / `>2^53` value risks the
  `integer/1` guard and BigInt-translation ambiguity.
- Hand-written `crosswordsmith.d.ts` (envelope + `ArrangeInput`/`Layout` + drops).
- **[RT] Decide the demo-client fate:** either rewrite `main.js`/`index.html` as the
  first `createCrosswordsmith()` consumer (doubles as the browser smoke test — the
  spike currently speaks the *old* `{type:"solve"}` protocol and would otherwise
  rot), or retire it. Pick one in the PR.

**Tests** (the red-team's headline: the envelope is only as good as "every resolver
either tags an Outcome or throws a *mapped* term — never plain-fails"):
- **Malformed-params invariant fuzz [RT]** — non-int/negative/huge `seed`, `size`
  0/neg/non-int, missing `clues`, non-string `answer`, duplicate answers,
  `bestEffort`+`seed` — each must yield a **typed `error`/`failure` envelope, never
  `"no layout (search failed)"`. (Include an invalid `mode` in the fuzz set.) This
  is the regression lock for the whole design.
- **Same-worker determinism** — `seed:42` then seedless; seedless == seedless-in-
  isolation; `seed:42`×2 identical (the reset lock).
- **wasm-side value-correctness golden [RT]** — run the entry under `node
  build.wasm/src/swipl.js`; **`JSON.parse` both sides and deep-equal** `result`
  against the native golden (seedless + seed + **`max` crop**) — NOT a byte diff
  (§1/§10). Still catches the `null`→`"null"` corruption (a value difference) that
  motivated it. Exposing `max` adds a crop-path golden fixture.
- **Edge-case matrix [RT]** — empty word list, single word, all-unplaceable,
  unicode/lowercase answers — asserting the *envelope*, not "contains placed".
- **Resolver-parity matrix [RT]** — invalid `seed`/`size`/combos reject *identically*
  on the CLI and browser paths (locks the accepted-for-now resolver duplication, §10).
- **Concurrency + cancel [RT]** — two overlapping `arrange` calls each resolve to
  their own result (by `id`); a cancel with a queued sibling behaves per §5.2.
- **Golden-satisfies-type** — §5.4 drift lock.

Explicitly, this PR closes deployment §10.4's *wasm output-correctness* and
*same-worker determinism* debts and **leaves open** concurrency-under-load and
cancel-mid-render hardening beyond the tests above.

**Verification** rebuilds `crosswordsmith.qlf` with the wasm/node binary and runs
the wasm golden + `headless.mjs`; depends on the local `~/src/swipl-devel/build.wasm`
being intact — confirm first.

---

## 8. Decisions & open questions

**Decided (this plan):**
- **DEC-1:** generic `browser_dispatch/3` spine.
- **DEC-2:** three-way discriminated envelope (`success`/`failure`/`error`), `id`-
  correlated.
- **DEC-3:** typed method-per-verb over a private RPC; no exposed `query(goal)`.
- **DEC-4:** Worker baked in; cancel = `AbortSignal` → `terminate()` + warm-spare,
  generation-fenced.
- **DEC-5:** hand-written `.d.ts` now + a golden-satisfies-type drift lock; JSON-
  Schema-SoT at lint/export.
- **DEC-6 [RT-corrected; `max` re-added by user]:** arrange slice = strict+`seed`
  randomisation, **`fixed`+`max` framing** (`max` = tight crop, default `fixed`),
  **no `checkTarget`**; params `size`/`mode`/`bestEffort`/`seed`.
- **DEC-7 [RT]:** the envelope seam is `arrange_best_layout/5` + `arrange_diag_
  layout_dict/5` (exported via a wrapper), **not** `arrange_solve/4`.
- **DEC-8 [RT]:** `result` is a live typed dict; correctness is **value-equality
  after `JSON.parse`**, not byte-identity.

**Open:**
- **OQ-1 — decided (2026-07-06, user):** in-repo module — `wasm/sdk/` with
  `crosswordsmith.mjs` + `crosswordsmith.d.ts`, imported directly (serve `wasm/` as
  the web root). npm packaging deferred; OQ-5/OQ-6 gate the eventual publish.
- **OQ-2 — decided (2026-07-06):** `seed` lives in `params` (landed with the
  arrange slice).
- **OQ-3:** fill's `fastrw`-index vs `.qlf` carrier under wasm (measurement, Phase 4).
- **OQ-4 [RT] — mostly landed (2026-07-06):** `capabilities()` is engine-sourced (a
  `capabilities` verb on the spine) and `cw.version = {sdk, engine}`; the remaining
  piece — a qlf build-hash in the capabilities reply — lands with OQ-7 at packaging.
- **OQ-5 [RT]:** CI qlf rebuild reproducibility (deployment §10.3 "not standalone-CI-
  runnable") + the asset-copy-into-tarball / package-size step.
- **OQ-6 [RT] — RESOLVED (2026-07-13):** `wasm/THIRD_PARTY_NOTICES.md` is now the
  self-contained redistribution notice: verbatim licence texts verified against the
  pinned sources, covering the engine, its vendored components (incl. LGPL
  library(isub) with a static-link compliance statement), zlib, PCRE2, the
  Emscripten runtime (musl, compiler-rt, dlmalloc, MiniLZ4), and crosswordsmith;
  `stamp-manifest.sh` copies it into `wasm/client/` and `build-manifest.json`'s
  `licenses` names the sibling file. The UKACD dict obligation stays recorded for
  the day `fill` browserifies.
- **OQ-7 [RT]:** stale-qlf param drift mitigation — engine echoes honoured params
  and/or JS↔qlf version assert at init (§4.1).
- **OQ-8 [RT] — RESOLVED (2026-07-06):** both halves closed. (a) Variety probed:
  50 consecutive seeds → **20 distinct** layouts (bundled-17, 6 words / 17×17),
  **50/50 distinct** (toc-demo, 25×25), 8 distinct (toy 5-word / 5×5 — a tiny
  solution space owns few layouts); no pathological collapse, "regenerate" is
  honest. (b) Reproducibility: the engine now **owns a portable PRNG** (splitmix64
  — see the §10 implementation finding), so a seed reproduces the same layout on
  every build/platform; the remaining scope is engine-*version* (heuristic changes
  re-map seeds — record `cw.version.engine` alongside a kept seed).

---

## 9. Relationship to the WASM deployment plan

This plan **supersedes and generalises** [`wasm-browser-deployment.md`](./wasm-browser-deployment.md)
§9.1: that arrange-only resolver-lift becomes **Phase 1 here**, reframed as the
first verb of a generic SDK spine. Two corrections carried over: §9.1's "parked
branch" footer is **stale** (the WASM work is on `main`, the fill campaign is
closed, the tree is clean); and its §10.2 tracked debt (per-request reset, error
channel, resolver duplication) is **absorbed** here as first-class prerequisites.
The deployment plan's §2–§8 (build recipe, dependency staging, gates #1/#2, the qlf
word-size provenance rule) stands unchanged as the substrate.

---

## 10. Red-team pass (2026-07-06)

Three read-only adversarial agents stress-tested the first draft on distinct lenses
— **source-grounding** (re-verify every codebase claim against source), **architecture**
(attack the design), **completeness** (find omissions/over-promises). All three
independently reproduced the two HIGH findings; the architecture agent empirically
serialized a nested dict to prove the byte-identity bust. Ledger (resolution folded
into the sections above, marked [RT]):

**High — folded (blocking, were load-bearing errors):**
- **The envelope seam.** The exported `*_solve/N` predicates emit to stdout and
  `fail` on the non-happy outcome; the envelope needs the unexported `arrange_best_
  layout/5` (outcome) + `arrange_diag_layout_dict/5` (dict). §3.3/§7/DEC-7 corrected;
  the slice now exports the real seam (a wrapper) rather than perpetuating the
  spike's `Module:Pred` leak.
- **Byte-identity is false.** `json_write_dict` pretty-prints; nesting shifts indent
  and flips spaces→tabs and drops the trailing `\n`. §4.3/§7/DEC-8: `result` is a
  live dict, tested by value-equality after `JSON.parse`.

**Medium — folded:**
- **Error taxonomy rewrite** — budget is an *atom* the engine already separates (not
  a throw to add); a bad seed *fails* (must validate-and-throw); `cancelled` is not a
  wire type; the union needs an `internal` fallback + a total throw→type map over the
  *real* existing terms. (§4.2.) The first draft's "branch-and-throw in the search"
  prescription pointed at the wrong layer.
- **No correlation id** → concurrent replies unrouteable; added `id` + a facade
  single-flight queue. (§4.1/§5.1.)
- **Stale-qlf param drift** — an old cached qlf silently ignores a new param. Added
  the honoured-params echo + version assert. (§4.1/OQ-7.)
- **`failure.detail` over-promise** — culprit words are usually absent; added a
  `reason` enum and mapped best-effort's "nothing placeable". (§4.1.)
- **Resolver duplication with no parity lock** — accepted-until-Phase-B + a parity
  test matrix. (§7.)
- **Edge-case tests** — the slice closed only 2 of 5 §10.4 debts; added the invariant
  fuzz + edge matrix + concurrency/cancel tests. (§7.)
- **OQ gaps** — version/capabilities provenance, CI qlf rebuild + asset packaging,
  redistribution licensing. (OQ-4/5/6.)

**Low — folded:** randomization variety-collapse + engine-build-scoped
reproducibility + seed hygiene `[0,1e9]` (OQ-8/§7); `resource_exhausted` best-effort
on mobile (§4.2); warm-spare backpressure + stale-message fencing (§5.2); the day-one
type-drift lock (§5.4); demo-client fate (§7); lint/export "easy" export-gap nuance
(§3.3/§6); §4.1 strict-vs-best-effort mechanism distinction (§4.1); §3.1 overstated
greenfield / §3.2 line-cite nits (corrected in place).

**Verified as holding (spend no worry):** the three-dynamics reset is complete (no
`nb_setval`/`b_setval`/`set_flag` anywhere) and sufficient (deterministic path is
fact-gated, RNG residue never read); the `cli_error → fail_with` channel never
throws; `doc_to_words/2` unexported, `mask_white_cells/3` exported, `lint_run/5`
pure, Exolve emits text, `fragment_dict_words/3` pure — all confirmed.

**Implementation finding (2026-07-06, arrange slice):** cross-VM *seeded* equality
is impossible through the VM RNG: `set_random(seed(N))` seeds GMP's RNG natively but
SWI's builtin RNG under the `USE_GMP=OFF` wasm build — the same seed draws a
different sequence (verified: first `random/1` after `seed(42)` is 0.7846… native
vs 0.3745… wasm). **Resolved by owning the algorithm** (same day, while the seed
feature was unreleased and the compat break free): the seeded path now draws from a
module-owned portable splitmix64 (`core.pl` — pure unbounded-integer arithmetic,
bit-identical under GMP and LibBF; known-answer locked against the published test
vectors in `tests/arrange.plt`), with `library(random)` retained only for
`--shuffle`'s entropy draw. OQ-8's reproducibility scope thus narrows from
engine-*build* to engine-*version* (search-heuristic changes still re-map seeds).
`wasm/test/value_golden.sh` deep-equals wasm `seed:42` against the CLI, plus
native seam parity, determinism, perturbation (seeded ≠ seedless), and provenance;
`crosswordsmith.d.ts` documents the scope on `seed`.

---

## Sources

Research provenance: a four-agent research sweep + a three-agent red-team, both
2026-07-06. Load-bearing external anchors:

- SWI WASM manual (version-matched, local): `docs/reference/swi-manual/wasm-{loading,calling,js-call}.md`;
  online https://www.swi-prolog.org/pldoc/man?section=wasm-calling
- `swipl-wasm` (substrate): https://github.com/SWI-Prolog/npm-swipl-wasm
- trealla-js (`status` envelope; generic `query<T>()`): https://github.com/guregu/trealla-js
- duckdb-wasm (Worker-behind-async-facade; `selectBundle`): https://github.com/duckdb/duckdb-wasm
- web-tree-sitter (typed OO API; dual `exports` map): https://github.com/tree-sitter/tree-sitter/tree/master/lib/binding_web
- Pyodide (`indexURL` versioned assets): https://pyodide.org/en/stable/usage/api/js-api.html
- `json-schema-to-typescript` (the schema-SoT route, §5.4): https://github.com/bcherny/json-schema-to-typescript
- In-repo grounding: `wasm/client/{worker.js,main.js,solve_browser.pl}`,
  `prolog/crosswordsmith/{core,arrange,lint,export,fill,stockgrid}.pl`, the
  `crosswordsmith` CLI script, `docs/json-{input,output}-spec.md`,
  `tests/golden/arrange_*.json`.
