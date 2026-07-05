# Research: SWI-Prolog VM cost model + WASM deployment notes

Distilled from a research sweep (2026-07-04) for the arrange performance
campaign (see `docs/experiments.md`, "Crosswordsmith arrange campaign").
Grounded in: installed swipl-devel 10.1.10 library sources, live profiles of
`arrange_best_layout` on the ladder fixtures, and the SWI docs/discourse
sources cited inline. Findings that became experiments are cross-referenced.

## VM cost model — what the profiles actually mean

- **`has_type/2` attribution correction.** Profiled `error:has_type/2` /
  `must_be/2` cost in the arrange search is NOT `aggregate_all`'s option
  checking — `aggregate.pl` and `solution_sequences.pl` contain zero
  `must_be` calls; `aggregate_all(count, G, N)` compiles to a lean
  `nb_setarg` fail-driven loop. The real source is `list_to_set/2`, whose
  first line is `must_be(list, List)`, called per candidate pair inside
  `find_intersecting_word/6`. Inclusive-time profiling attributed it to the
  enclosing `aggregate_all` frame. Fix: an order-preserving private dedup
  (or eliminate the intersection machinery entirely via precomputed
  crossing triples — E-H3). Do NOT substitute `sort/2`: it reorders
  candidates and breaks golden byte-identity.
- **Counter replacement is a small lever.** Hand-rolled `nb_setarg` counter
  vs `aggregate_all(count, capped(2, ...))`: ~8-9% faster on the counting
  machinery in a 200k-rep micro-bench (SWI 10.1.10). `findall`+limit ~17%
  slower; `call_nth/2` ~28% slower. (Project history F007 already found
  aggregate_all beats findall+length in situ.)
- **assoc vs compound term.** `get_assoc/3` = O(log N) with `compare/3` per
  AVL node; `put_assoc/4` rebuilds the root-to-leaf path (allocation → GC
  pressure). Live profile: `assoc:adjust/5` + `$btree_find_node/5` ≈ 15% of
  self-time on a mid rung. `arg/3` is O(1), allocation-free. A grid-sized
  compound (441 args at 21x21; 2,500-arity verified fine) with either
  trailed `setarg/3` writes or unbound-var cells bound by unification is
  the right representation for a backtracking search (never `nb_setarg` —
  not undone on backtrack). → E-H2.
  - Term representation: https://eu.swi-prolog.org/download/publications/swi7.pdf
  - setarg semantics: https://www.swi-prolog.org/pldoc/man?section=setarg
- **Dicts cost ~1.5x compound-term access** (community-corrected figure; the
  circulated 5x number came from flawed methodology). SWI guidance: compound
  terms for low fixed-arity performance-sensitive records. The placed-word
  dicts (~94k `get_dict` calls, ~3.5% self-time on a mid rung) are a
  moderate, second-order target after the grid.
  https://swi-prolog.discourse.group/t/performance-cost-of-dictionaries/7779
- **`atom_chars/2` lists cost 2-3 cells/char** (SWI docs); `entry_letters/2`
  re-derives per word per recount node. Precompute at load. → E-H3.
- **Checked and NOT hot spots — don't spend effort there:**
  `call_with_inference_limit/3` (happy path = one catch frame + counter
  push/pop, called once per corner), `keysort`/`library(pairs)` (native
  stable sort + thin linear wrappers), `nth1/3` (C builtin, already adopted
  as F010), `intersection/3` at word-length scale, `statistics/2`, the one
  yall lambda (rescore path, once per corner candidate).
- **JIT clause indexing** is automatic (inspect via `jiti_list/0`), only
  relevant to clause-based predicates; irrelevant to an `arg/3` grid.
  https://www.swi-prolog.org/pldoc/man?section=jitindex

## WASM engine — deployment risks and guidance

- **Slowdown vs native:** ~4x (2022 self-report), regressed toward ~10x with
  newer Emscripten, fixed Sept-2025 by moving `setjmp()`/`longjmp()` off the
  hot VM path (swipl-devel `7b4b138ba`, 2025-09-29), back to ~2.7-3x.
  **CORRECTION (2026-07-05):** our native pin is NOT the bare 10.1.10 tag — it
  is `V10.1.10-17-gaa6289399` (commit 2026-07-01), and
  `git merge-base --is-ancestor 7b4b138ba aa6289399` confirms **the pin already
  contains the fix.** So building the WASM target *from our own tree* gives the
  fast path automatically; no separate post-fix checkout is needed. (Earlier
  text here assumed the bare tag and wrongly said the pin predates the fix.)
  Do NOT enable `WASM_EXCEPTIONS` chasing more speed: it is an *alternative* to
  this fix (`-fwasm-exceptions`), ~30% slower, and superseded — leave it OFF.
  https://swi-prolog.discourse.group/t/wasm-performance/9320 · see
  `docs/plans/wasm-browser-deployment.md`
- **Backtracking is disproportionately expensive under WASM**: the
  setjmp/longjmp emulation sits exactly on the backtracking/exception
  control path a DFS exercises. Tree-size reductions pay more than
  inference-count parity suggests.
- **Inference-count portability is NOT certified.** No documented guarantee
  the same program yields identical counts under WASM (GC and setjmp
  overhead are invisible to counts; indexing-driven choicepoint survival
  could in principle shift redo counts, though the indexing code has no
  `__EMSCRIPTEN__` conditionals). Cheap validation: run the existing ladder
  (`benchmarks/workloads.pl`) once under swipl-wasm via Node and diff
  against `benchmarks/baseline.json`. Do this before trusting the ratchet
  as a WASM proxy.
- **Startup:** `.qlf` is the documented fast path — `qcompile/2` with
  `include(user)` produces a single file loadable via `Prolog.consult()`
  from a URL. Author-measured: ~20x faster loading than source, ~50%
  smaller (s(CASP) 0.45s→0.17s; clpfd 0.11s→0.006s). Full `qsave_program`
  saved states: no evidence of WASM support (explicit negative finding).
  Ship qlf-only (`-DINSTALL_QLF=ON -DINSTALL_PROLOG_SRC=OFF`, ≥9.3.22);
  trade-off: no source-level tracing. No published numbers exist for WASM
  fetch/instantiate time.
  https://swi-prolog.discourse.group/t/support-precompiled-qlf-libraries/8887
  https://www.swi-prolog.org/pldoc/man?section=wasm-consult
- **Missing under WASM:** threads (`--disable-mt`; use multiple async
  engines on one core), ~21 libraries (socket, ssl, crypto, process,
  archive, bdb, janus, ...). GMP optional; LibBF fallback has considerably
  worse rational arithmetic (irrelevant to arrange's integer arithmetic; the
  WASM build sets `USE_GMP=OFF` as standard). **Tabling is available** (a core
  engine feature, not a package — CORRECTED 2026-07-05 from "unconfirmed"; the
  emscripten package set is `clpqr plunit chr clib http semweb pcre utf8proc`
  and `library(http/json)` resolves via a 2025 compat shim to the standalone
  `json` package). A smoke test is still worth it before depending on tabling
  (E-H3's variant), but it is not disabled.
  https://swi-prolog.discourse.group/t/notes-on-building-swi-prolog-for-webassembly-wasm/9488
- **Long searches block the browser event loop**: auto-yield happens only at
  the exit port ("recursive loops that call no other predicates only exit
  when all is done") — mitigate with `set_prolog_flag(heartbeat, N)`.
  Stack limits are set at runtime from JS
  (`Prolog.call("set_prolog_flag(stack_limit, N)")`), not Emscripten flags.
- **Build flags that matter:** `-s STACK_SIZE=1048576
  -s STACK_OVERFLOW_CHECK=1` (default 64k Emscripten stack has caused
  crashes), `ALLOW_MEMORY_GROWTH=1`, `WASM_BIGINT`, `LZ4`, Release build.
  **CORRECTED 2026-07-05:** these are already the in-tree defaults in
  `cmake/EmscriptenTargets.cmake` — no manual tuning. And **drop the old
  `-DVMI_FUNCTIONS=ON` (">10% faster")** advice: the tree now forces
  `DEFAULT_VMI_FUNCTIONS OFF` under `if(EMSCRIPTEN)`, so leave it at the
  default. https://swi-prolog.discourse.group/t/wasm-updates/6516
- **Word-size provenance (NEW 2026-07-05, adversarial review).** WASM is 32-bit;
  two build-time artifacts must be produced by a *same-pointer-size* Prolog and
  are easy to get wrong: (a) **no `-DSWIPL_NATIVE_FRIEND`** pointing at our x86-64
  build — a friend must match the *target* pointer size, so a 64-bit friend emits
  a mismatched `boot.prc`/`.qlf`; leave it unset and let the wasm binary self-
  bootstrap under node (`cmake/Ports.cmake`: "compatible ⇒ same pointer size").
  (b) The **app `.qlf` must be built with the wasm swipl** (`node
  build.wasm/src/swipl.js -g "qcompile(...)"`), not native `swipl`; qlf encodes
  word-sized instruction variants/offsets and the load check is
  version+VM-signature only (word-size-independent), so a native x86 qlf loads-
  and-misbehaves rather than failing cleanly. Native qlf validation proves 64→64
  packaging only, not the crossing.
  **BUILT 2026-07-05:** `swipl-web` was built from the pin with emsdk 6.0.1 (no
  friend); the wasm binary self-bootstrapped `boot.prc`+library qlf under node,
  confirming (a). `library(http/json)` + `json_write_dict` verified under `node
  src/swipl.js`, confirming `json` is pulled via the `http` dep (configure log:
  "Added packages json: Required by package http") — refuting a review claim that
  it was missing. Build gotchas: `scripts/configure` needs `-y` as its first arg
  (else it hangs on a tty-less confirm prompt); cmake 4.2.0–4.3.3 warns
  "does not support emscripten shared libraries" on every probe but the static
  wasm kernel builds fine. **RAN IN-BROWSER (headless Chrome, same day):** the
  spike qcompiles a self-contained app qlf (via node/wasm swipl), loads it in a
  Web Worker, and renders a solved grid. Three Worker-only fixes were required —
  `self.window=self` shim (SWI's URL helpers assume a browser `window`),
  absolute-URL `Prolog.consult`, and **JSON-in** (the JS→Prolog binding delivers
  JS strings as *atoms*, breaking `string`-typed checks — the input-side twin of
  the null/bool output quirk). In the web image the `library(http/json)` *alias*
  doesn't resolve (only `library(json)`), so explicit `use_module` warns, but json
  predicates autoload — so a genuine web/node divergence exists, just not a fatal
  one. See `docs/plans/wasm-browser-deployment.md` §5/§6.
- **Gotchas:** `statistics(cputime)` returns walltime under WASM (display
  only; the arrange budget is inference-based so cutoffs are unaffected);
  multiple SWIPL instances leak (no destroy(),
  https://github.com/SWI-Prolog/npm-swipl-wasm/issues/23); ~300MB memory
  ceiling observed on Android Chrome.

## Standing recommendations for the WASM milestone (outside this campaign)

The full strategy now lives in `docs/plans/wasm-browser-deployment.md` (with a
runnable spike under `spikes/wasm-browser/`). In brief:

1. Build `swipl-web` **from our own tree** (`~/src/swipl-devel`, exact-version
   parity — the pin already has the setjmp fix); qlf-only, `WASM_EXCEPTIONS` OFF.
2. Validate ladder inference-count parity under Node/swipl-wasm once.
3. Set heartbeat + stack_limit at init; single instance, reused; per-solve
   throwaway engine (`forEach({engine:true})`).
4. Input via query binding; output via JSON (`json_write_dict` + `JSON.parse`),
   not a bound dict — Prolog `null`/bool atoms come back as JS strings.
5. Run the solver in a Web Worker; cancel via `worker.terminate()`.
