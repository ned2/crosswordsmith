# Plan: smaller and faster browser payload

Status: **in progress** · Drafted 2026-07-14. Landed: Phase 0 (benchmark
harness + committed baseline, `make bench-wasm-payload*`), Phase 1
(browser-specific load root; qlf 71,580 → 51,068 bytes raw, fill/stockgrid/
sha/fastrw out of the closure and gated in `run_all.sh`, the historic
`http/json` source_sink load warnings gone), and Phase 2 (minimal preload
image: strace-derived 22-file keep-list in `wasm/build/preload-profile.txt`,
deterministic link-command replay in `build-wasm.sh` step 2.6 producing
`src/profile-crosswordsmith-web/`, provenance in the manifest's
`preloadProfile`; swipl-web.data 1,646,186 → 195,323 bytes raw. Trio brotli
1,680,777 → 827,286 (−50.8%) and sdk readiness median 118.9ms → 95.9ms (−19%)
— both headline targets already met at Phase 2) — all 2026-07-14.

Goal: make the crosswordsmith browser engine substantially cheaper to fetch,
compile, and initialise without changing solver semantics, native behaviour, or
the SDK's default cancellation guarantees. The immediate consumer is ned.sh;
its site-specific rollout and placeholder-removal gate live in
[`ned2/ned.sh: docs/09-crossword-load-performance.md`](https://github.com/ned2/ned.sh/blob/main/docs/09-crossword-load-performance.md).

This is the canonical technical plan. The existing
[`wasm-browser-deployment.md`](wasm-browser-deployment.md) remains the record of
how the browser port works, and [`wasm-sdk-strategy.md`](wasm-sdk-strategy.md)
remains the SDK contract. This plan owns only payload/startup optimisation.

---

## 1. Outcome and boundaries

The production browser build should become a crosswordsmith-specific SWI-Prolog
distribution rather than the current general-purpose SWI WASM image.

Target outcomes:

- at least **50% less Brotli transfer** for a correct cold load;
- at least **15% faster median cold engine readiness** on the reference desktop
  benchmark, with no regression on a representative lower-powered device;
- the same `arrange`, `lint`, `export`, and `capabilities` envelopes and values;
- the existing cancellation, recovery, provenance, and notice guarantees;
- a deterministic build whose manifest describes the reduced package surface.

Non-goals:

- changing the arrange algorithm or its value/performance characteristics;
- changing the native CLI package surface;
- removing failure fallback or `AbortSignal` support;
- weakening the generic SDK's default eager-spare cancellation behaviour;
- hiding latency with a more elaborate loading animation. Presentation belongs
  to each consumer.

## 2. Measured baseline

Measurements were taken on 2026-07-13 from the pinned production artifact set
currently vendored by ned.sh. Raw sizes are exact bytes. Gzip level 9 and
Brotli quality 11 are diagnostic comparisons, not a promise that a production
host will emit byte-identical streams.

| Artifact | Raw | Gzip | Brotli |
|---|---:|---:|---:|
| `worker.js` | 8,635 | 3,797 | 3,251 |
| `build-manifest.json` | 1,929 | 993 | 915 |
| `swipl-web.js` | 192,096 | 51,258 | 43,682 |
| `swipl-web.wasm` | 2,195,615 | 796,588 | 628,861 |
| `swipl-web.data` | 1,646,186 | 1,179,939 | 1,008,234 |
| `crosswordsmith.qlf` | 71,580 | 32,377 | 26,904 |

One copy of the six files is 4.12 MB raw or 1.71 MB Brotli. Browser request
behaviour complicates the real total:

- `createCrosswordsmith()` eagerly starts an active worker and a warmed spare;
- without explicit cacheability, a counting server observed both workers
  transferring every large engine artifact: about **3.48 MB Brotli**;
- with `public, max-age=31536000, immutable` on hashed artifacts, Chrome fetched
  the hashed JS/WASM/data once, but still requested `worker.js` and the manifest
  twice and produced four QLF requests (the known duplicate consult/probe per
  worker): about **1.80 MB Brotli**;
- a single worker with the current duplicate QLF behaviour is about **1.74 MB
  Brotli**.

On localhost, fresh Chrome processes measured the current single-worker
`capabilities` readiness at a 109 ms median. A first five-word arrange is only
about 7–11 ms after boot, so payload/runtime startup is the meaningful engine
cost for small browser consumers.

## 3. Feasibility prototype and what it proves

A disposable prototype was assembled under `/tmp`; neither repository nor the
pinned SWI checkout was modified. It used:

1. a 225 KB allowlisted preload tree instead of copying the entire SWI library;
2. only the foreign libraries exercised by the current browser QLF, with no-op
   registration stubs standing in for 15 unrelated statically linked plugins;
3. the existing app QLF and unchanged worker/SDK protocol.

The resulting engine trio was:

| Artifact | Current raw | Prototype raw | Current Brotli | Prototype Brotli |
|---|---:|---:|---:|---:|
| JS | 192,096 | 153,007 | 43,682 | 35,330 |
| WASM | 2,195,615 | 1,405,233 | 628,861 | 447,045 |
| data | 1,646,186 | 198,930 | 1,008,234 | 165,004 |
| **Total** | **4,033,897** | **1,757,170** | **1,680,777** | **647,379** |

That is a 61% reduction in the compressed engine trio. Under the current
two-worker request pattern with correct immutable caching, the estimated cold
transfer falls from 1.80 MB to **763 KB Brotli** (57% smaller).

Ten alternating fresh-Chrome samples measured median engine readiness at
109.1 ms current versus **88.5 ms** for the fully trimmed prototype (19%
faster). Localhost suppresses the network benefit, so real cold-network gains
should primarily come from the removed megabyte of transfer.

Most importantly, the prototype passed the full `wasm/test/headless.mjs`
browser suite:

- smoke solve and rendered grid;
- typed envelopes and validation errors;
- deterministic seeded results;
- overlapping request routing;
- queued and in-flight cancellation plus spare recovery;
- heavy solve completion;
- lint, IPUZ export, Exolve export, and capabilities;
- random-seed range hygiene.

This proves the reduced runtime surface is viable. It does **not** make the
temporary allowlist or no-op stubs an acceptable production build; the phases
below replace them with tracked build inputs and full provenance.

## 4. Phase 0 — make the benchmark reproducible

Before changing production artifacts, add a benchmark/report command under
`wasm/test/` that records:

- raw, gzip, and Brotli size of every shipped artifact;
- the aggregate one-worker and eager-two-worker transfer models;
- actual request counts through a local cache-aware counting server;
- fresh-browser time from worker construction to `capabilities`;
- time for the first small arrange after readiness;
- peak/settled process memory where Chrome exposes a stable measurement;
- toolchain, SWI commit, package list, and build ID.

Required modes:

- current SDK (`eager` spare);
- raw single worker;
- cold browser profile per sample, not repeated calls in a shared HTTP/WASM
  cache;
- cache headers absent and correctly configured;
- at least one representative lower-powered machine or device before release.

Commit the benchmark script and a machine-readable baseline, but do not turn a
noisy wall-clock threshold into CI until variance is characterised. CI should
immediately gate deterministic artifact sizes and regression-suite correctness.

## 5. Phase 1 — browser-specific Prolog load root

The current [`solve_browser.pl`](../../wasm/client/solve_browser.pl) reaches the
global `load.pl`, which deliberately imports the whole product, including
`stockgrid` and `fill`. The browser contract exposes only `arrange`, `lint`,
`export`, and `capabilities`.

Create a browser load root that imports exactly:

- `core`;
- `metrics`;
- `arrange`;
- `lint`;
- `export`;
- `browser`.

It must remain location-independent for QLF creation without pulling in
file/path helpers solely to locate the general loader. Prefer explicit
module-relative loads owned by the browser root. Keep `load.pl` unchanged for
the CLI, tests, benchmarks, `stockgrid`, and `fill`.

Acceptance:

- QLF creation still uses the wasm32 node build, never native SWI;
- all four browser verbs pass value goldens against the CLI;
- `make test` remains unchanged and green;
- the complete browser suite remains green;
- the QLF dependency inventory no longer includes fill-only `sha`, `fastrw`,
  dictionary, or stock-grid paths.

The QLF itself is only about 72 KB raw, so direct savings are modest. This
phase matters because it defines the dependency closure for the larger preload
and package reductions.

## 6. Phase 2 — deterministic minimal preload image

SWI's current `EmscriptenTargets.cmake` copies the complete
`SWIPL_BUILD_LIBRARY` into `wasm-preload`. Introduce a supported way for this
build to supply a reduced preload manifest/profile while preserving the current
full image as SWI's default.

Preferred design:

- a CMake cache input or tracked manifest enumerates the files/directories to
  preload;
- the default remains the full SWI library tree;
- crosswordsmith's build selects a named `crosswordsmith-web` profile;
- the profile is assembled into a clean staging directory and passed to
  `--preload-file` deterministically;
- missing entries fail at build/test time rather than falling back to source or
  the host filesystem;
- the build manifest records the preload profile and a digest of its file list.

The successful prototype needed `boot.prc`, `ABI`, the WASM runtime helpers,
the JSON library, and the core libraries reached by the current QLF. Treat that
list as discovery evidence, not the final manifest: regenerate it after Phase 1
and validate every autoload/dynamic path through tests.

Supply-chain constraint: do not silently dirty a shared pinned SWI checkout.
Prefer an upstreamable CMake option in the pinned source. If a local overlay is
unavoidable, keep it as a tracked patch, apply it only to a disposable checkout,
record its digest in `build-manifest.json`, and teach the pin guards about that
explicit state.

Acceptance:

- `.data` is below **250 KB raw** on the current pin;
- no `source_sink`/autoload warning occurs in the full browser suite;
- error, cancellation, and recovery probes remain green;
- a deliberately removed required QLF makes the build/test fail clearly;
- identical inputs produce byte-identical data/manifest outputs.

## 7. Phase 3 — minimal compiled package/plugin profile

The generic Emscripten package set currently brings in 13 packages and links
plugins for CHR, CLP(Q/R), HTTP, semweb/RDF, SGML, PCRE, NLP, PlUnit, zlib, and
others. The browser API needs a much narrower runtime.

Add a named package profile rather than changing the native build or SWI's
general WASM default. After Phase 1, derive the exact requirement from tests;
the expected irreducible external surface is approximately:

- JSON foreign library;
- URI support required by `library(wasm)`;
- possibly the small `files`/`readutil` subset if the final browser load root
  still needs it;
- core zlib required by SWI itself, distinct from the optional Prolog zlib
  plugin.

Package-level selection is not sufficient for `clib`, whose CMake file builds
multiple plugins. Add an explicit WASM clib-plugin list or equivalent option so
the profile does not compile/link unrelated crypto, memfile, stream, and hash
plugins.

Acceptance:

- the WASM is below **1.5 MB raw** and **475 KB Brotli** on the current pin;
- the manifest and third-party notices enumerate the actual reduced link
  surface—removing packages must remove stale notice claims without removing
  licences still required by SWI core;
- `make test`, type/golden checks, the complete browser suite, error probes, and
  cancellation probes all pass;
- native artifacts and CLI behaviour are byte/semantics unchanged where their
  existing contracts require it.

## 8. Phase 4 — remove avoidable bootstrap requests

Once the dominant size work is done, address request overhead in measured
order:

1. **QLF duplicate fetch/probe.** Determine why `Prolog.consult(URL)` produces
   the known aborted duplicate. Preferred solutions are either one explicit
   fetch into MEMFS followed by a local consult, or shipping the app QLF in the
   preload image. Preserve content hashing and provenance either way.
2. **Manifest round trip.** The synchronous worker-side XHR is safe off-main-
   thread but serialises startup. Consider a versioned worker/bootstrap with
   stamped artifact names only if measurement shows meaningful RTT savings;
   do not weaken partial-deploy safety for one tiny request.
3. **Capabilities round trip.** Keep capabilities engine-sourced. It may be
   fused into the warmed/ready response if that preserves the contract and
   removes a measurable serial step.

Do not merge these changes on request-count aesthetics alone. Record timings on
a non-local network profile.

## 9. Phase 5 — configurable spare-worker policy

The generic SDK's eager spare is intentional: terminating a CPU-bound solve and
promoting a warm worker gives prompt recovery. Preserve that as the default.

Add an explicit policy such as:

- `eager` — current behaviour and default;
- `idle` — create the spare after the active engine's first successful request,
  using `requestIdleCallback` where available and a bounded timer fallback;
- `lazy` — create a replacement only when cancellation or a queued request
  requires one.

Requirements:

- no replacement worker is started after cancellation when the instance is
  being disposed and no queued work exists;
- queued work still proceeds correctly under every policy;
- in-flight cancellation remains a prompt `CancelledError`;
- recovery latency is measured and documented for each policy;
- types and SDK docs expose the option without consumer knowledge of workers.

This is primarily a memory/CPU and badly-configured-cache win. With correct
immutable caching, the measured cold-network difference between current eager
and single-worker startup was only about 58 KB Brotli; without cacheability it
was close to half the payload.

## 10. Phase 6 — compression-format gate

The current build uses Emscripten `LZ4=1` inside the data package and then relies
on HTTP gzip/Brotli outside it. On the untrimmed preload, a temporary no-LZ4
link saved about 276 KB Brotli with statistically similar local readiness
(101.1 ms versus 102.0 ms median), at the cost of a larger raw response and a
less lazy in-memory filesystem.

Re-run this comparison only after the preload is minimal. Adopt no-LZ4 only if:

- Brotli and gzip both win materially for the final profile;
- raw/uncompressed delivery is prevented or remains acceptable;
- startup and memory on lower-powered devices do not regress;
- hosting requirements explicitly guarantee content encoding.

This is a secondary gate, not a substitute for removing unused files.

## 11. Explicitly rejected low-yield paths

The feasibility pass measured two tempting build tweaks:

- post-link `wasm-opt -Oz` saved about **627 bytes Brotli**;
- reducing the exported C API from 225 symbols to the 71 bound by `prolog.js`
  saved about **5 KB Brotli** and still required the generic bridge.

Do not pursue either unless a future toolchain changes the result. A full
compile with `-Oz` is a separate speed/size experiment and must be benchmarked
against heavy arrange/fill workloads; the production `-O3` choice remains.

## 12. Verification and release gate

Every phase that changes artifacts must run, in order:

1. `make test` and existing native deterministic/golden checks;
2. SDK type/golden checks;
3. `wasm/test/headless.mjs`;
4. error, yield, cancellation, and blocked/failure probes;
5. the new payload/request/startup benchmark;
6. a from-clean `wasm/build/build-wasm.sh` run with pin/provenance guards;
7. notice and manifest validation against the actual link line;
8. consumer validation in ned.sh after re-copying the complete hashed set,
   worker, manifest, SDK facade/types, and notices.

Release target on the current SWI/toolchain pin:

| Measure | Current | Target |
|---|---:|---:|
| Engine trio, Brotli | 1.68 MB | ≤700 KB |
| Correct eager cold request model, Brotli | 1.80 MB | ≤800 KB |
| Median fresh-Chrome readiness, reference desktop | 109 ms | ≤90 ms |
| Browser regression | green | green |
| Native/value parity | green | green |

If size targets pass but lower-powered startup or memory regresses, keep the
smaller preload and back out the responsible package/compression step rather
than treating transfer bytes as the only performance metric.

## 13. Suggested landing sequence

Keep commits independently measurable and reversible:

1. benchmark harness + checked-in baseline;
2. browser-specific Prolog load root and QLF parity;
3. configurable preload profile;
4. minimal compiled package/plugin profile;
5. QLF/bootstrap request cleanup;
6. SDK spare policy;
7. final LZ4 decision;
8. provenance/notices/docs sweep and consumer re-copy.

Do not update ned.sh's presentation on an intermediate artifact. First publish
one fully verified hashed bundle; the consumer plan then measures the complete
page and decides whether the loading placeholder can be removed.
