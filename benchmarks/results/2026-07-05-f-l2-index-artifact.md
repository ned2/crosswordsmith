# F-L2 — precomputed index artifact — results

**Experiment:** F-L2 (Fill campaign, Phase 3). An ADDITIVE, explicit, off-by-default
precomputed dictionary-index artifact that removes the dominant startup tax
(dictionary parse + index build) from fill's interactive path.

**Branch:** `experiment/f-l2-index-artifact` (base `f9e5a83`). Host: lana x86_64.
SWI-Prolog 10.1.10. Nothing merged; nothing recorded; raw baseline untouched.

**Hypothesis (confirmed):** the parse + `build_index` is a pure function of the
frozen dict file and dominates end-to-end latency (P-F1: dict_load 58-84% of CLI
wall on 10/11 rungs; F-L1 cut the parse, leaving the inference-blind `build_index`
keysort + GC residue). Precompute once, serialize the EXACT runtime structures
(the `DictByLen` buckets + the assoc `Index`), load them back instead of
recomputing. Loaded structures are `==`-identical to a fresh build (equivalence
trivial), so fill output is byte-identical.

---

## (a) Format comparison — A (.qlf) vs B (fastrw)

Both formats serialize the term `t(DictByLen, Index)` (production adds a Meta
keylist, ~+270 bytes). Load = warm, median of 3 after 1 warmup, in-process.
"eq" = loaded term `==` freshly-built term.

| format | scale | serialize/qcompile wall | load wall | load inf | file bytes | eq |
|--------|-------|------------------------:|----------:|---------:|-----------:|:--:|
| **B fastrw** | 10k  | 0.011s | **0.007s** |   4 |    918,104 | yes |
| **B fastrw** | 50k  | 0.055s | **0.030s** |   4 |  4,292,790 | yes |
| **B fastrw** | 172k | 0.199s | **0.099s** |   4 | 14,602,081 | yes |
| A .qlf | 10k  | 0.376s | 0.014s | 710 |    930,535 | yes |
| A .qlf | 50k  | 2.380s | 0.062s | 710 |  4,348,188 | yes |
| A .qlf | 172k | 9.295s | 0.269s | 710 | 15,480,981 | yes |

Reference: post-F-L1 **warm raw** load_dict at 172k ≈ **1.96 s** / 10.7M inf.

- fastrw loads **~20x** faster than raw at 172k (0.099 s), **~2.7x** faster than
  .qlf, in a smaller file, at ~0 inferences.
- .qlf loads ~7.3x faster than raw at 172k — still under the <0.5 s target, but
  its qcompile build is ~9.3 s (vs fastrw's 0.2 s serialize) and its file is
  larger. The .qlf author-measured "~20x" is for *code*; a data-heavy single
  fact lands at ~7x, exactly the spread this experiment was to measure.
- **No qlf clause-size / qlf limit hit** at 172k — the single fact compiled and
  loaded cleanly (sharding was not needed).

**Shipped: B (fastrw)** — the fastest candidate that meets every gate. All
numbers reported; the .qlf/WASM consideration is surfaced under Risks (f) and
the verdict (g) for orchestrator adjudication.

**Production one-off build cost** (`fill_save_index/2` = load_dict + SHA-256 +
fast_write, the honest end-to-end build): 10k 0.11 s · 50k 0.62 s · 172k 2.40 s.
Dominated by the load_dict it wraps; serialization itself is the ~0.2 s tail.

---

## (b) End-to-end CLI wall — raw (`--dict`) vs artifact (`--index`)

Full `crosswordsmith fill` process, median of 3 after 1 warmup, /usr/bin/time
`%e`. Full-ENABLE core rungs + g09_full (all use enable1.txt, 172,823 words).

| rung | raw ms | artifact ms | speedup |
|------|-------:|------------:|--------:|
| sq04_full     | 2530 | 230 | 11.00x |
| g11_full      | 2590 | 300 |  8.63x |
| g11_full_seed | 2520 | 240 | 10.50x |
| sq05_full     | 2560 | 250 | 10.24x |
| g17_full      | 2620 | 410 |  6.39x |
| g21_full      | 2850 | 500 |  5.70x |
| g13_full      | 2640 | 400 |  6.60x |
| g09_full      | 3220 | 870 |  3.70x |

The ~1.96 s warm dict-load tax is removed; the artifact-mode residual
(230-870 ms) is SWI startup + artifact load (~0.14 s) + grid + search + emit.
Search-light rungs speed up most (11x); the top search rung g09_full (14.87M
search inf) speeds up least (3.7x) because search, not load, now dominates it —
the expected signature.

**Raw-path +0.00% evidence (`check_fill_baseline.pl --heavy`):** all 11 rungs,
BOTH gated layers, delta +0.00% ok, `RESULT: PASS (no change; 0 regressions)`.
search_inf and load_inf are measured on `fill_attempt/8` and `load_dict/3` in
isolation (fill_subjects.pl), which this change does not touch; the `fill_solve`
refactor only reorganises the CLI glue around them.

```
search_inf: all 11 rungs +0.00% ok  (sq04_full 371,343 … g09_full 14,870,446)
load_inf:   all 11 rungs +0.00% ok  (full-ENABLE 10,724,707; 50k 3,164,117)
```

---

## (c) Design

### CLI surface (additive, explicit, off by default)

Two flags on the existing `fill` verb (house pattern: mutually-exclusive flag
guards, exactly like arrange's `--size`/`--max-size`). The raw path is unchanged
when neither is given.

- `fill --dict D --save-index FILE` — BUILD mode: serialize the index from `D`
  (or the bundled sample) to `FILE`, then exit. No `--grid`, no fill.
- `fill --grid G --index FILE [--seeds S]` — CONSUME mode: fill using the
  artifact instead of parsing a dict. If `--dict D` is ALSO given, verify the
  artifact's embedded SHA-256 matches `D` (else refuse).
- `--save-index` + `--index` together → clear "mutually exclusive" error.

### Artifact term (versioned + extensible)

```prolog
fill_index(Version, Meta, DictByLen, Index)
```

- `Version` — integer artifact-SCHEMA version (this build: `1`). The loader
  refuses an unrecognised version.
- `Meta` — a keyed list `[Key(Value), ...]`:
  - `dict_sha256(Hex)` — SHA-256 of the source dict bytes (matches coreutils
    `sha256sum`; staleness key)
  - `swi_version(V)` — the SWI-Prolog that built it (binary-format guard)
  - `words(N)`, `source(Path)`, `built_epoch(E)` — informational
- `DictByLen`, `Index` — the exact runtime structures, serialized verbatim.

**On-disk format:** `fast_write`/`fast_read` (library(fastrw)), single term.

### F-H2 mask extension point

The F-H2 bitset masks slot into **`Meta`** as an added `masks(...)` key under a
**Version bump** (1 → 2). This is a version bump, NOT a format break: v1 readers
never look for the key; a v2-aware reader gates on `Version >= 2` before using
it. F-H2's binding condition ("masks shipped inside the F-L2 artifact,
construction amortized to ~0") is thus satisfied by construction — the builder
would compute masks once at `fill_save_index` time and append them to Meta.

### Staleness / integrity policy (refuse, never silently rebuild)

`fill_load_index/4` checks, in order, throwing a clear `fill_index_*` error on
any failure (rendered by `prolog:error_message//1`; CLI exit 1):

1. file exists (`fill_index_missing`)
2. `fast_read` succeeds (`fill_index_unreadable` — catches a cross-SWI binary)
3. term is `fill_index/4` (`fill_index_malformed`)
4. `Version == 1` (`fill_index_version`)
5. `swi_version` matches the running SWI (`fill_index_swi`)
6. iff `--dict` given: `dict_sha256` matches that file (`fill_index_hash`)

Three independent version axes are kept separate: our schema version (arg 1),
the builder's SWI version (Meta, binary-compat), and the file's own fastrw
framing (SWI's loader). No transparent auto-caching; no staleness guesswork on
the default path.

### Product-code changes (load_dict/build_index untouched)

`fill.pl`: added exports `fill_solve_index/5`, `fill_save_index/2`; refactored
`fill_solve/4` to share `fill_prepare/5` (grid+seeds) and `fill_place_and_emit/6`
(search+emit) with the new artifact path — so both paths converge on identical
search + emit, and the raw path's behaviour is byte-for-byte preserved (proven
by the raw identity oracle 11/11 and baseline +0.00%). `load_dict/3` and
`build_index/2` are NOT modified: `fill_save_index` calls `load_dict`, and
`fill_load_index` reconstructs `==`-identical terms.

### Variants tried

- **.qlf (format A)** — measured, works at all scales, no clause-size limit;
  ~2.7x slower load, larger file, ~9.3 s build at 172k. Not shipped (fastrw
  strictly faster natively) but fully measured; retained as the WASM-aligned
  alternative for the browser milestone (see f/g).
- Format C (generated .pl, plain-consulted) — not needed; A and B both cleared
  the gates.
- Meta as a Prolog dict vs a keyed list — chose the keyed list: robust under
  `fast_write`/`fast_read` and trivially extensible (add a key = version bump).

---

## (d) Equivalence evidence

- **Round-trip lemma (production predicates), loaded == built:**
  | scale | DictByLen == | Index == | build wall | load wall | bytes |
  |-------|:---:|:---:|--:|--:|--:|
  | 10k  | yes | yes | 0.112 s | 0.010 s |    918,377 |
  | 50k  | yes | yes | 0.618 s | 0.033 s |  4,293,064 |
  | 172k | yes | yes | 2.398 s | 0.141 s | 14,602,352 |
- **CLI identity, RAW mode:** 11/11 OK vs `fill_identity.sha256` (raw path
  untouched).
- **CLI identity, ARTIFACT mode:** 11/11 OK vs the SAME
  `fill_identity.sha256` — new script
  `benchmarks/check_fill_identity_artifact.sh` (mirrors the raw oracle with
  `--index`, builds artifacts into a scratch dir, deletes on exit). **THE gate —
  passed.**
- **Tests:** `run_tests.sh` green; plunit 212 passed / 0 failed (+8 F-L2 tests:
  round-trip `==`, fill-matches-raw, hash-ok, hash-mismatch refusal, version
  refusal, swi refusal, malformed refusal, missing refusal).
- **Baseline ratchet:** `--heavy`, +0.00% both layers, all 11, PASS.

---

## (e) Files + commit + branch

- Branch `experiment/f-l2-index-artifact` (base `f9e5a83`).
- `prolog/crosswordsmith/fill.pl` — artifact predicates + shared-body refactor.
- `crosswordsmith` — `--save-index` / `--index` flags + 3-mode dispatch.
- `tests/fill.plt` — 8 F-L2 tests.
- `benchmarks/check_fill_identity_artifact.sh` — artifact-mode identity oracle (new).
- `benchmarks/results/2026-07-05-f-l2-index-artifact.md` — this doc.
- NOT touched: `fill_baseline.json`, `fill_history.jsonl`, `fill_identity.sha256`,
  `fixtures/dict/`, goldens. No `--record`.

---

## (f) Risks / anomalies

- **fastrw cross-version portability (staleness distribution).** `fast_read`'s
  binary format is SWI-version-bound (documented). The artifact is therefore a
  BUILD-TIME cache, not a distribution format across SWI versions: embedding
  `swi_version` + refusing a mismatch handles it safely (a `fill_index_swi` /
  `fill_index_unreadable` error → rebuild), but it means an artifact shipped to a
  machine on a different SWI must be rebuilt there. Fine for the native CLI
  (rebuild once locally); a consideration for any prebuilt-artifact distribution.
- **WASM / browser milestone (the ultimate consumer).** The documented WASM
  fast-path is `.qlf` + `Prolog.consult()` from a URL; fastrw's availability
  under swipl-wasm is unverified here. Shipped fastrw wins the NATIVE CLI
  decisively; the browser milestone should either (i) verify fastrw under
  swipl-wasm, or (ii) emit `.qlf` (measured here: 7.3x load win, WASM-aligned).
  Because the format is a single builder/loader predicate pair behind a versioned
  term, adding a `.qlf` emitter later is localized — the artifact schema does not
  change.
- **No RSS anomaly** (unlike F-H1's g21_full): artifact-mode cmd_rss tracks raw
  mode; the change is startup, not search.

---

## (g) Verdict recommendation + draft ledger entry

**Recommendation: ACCEPT (ship fastrw).** All gates pass: raw path +0.00% both
layers ×11; identity 11/11 in BOTH modes; equivalence `==` at 3 scales; tests
green. Delivers 3.7-11x end-to-end CLI speedup by removing the ~1.96 s dict-load
tax, exactly the dominant term P-F1 identified. Off by default, additive,
explicit; raw path byte-for-byte untouched. Orchestrator decision point: whether
to also require/emit `.qlf` for the browser milestone now, or defer to that
milestone (fastrw shipped for native; `.qlf` numbers in hand).

Draft ledger entry:

> ### F-L2 — precomputed index artifact — ACCEPTED (recommended)
> - **Where:** branch experiment/f-l2-index-artifact (base f9e5a83).
>   prolog/crosswordsmith/fill.pl (new fill_save_index/2, fill_solve_index/5;
>   fill_solve refactored to share fill_prepare/5 + fill_place_and_emit/6;
>   load_dict/build_index UNTOUCHED); crosswordsmith CLI (--save-index / --index,
>   3-mode dispatch); benchmarks/check_fill_identity_artifact.sh (new artifact
>   identity oracle); tests/fill.plt (+8); results doc
>   benchmarks/results/2026-07-05-f-l2-index-artifact.md.
> - **What + mechanism:** serialize the exact DictByLen buckets + assoc Index as
>   a single fastrw term `fill_index(Version, Meta, DictByLen, Index)` (Meta =
>   keyed list: dict SHA-256 + SWI version + provenance; the F-H2 mask extension
>   point, added under a version bump). `--save-index` builds; `--index`
>   consumes, verifying schema-version + SWI-version (+ dict hash iff --dict
>   given), REFUSING every mismatch (no silent rebuild).
> - **Result:** artifact load 1.96 s → ~0.10 s warm at 172k (~20x, ~0 inf);
>   end-to-end CLI 3.70x (g09_full, search-dominated) to 11.0x (sq04_full,
>   load-dominated). Raw path +0.00% both gated layers on all 11 rungs
>   (core+heavy); CLI identity 11/11 in BOTH raw and artifact modes vs the
>   pinned fill_identity.sha256; loaded==built (==) at 10k/50k/172k; 212 tests
>   green. Format bake-off: fastrw beat .qlf (0.099 s vs 0.269 s load, smaller
>   file, 0.2 s vs 9.3 s build) — .qlf works at 172k with no clause-size limit
>   and is retained as the WASM-aligned alternative.
> - **Ratchet:** artifact-mode load is C-level (fastrw) — its load_inf is ~4
>   (near-zero) while WALL is the honest metric; do not read the inference
>   collapse as the win. Proposed artifact-mode rungs are a PROPOSAL (see below),
>   not recorded.

---

## (h) Proposed bench-manifest rows for artifact-mode rungs (PROPOSAL ONLY)

Not recorded — the orchestrator decides. Artifact-mode load is C-level, so its
**load_inf is ~4 (near-zero) and gates nothing meaningful**; the honest ratchet
metric for artifact mode is **load WALL**, which is host-specific
(reporting-only under the current schema). search_inf is unchanged from raw mode
(same fill_attempt). Options:

1. **No new gated rungs (recommended).** Artifact mode adds no NEW deterministic
   inference count to gate (search_inf already gated via raw; artifact load_inf
   ≈ 4 is trivially constant). Track artifact load WALL in the results doc /
   history as reporting-only, as wall is treated everywhere else.
2. If a regression tripwire on the artifact path is wanted, gate the **artifact
   file byte size** per dict (a deterministic, host-independent number:
   918,377 / 4,293,064 / 14,602,352 at 10k/50k/172k) — it changes iff the
   serialized structure changes, catching an accidental format/shape drift
   without depending on wall. This would need a new column in the schema, not a
   search_inf/load_inf rung.

Suggested (if option 2 adopted) — one row per frozen dict, keyed by dict:

| artifact_rung | dict | schema_ver | file_bytes (gate) | build_wall_ms (info) | load_wall_ms (info) |
|---|---|---|---|---|---|
| art_172k | enable1.txt    | 1 | 14,602,352 | ~2400 | ~100 |
| art_50k  | enable_50k.txt | 1 |  4,293,064 |  ~620 |  ~33 |
| art_10k  | enable_10k.txt | 1 |    918,377 |  ~112 |  ~10 |
