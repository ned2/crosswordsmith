# F-H2 — bitset counting via artifact v2 (2026-07-05)

Phase 2 final experiment of the fill performance campaign
(`docs/plans/fill-perf-campaign.md`). Branch `experiment/f-h2-bitset-count`
(base 7181337, post-F-L2). Engine file changed: `prolog/crosswordsmith/fill.pl`
only (+ its plunit `tests/fill.plt`). No goldens / baseline / history / identity
manifest / fixtures written; no `--record`. Instruments live branch-only under
`benchmarks/probe_fh2/` (P-F1/gate-probe precedent). Machine: host `lana
x86_64`, SWI-Prolog 10.1.10, GMP 6. Inference counts deterministic; wall
host-specific (medians reported, 3x reproduced).

F-H2 replaces the MRV counting kernel — `ord_intersection` chains + `length/2`
over the ordset index sets — with bignum bitset `/\` chains + `popcount`, shipped
as precomputed masks inside a **schema-v2 index artifact** (the gate probe's
binding condition: masks built at `--save-index` time only, never at load/fill).
Enumeration (`candidates/4`) stays on ordsets everywhere (bit iteration is ~7x
slower than an ordset walk — scope cut, plan §143). The raw text-dict path is
untouched and stays inference-identical.

---

## PHASE A — re-attribution gate (measure first, then decide)

F-H2's economics were set by P-F1 (counting = 59-90% of `search_inf`) but two
accepted experiments moved them: F-H1 collapsed the recount VOLUME (≤5 crossing
recounts/placement instead of a full per-node recount) and F-L2 removed the load
tax. So the gate had to be re-measured on the CURRENT engine before building.

**Method.** Artifact-consuming mode (v2 artifacts built with `--save-index`
first), post-F-H1 engine. A verbatim copy of the search hot path routes the two
count/materialize seams through `call_time/2` accumulators (the engine's own
`candidate_count/4` and `candidates/4`, so the tree is identical — verified
term-identical filled grid). Denominator = the CLEAN fill-phase wall (median of
5 WARM reps). The profiler is unusable here: P-F1 §D showed LCO buries the
`ord_intersection` self-time and surfaces `length/2`'s `$skip_list` — so the wall
attribution needs call-site region timing, not `profile/2`.
(`benchmarks/probe_fh2/phase_a.pl`.)

**Attribution — 3 rungs × {counting kernel, materialization, other} × {wall, inf}
(median of 5 warm reps; inf counts exact/deterministic):**

| rung | fill wall | search_inf | counting kernel `candidate_count/4` | materialization `candidates/4` | other |
|------|----------:|-----------:|:-----------------------------------:|:------------------------------:|:-----:|
| g09_full | 0.55-0.58s | 14,869,949 | **47.3% wall** / 43.4% inf | 54.0% wall / 43.2% inf | ~0% wall / 13.4% inf |
| g17_50k  | 0.47-0.52s | 13,772,319 | **44.1% wall** / 44.9% inf | 21.5% wall / 19.4% inf | 34.4% wall / 35.7% inf |
| g21_full | 0.072-0.082s | 1,851,286 | **37.8% wall** / 47.3% inf | 41.3% wall / 27.7% inf | 20.9% wall / 25.0% inf |

- The `search_inf` here (fill_search-only) is ~0.003% below the baseline
  `fill_attempt` figures (the call_with_inference_limit + emit wrapper); shares
  are computed against the probe's own consistent total.
- The ~0% / negative "other" wall on g09 is the ~1-2% instrumentation inflation
  of the two `call_time`-measured numerators (they sum to slightly more than the
  clean denominator); harmless for the gate. By inference the true "other" is
  13-25%.

**candidate_count sub-decomposition (what F-H2 actually replaces).** F-H2 swaps
the `index_intersection` + final `length` inside `candidate_count`; it does NOT
remove `slot_bucket`'s `Bound` findall (the mask still needs the bound positions).
Median of 5 warm reps:

| rung | intersection (replaced) | count length (replaced) | Bound findall (survives) | **F-H2-replaceable share of candidate_count** |
|------|---:|---:|---:|:---:|
| g09_full | 0.214s | 0.040s | 0.069s | **78.7%** |
| g17_50k  | 0.180s | 0.024s | 0.042s | **82.9%** |
| g21_full | 0.023s | 0.004s | 0.003s | **89.8%** |

(The `count length` sub-part is ~-1 inf/call — the inference-blind `$skip_list`
walk P-F1 §D flagged — yet real wall, e.g. 0.040s on g09. `popcount` deletes it too.)

**GATE DECISION (pre-registered 20% rule).** The counting kernel is **37.8% /
44.1% / 47.3%** of fill-phase wall — ≥ 20% on ALL THREE, not borderline. **GATE
PASSES → proceed to Phase B.** Post-F-H1, counting no longer *dominates*
(materialization is now 21-54%, other 21-34%), but it is still the single
largest F-H2-addressable slice; addressable fill-phase wall ≈ counting-share ×
replaceable-share ≈ 34-37%.

---

## PHASE B — build

### Design

- **Threaded counting context `Masks`** through the search: `none` (raw text /
  any non-v2 path → ordset kernel, byte-identical to pre-F-H2) or
  `masks(MaskAssoc)` (a v2 artifact → bignum kernel). Threading a ground
  atom/compound through the recursion is **inference-neutral** (measured: a
  2-clause dispatch on the threaded arg over 100k iterations = delta 0 inf; SWI
  clause indexing has no penalty), so the `none` path is inference-identical.
- **`candidate_count/5` dispatch.** `candidate_count/4` (ordset) is retained as
  the reference for the white-box tests (P3) and `select_mrv/6`. `/5` has a
  `none` clause (old body verbatim, with a leading cut for determinism — the
  reference selector must leave no choicepoint, P13) and a `masks(M)` clause
  (bignum). The two Masks values are disjoint (atom vs compound), so the cut
  changes no answer and is not a logical inference.
- **Bignum kernel** mirrors `index_intersection` exactly: seed with the first
  bound cell's mask, `/\`-fold the rest, `popcount`. An absent key is mask 0
  (mirrors `index_set`'s `S = []` → dead cell zeroes the chain → count 0). The
  `all` branch (no cell bound) still counts the whole length bucket — semantics
  preserved.
- **Artifact schema v2.** `fill_index_format_version 1 → 2`. Masks are a new Meta
  key `masks(MaskAssoc)` (assoc `k(Len,Pos,Char) → bignum`, bit i set iff bucket
  index i is in Index's ordset). Derived FROM the ordsets at `--save-index` time
  (`build_masks/2`), so `popcount(chain) == |ord_intersection(chain)|` by
  construction. Built ONLY inside `fill_save_index` (amortized). The v1→v2 bump
  is a hard break: a v1 (or any Version ≠ 2) artifact is refused by the existing
  loader gate. `fill_load_index/5` returns the masks; `/4` is kept (discards
  them) for the roundtrip test.
- **`fill_attempt/8` (the gated bench seam) stays a pristine ordset entry**
  (body edited only to pass `none` to `fill_search/5`); the product artifact path
  uses a separate `fill_attempt_masked/9`. So the ratchet's measured predicate is
  byte-identical.

Diff: `fill.pl` +179/−31 (net +148); `tests/fill.plt` +33/−… (one v2 fix + the
new mask-count test).

**Variants considered.** (a) Duplicating the whole masks search driver (safest
for byte-identity but ~4 redundant predicates) — REJECTED once the threading was
proven inference-neutral. (b) Shape-dispatching on the Index term
(`masked(O,M)`) — REJECTED: it pollutes the shared `index_set`/`candidates`
paths with a dispatch, perturbing raw inferences. (c) Storing per-length bucket
sizes to also kill the `all`-branch `length` walk — dropped as unnecessary
(root-only post-F-H1); the `all` branch keeps `length(bucket)`.

### Correctness / equivalence

- **New plunit test** `fill_index_v2_mask_count_matches_ordset` (10k dict): for a
  broad deterministic pattern sample (lens 3-6, single- and two-cell, incl.
  rare Q/Z dead cells), `mask_count == length(index_intersection)`. GREEN.
- **Raw identity oracle 11/11** (`check_fill_identity.sh`) — raw path untouched.
- **Artifact identity oracle 11/11** (`check_fill_identity_artifact.sh`) — this
  script builds v2 artifacts via `--save-index` (now mask-bearing) and consumes
  them via `--index` → `fill_solve_index` → **bignum kernel active** (confirmed:
  masks loaded, 5,787 keys), yet every rung's stdout is byte-identical to the
  same pinned `fill_identity.sha256`. This is the end-to-end masks-active
  equivalence gate.
- **Baseline `--heavy`: +0.00% on BOTH `search_inf` and `load_inf`, all 11
  rungs** — the raw path is inference-identical (wall/RSS deltas are noise).
- **Full suite** `run_tests.sh`: 213 plunit passed / 0 failed (+1 new test),
  goldens + CLI/stderr contracts all green.

### Result — fill-phase wall (the honest artifact-mode metric)

Controlled kernel ablation (`phase_b.pl`): load a v2 artifact once, then time the
fill phase with the ordset kernel (`none` = v1 behavior) vs the bignum kernel
(`masks`) on the SAME dict + fresh slots. Median of 9 warm reps; 3x reproduced,
counts exact, equivalence IDENTICAL every run:

| rung | ordset fill wall (v1) | bignum fill wall (v2) | **fill-phase wall win** | search_inf ordset→bignum |
|------|----------------------:|----------------------:|:-----------------------:|:------------------------:|
| g09_full | 0.548-0.556s | 0.418-0.428s | **23.1-23.7%** | 14,869,949 → 9,244,437 (1.61x) |
| g17_50k  | 0.458-0.462s | 0.324-0.329s | **28.7-29.9%** | 13,772,319 → 8,102,307 (1.70x) |
| g21_full | 0.069-0.073s | 0.049-0.051s | **29.2-30.7%** | 1,851,286 → 1,006,317 (1.84x) |

- **The win is wall, not inferences.** `search_inf` collapses only 1.6-1.8x (not
  the 10x+ a pure-kernel reading suggests) because bignum counting is
  inference-blind (~few inf/count) while the surviving Bound findall +
  materialization inferences stay. `search_inf` is NOT the metric here (standing
  ledger note); the wall win (23-29%) is, and it tracks the Phase-A prediction
  (~34-37% addressable × the kernel not being infinitely fast).

### Result — build cost + artifact size (the cost side)

In-process, warm (`build_probe.pl`); real v1 artifacts rebuilt on the base engine
for the paired comparison:

| scale | v1 size | v2 size | **Δ size** | mask-construction tax (build_masks) | full v2 save |
|------:|--------:|--------:|:----------:|:-----------------------------------:|:------------:|
| 10k   | 918,301 | 1,218,577 | **+32.7%** | 0.020s | 0.137s |
| 50k   | 4,292,988 | 5,654,080 | **+31.7%** | 0.170s | 0.758s |
| 172k  | 14,602,276 | 19,279,570 | **+32.0%** | 0.759s | 3.662s |

- **The construction tax lands EXACTLY on the gate probe's prediction** (§C:
  0.021/0.167/0.922s native at 10k/50k/172k; measured 0.020/0.170/0.759s) — and
  it is paid once, at `--save-index`, off the latency path. That part of the gate
  probe's binding condition holds.
- **But the masks inflate the artifact by +32% at every scale**, and a bigger
  artifact costs more to `fast_read` on EVERY load: 172k 0.101s → 0.132s
  (+31%), 50k 0.031s → 0.040s, 10k 0.007s → 0.009s.

### Result — end-to-end CLI wall (v1 vs v2 artifact, paired, 7-sample median)

| rung | v1 CLI (ordset) | v2 CLI (bignum) | **end-to-end delta** |
|------|----------------:|----------------:|:--------------------:|
| g09_full (172k art) | 910ms | 709ms | **−22%** (search-dominated → clear win) |
| g17_50k (50k art)   | 710ms | 611ms | **−14%** (win; small artifact) |
| g21_full (172k art) | 510ms | 611ms | **+20% SLOWER** (load-dominated + tiny fill) [CORRECTED below: true net ≈ +2-3%; CLI-timer quantization] |

**The size tax is the story.** On g21_full the fill phase is ~0.07s of a ~0.5s
CLI; the bignum kernel saves ~0.02s but the +32% mask-bearing artifact costs
more on every load (fast_read +0.031s at 172k). Net end-to-end LOSS — though
the follow-up's paired re-measurement corrects the magnitude to ~+10-16ms
(~+2-3%), not +20% (see "Measurement correction" in the follow-up section: this
host's CLI wall is ~100ms-quantized and g21's true delta straddles a quantum
boundary). On search-dominated rungs (g09/g17_50k) the fill saving exceeds the
load tax → net win.

---

## Risks / anomalies / the key finding

1. **The gate probe's "artifact masks = pure win everywhere" is REFUTED for
   load-dominated rungs.** Construction amortizes to build time (confirmed), but
   the mask BYTES do not: every `--index` load pays for a +32% artifact
   (fast_read +0.031s at 172k). So end-to-end artifact-mode wall REGRESSES on
   load-dominated rungs (the ladder majority: the small-grid `*_full` rungs all
   load the 172k artifact and have sub-0.1s fill phases) and WINS only on the
   search-dominated rungs (g09_full, g17_50k, and g15-class). Masks are amortized
   in COMPUTE, not in LOAD SIZE — the correction this experiment supplies.
   (Magnitude corrected in the follow-up section: the g21 regression is
   ~+10-16ms ≈ +2-3%, not +20% — CLI-timer quantization inflated the original
   reading; direction unchanged.)
2. **RSS.** g21_full's F-H1 RSS anomaly is unrelated; the bignum masks add ~4.7MB
   resident to the loaded artifact at 172k (reporting-only, never gated).
3. **No inference regression anywhere** (raw +0.00%); the artifact-mode
   `search_inf` "collapse" (1.6-1.8x) is the inference-blind measurement blind
   spot, reported alongside the wall as required, never presented as the win.

## Verdict recommendation

**ACCEPT-WITH-CONDITION** (the code is correct and the win is real; the *default*
shipping policy is the open question for the orchestrator).

- The bignum counting kernel is correct (equivalence-proven end-to-end, new
  plunit + both identity oracles + baseline +0.00%), and delivers a genuine
  **23-29% fill-phase wall** win — exactly the search-dominated / WASM-interactive
  regime the gate probe pre-authorized as F-H2's unconditional win.
- BUT shipping masks **unconditionally** in the v2 artifact regresses end-to-end
  artifact-mode wall on load-dominated rungs (+20% on g21_full) because the +32%
  size taxes every load. On the ladder's load-dominated majority this is a net
  end-to-end LOSS.
- **Recommendation:** accept the kernel + v2 schema, but **do not make masks the
  default** — gate them behind an opt-in (`--save-index --with-masks`) or the
  WASM/search-dominated deployment, so the common load-dominated path keeps the
  smaller v1 artifact. That converts F-H2 into a pure win where it helps with no
  regression elsewhere. If the orchestrator gates purely on end-to-end
  artifact-mode wall across all 11 rungs with masks default-on, the honest call
  is HOLD/REJECT-as-default; the fill-phase kernel itself is sound and worth
  keeping.

### Draft ledger entry (superseded by the follow-up variant below)

> ### F-H2 — bitset counting via artifact v2 — ACCEPT-WITH-CONDITION (masks not default)
>
> Bignum `/\`+popcount counting kernel shipped as precomputed masks in a schema-v2
> index artifact (masks built at `--save-index` only; enumeration stays on
> ordsets). Raw path threads `Masks=none` → inference-identical (+0.00% both
> layers, 11/11); v2 artifacts thread `masks(_)` → bignum kernel.
> **Phase-A gate (re-measured post-F-H1): counting kernel 37.8/44.1/47.3% of
> fill-phase wall on g21/g17_50k/g09 — PASS (≥20%).** Built. **Fill-phase wall
> −23-29%** (g09/g17_50k/g21, equivalence IDENTICAL, 3x reproduced; `search_inf`
> −1.6-1.8x, NOT the metric). Correctness: new `mask_count==ordset` plunit (10k),
> raw + artifact identity 11/11 (masks active), baseline +0.00% both layers.
> **Cost: masks inflate the artifact +32% (mask-construction tax 0.76s at 172k —
> matches the gate probe; amortized to build). The size is NOT amortized: every
> load pays ~+0.09s at 172k → end-to-end CLI wins search-dominated (g09 −22%,
> g17_50k −14%) but REGRESSES load-dominated (g21 +20%).** The gate probe's
> "artifact masks = pure win everywhere" holds for construction, fails for load
> size. Recommend masks opt-in / WASM-scoped, not default-on.

---

## FOLLOW-UP VARIANT — masks optional within v2 (adjudicated; strict bar)

The orchestrator adopted ACCEPT-WITH-CONDITION and authorized exactly one
follow-up: make masks OPTIONAL within schema v2 — default build has no masks
(no size tax, no regression anywhere), `--masks` opts in for search-heavy /
WASM deployments. Pre-declared bar: (1) default v2 == v1 in size/load/CLI wall,
(2) masks mode retains the fill-phase win, (3) raw +0.00% both layers,
(4) identity 11/11 in three modes, (5) suite green + a no-masks-default test,
(6) 3x reproduction. All six met; evidence below.

### Change

- `fill_save_index/3` (new export): `masks(true)` in Options embeds the masks
  Meta key; `fill_save_index/2` = `/3` with `[]` (default: NO masks). Masks are
  still derived from the ordsets, still build-time only.
- Loader unchanged in shape: `fill_load_index/5` already returns `masks(_)` iff
  the key is present, else `none` → ordset kernel. v2-without-masks is a valid
  artifact; v1 files still refuse with the rebuild message (policy unchanged).
- CLI: new boolean `--masks` (help names the search-heavy/WASM use case and the
  size/load cost), passed as `masks(Bool)` to `fill_save_index/3`. House style
  (optparse boolean, longflags, same spec list).
- Tests: the mask-equivalence test now builds with `masks(true)`; NEW test
  `fill_index_v2_default_no_masks_counts_via_ordsets` pins the default shape
  (Masks == none; /5 `none` kernel agrees with the /4 reference on both
  branches; artifact-mode fill completes). 214 passed / 0 failed.
- `check_fill_identity_artifact.sh`: optional `FILL_SAVE_INDEX_FLAGS` env is
  word-split into the build command (empty = default no-masks build; `--masks`
  = masks-active build). Default invocation unchanged.

### Bar 1 — default v2 == v1 (no regression on the default path): **MET**

- **Size: byte-identical** to v1 at all three scales (CLI builds):
  918,301 / 4,292,988 / 14,602,276 bytes (10k/50k/172k) — the default Meta has
  exactly the v1 keys; only the Version integer differs.
- **Load wall (fast_read, median of 5, 3 batches): within noise** —
  172k: v1 0.1020-0.1069s vs default-v2 0.1013-0.1086s; 50k: 0.0308-0.0329 vs
  0.0306-0.0325; 10k: 0.0068-0.0070 vs 0.0068-0.0072 (deltas ≤2%, both signs).
- **End-to-end CLI wall (paired, interleaved batches, 7-sample medians, 3
  batches, exit + output checked every sample): identical to ≤1ms** —
  g21_full: v1 511/511/511ms vs default-v2 511/511/510ms;
  g11_full: v1 412/411/413ms vs default-v2 412/412/413ms. Zero failures.
- Default build time ≈ v1 build time (2.6s vs 2.4s at 172k in-process, no mask
  step; the delta is run-to-run load_dict variance, the mask step is absent).

### Bar 2 — masks mode retains the fill-phase win: **MET**

Kernel ablation re-run 3x on artifacts built by the NEW `--masks` path:
fill-phase wall −24.1/−24.7/−25.9% (g09_full), −29.6/−30.5/−32.7% (g17_50k),
−29.2/−29.5/−30.4% (g21_full); `search_inf` EXACTLY reproduced every run
(9,244,437 / 8,102,307 / 1,006,317); equivalence IDENTICAL every run. Nothing
lost in the refactor.

### Bar 3 — raw path: **MET**

`check_fill_baseline.pl --heavy` after the follow-up: exit 0, PASS, every rung
+0.00% on BOTH `search_inf` and `load_inf` (the only non-zero rows are the
informational cmd_wall/rss lines, never gated).

### Bar 4 — identity 11/11 in THREE modes: **MET**

1. **Raw**: `benchmarks/check_fill_identity.sh` — 11/11 OK.
2. **Artifact-default (no masks)**: `benchmarks/check_fill_identity_artifact.sh`
   (no env) — builds default v2, ordset kernel — 11/11 OK.
3. **Artifact-masks**: `FILL_SAVE_INDEX_FLAGS=--masks
   benchmarks/check_fill_identity_artifact.sh` — builds masks v2, bignum kernel
   active — 11/11 OK. All against the same pinned `fill_identity.sha256`.

### Bar 5 — suite green: **MET** (214 passed / 0 failed; goldens + CLI contracts OK)

### Bar 6 — 3x reproduction: **MET**

Ablation 3x (counts byte-exact, above); paired CLI medians 3 interleaved
batches per rung (bar 1 and the masks-mode table below); fast_read 3 batches.

### Measurement correction (honesty note on the earlier +20% figure)

The paired re-measurement exposed ~100ms quantization in this host's CLI wall
(samples cluster tightly at ~510/611/712/813/913ms; g21_full v1 is bimodal
across the 510/611 boundary between sessions). Contemporaneous interleaved
pairing (3 batches, 7 samples):

| rung | v1 CLI | masks-v2 CLI | net |
|------|-------:|-------------:|:----|
| g09_full | 913/913/912ms | 813/812/813ms | **−100ms (−11%)** win |
| g17_50k  | 713/712/712ms | 613/613/612ms | **−100ms (−14%)** win |
| g21_full | 611/613/612ms (bimodal 510-613) | 613/612/612ms | **~equal at quantum resolution** |

The g21_full masks-mode net is truly **≈ +10-16ms (~+2-3%)** — from in-process
components: fast_read +0.031s (172k masks) minus fill-phase −0.021s — not the
+20% the original session's quantized CLI timer suggested (v1 landed in the
510ms quantum, v2 in the 611ms one). The DIRECTION stands (masks are a net
end-to-end loss on load-dominated rungs — the size tax exceeds the tiny fill
win), the magnitude is corrected on the record. Fine-grained deltas in this doc
rest on the in-process measurements (fast_read + ablation), not the CLI quanta.
The default-off decision is unaffected either way: default v2 is exactly v1.

### Revised verdict + draft ledger entry

**ACCEPT (kernel + optional masks).** The default path is provably untaxed
(byte-identical size, identical CLI wall, +0.00% inference); the opt-in delivers
a reproducible 24-33% fill-phase wall win where fills are search-heavy.

> ### F-H2 — bitset counting via artifact v2 (masks opt-in) — ACCEPTED
>
> Bignum `/\`+popcount counting kernel, shipped as OPTIONAL precomputed masks in
> the schema-v2 index artifact: `fill --save-index` (default) emits v2 WITHOUT
> masks — byte-identical size to v1 (918,301/4,292,988/14,602,276 at
> 10k/50k/172k), load + end-to-end CLI wall within noise of v1 (g21_full/g11_full
> medians equal to ≤1ms, paired interleaved) — and `--masks` opts in for
> search-heavy/WASM fills. Counting dispatches on the loaded context: masks(_) →
> `/\`+popcount, none → the unchanged ordset kernel; enumeration stays on
> ordsets. Phase-A re-attribution gate (post-F-H1): counting = 37.8-47.3% of
> fill-phase wall, PASS ≥20% on all three hard rungs. Masks mode: fill-phase
> wall −24-33% (3x, counts exact: g09 14.87M→9.24M, g17_50k 13.77M→8.10M, g21
> 1.85M→1.01M; inference collapse is measurement-blind, wall is the metric);
> end-to-end −11%/−14% on g09_full/g17_50k, ~+2-3% on load-dominated g21_full
> (the +32% artifact-size load tax; original +20% figure corrected — CLI timer
> quantization). Mask-construction tax 0.76-0.78s at 172k, build-time only
> (matches the gate probe). Raw path +0.00% BOTH layers, 11/11; identity 11/11
> in raw, artifact-default, and artifact-masks modes; 214 tests green incl.
> mask==ordset equivalence + default-no-masks pins. Gate-probe correction
> stands on the ledger: artifact masks amortize CONSTRUCTION, not LOAD SIZE —
> hence opt-in, never default-on.
