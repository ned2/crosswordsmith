# F-H2 gate probe — bignum AND+popcount under WASM (2026-07-05)

MEASUREMENT ONLY. Gates whether F-H2 (replace the MRV counting kernel —
`ord_intersection` chains + `length/2` over ordset index sets, fill.pl:148-184 —
with bignum bitset `A /\ B` + `popcount`) is ever built. Everything here lives on
`probe/f-h2-wasm-bignum`; `prolog/crosswordsmith/fill.pl` is untouched. Scripts:
`benchmarks/probe_fh2/{kernel_bench.pl,index_stats.pl}`. Logs teed to
`/tmp/claude-1000/f-h2-probe/`.

The plan flags this as **"the one change most likely to be a native win and a
WASM loss"** (fill-perf-campaign.md:143-154). It is not. Under WASM/LibBF it is a
win too. The real constraint is *load-time construction*, not the kernel.

## Runtimes / bignum backend

| runtime | how run | `bounded` | bignum backend | provenance |
|---------|---------|-----------|----------------|------------|
| native | `swipl` 10.1.10 (pinned) | false | **GMP 6** (`gmp_version=6`) | host x86-64 |
| wasm | `node ~/src/swipl-devel/build.wasm/src/swipl.js` (NODERAWFS) | false | **LibBF** (`gmp_version=none`, `USE_GMP=OFF`) | swipl-devel `aa6289399`, **post-setjmp-fix** (`7b4b138ba` is ancestor), reports 10.1.10 |

Both compute unbounded integers correctly (`1<<200` exact on both). The WASM
build is the spike's (577f06a recipe) — read-only reused, not rebuilt.

## A. Kernel micro-bench — headline "same-N" view

Per-op cost of one pairwise combine + size: `ord_intersection(A,B,C),length(C,_)`
vs `X is A/\B, popcount(X)`. Operand size N = element count (ordset) = universe
bits with N bits set (bignum), i.e. the task's {10^2, 10^3, 10^4, ~2.4·10^4}.

| N | nat ord µs/op | nat big µs/op | **nat ratio** | wasm ord µs/op | wasm big µs/op | **wasm ratio** | ord inf/op | big inf/op |
|------:|------:|------:|------:|------:|------:|------:|------:|------:|
| 100   | 10.95 | 0.232 | 47.2x | 42.54 | 0.727 | **58.5x** | 405 | 3 |
| 1000  | 115.20 | 0.233 | 494.3x | 411.52 | 0.940 | **437.8x** | 4005 | 3 |
| 10000 | 1111.79 | 0.413 | 2691x | 4077.60 | 1.800 | **2265x** | 40005 | 3 |
| 24000 | 2608.42 | 0.672 | 3881x | 10303.67 | 2.000 | **5152x** | 96005 | 3 |

- **Red-team's 250-600x native anchor: REPRODUCED** at the 1000-element scale
  (494x native). It is *size-dependent* — 47x at 100, ~3900x at 24000 elements —
  because `ord_intersection` is O(|A|+|B|) Prolog-list traversal while the bignum
  op is O(bits/64) limb-wise C/Wasm.
- **The WASM ratio does not collapse.** LibBF bignum slows 3.0-4.4x vs native
  GMP; the Prolog VM's own ordset path slows 3.6-4.0x under WASM. The two
  penalties roughly cancel, so the WASM ratio stays the same order as native
  (438-5152x) and at N=100/24000 is even *higher* than native.
- **Bignum is inference-invisible: 3 inf/op, flat, vs 405-96005 for the ordset**
  (identical native and WASM). The win is **wall-only and inference-blind** —
  exactly the keysort/length pattern P-F1 §D/§G flagged. F-H2 would make
  `search_inf` (the gated metric) plummet, but the replacement work does not
  show in inferences; the true improvement is the ord-µs → big-µs delta.
  **Wall medians must be reported alongside any F-H2 ratchet number.**

## A'. Realistic pairing — grounded in the measured set-size distribution

The "same-N" view assumes ordsets as wide as the bucket. They are not: real
index sets are **small** (see §B distribution). A per-slot count intersects sets
of size `setS` (median 7-38), but the bignum still pays the full bucket width
`bucketU`. These are the decision-relevant pairings:

| pairing (setS, bucketU) | nat ord µs | nat big µs | **nat ratio** | wasm ord µs | wasm big µs | **wasm ratio** |
|---|---:|---:|---:|---:|---:|---:|
| 50k-median (19, 8306)    | 3.24 | 0.244 | 13.3x | 11.98 | 1.26 | **9.5x** |
| 50k-q3 (94, 8306)        | 15.58 | 0.239 | 65.2x | 57.62 | 1.30 | **44.3x** |
| 172k-median (38, 28420)  | 6.31 | 0.333 | 19.0x | 42.68 | 4.72 | **9.1x** |
| 172k-q3 (247, 28420)     | 40.06 | 0.439 | 91.3x | 180.73 | 2.63 | **68.6x** |
| 172k-p99 (2649, 28420)   | 432.49 | 0.313 | 1381x | 1608.50 | 2.67 | **603x** |
| 172k-max (10253, 28420)  | 1717.03 | 0.289 | 5948x | 6219.60 | 2.00 | **3110x** |

- **Worst realistic case = median count at full-ENABLE under WASM = 9x** bignum
  win. Everything above the median is bigger. The kernel wins on *every*
  realistic count on *both* runtimes.
- The ratio *narrows* under WASM at the small end (172k-median 19.0x → 9.1x)
  because LibBF's per-op AND cost scales ~linearly with bucket-width limbs
  (~2-5µs at the 445-limb mode bucket) while GMP is near-constant (~0.3µs). It
  narrows; it never approaches 1.
- Chain caveat: a real count folds (k-1) combines for k bound cells (k≤4 on the
  short-light bench). Both representations scale linearly in k, so the ratio is
  preserved; the table's pairwise unit is the per-combine cost.

## B. Index-set distribution + bucket sizes (grounds the operand sizes)

Loaded via the engine's own `load_dict/3`; `load_inf` reproduces the P-F1 /
baseline figures exactly (172k 26,602,936; 50k 7,757,023; 10k within +2.4k warmup).

| scale | nkeys | max bucket | set q1 | med | q3 | p90 | p99 | max | mean |
|------:|------:|-----------:|----:|----:|----:|----:|----:|----:|-----:|
| 10k   | 3871 | 1691  | 2 | 7  | 28  | 66  | 177  | 597   | 23.4 |
| 50k   | 4954 | 8306  | 4 | 19 | 94  | 272 | 785  | 2918  | 91.7 |
| 172k  | 5787 | 28420 | 5 | 38 | 247 | 822 | 2649 | 10253 | 271.4 |

Mode bucket 28420 ≈ the task's "~2.4·10^4 bits". **Half of all counts intersect
sets ≤38 elements**; the 250-600x kernel regime lives only in the p99/max tail.

## C. Construction tax — building the parallel bignum index from the ordsets

The load-time cost F-H2 adds: `foldl(\I^A^A1^(A1 is A \/ (1<<I)), Set, 0, Mask)`
per (len,pos,char) key (the natural construction). `build_inf` is **identical**
native↔WASM (inference-portable). % is of the *current* `load_dict/3` cost.

| scale | load_inf | build_inf | **% load_inf** | nat build wall | nat % load_wall | wasm build wall | wasm % load_wall | wasm/nat wall |
|------:|---------:|----------:|-----:|-----:|-----:|-----:|-----:|-----:|
| 10k   | 1,605,532 | 314,106   | 19.6% | 0.021s | 11.9% | 0.075s | 18.4% | 3.6x |
| 50k   | 7,757,023 | 1,415,657 | 18.2% | 0.167s | 18.2% | 0.441s | 22.5% | 2.6x |
| 172k  | 26,602,936 | 4,771,904 | 17.9% | 0.922s | 27.1% | 2.493s | 35.9% | 2.7x |

(load_wall native 0.177/0.917/3.406s; wasm 0.407/1.958/6.950s — WASM load penalty
only ~2.0x, milder than the doc's ~3x, consistent with the setjmp fix present.)

- **NOT the "small" the plan pre-registered.** Construction is ~18% of load_inf
  and, in wall, **27% (native) / 36% (WASM) of the load** at 172k — and WASM
  makes it *proportionally worse* because LibBF's OR-on-growing-bignum is 2.7x
  native. Load already dominates end-to-end (P-F1: 58-84% on 10/11 rungs), so a
  naive load-time build gives back a real slice of *perceived* latency, more so
  in the browser.
- **Byte-identity holds**: `popcount(mask) == length(ordset)` verified at every
  scale (`count_identical`). F-H2's count equals the current count exactly →
  MRV selection identical → search tree preserved → output byte-identical.
  Ratchet-safe by construction.

## GATE VERDICT: **GO** — kernel gate cleared; one binding construction condition

**The specific gate the plan set — "is bignum AND+popcount a WASM/LibBF loss?" —
is REFUTED.** Under LibBF the kernel is 9x (worst realistic: median count, full
ENABLE) to 603-3110x (tail) faster than the ordset kernel, and 438-5152x in the
same-N view. By the pre-registered rule (stays clearly >1 = GO) this is **GO**.
The LibBF penalty is real (~3-14x vs GMP, limb-linear) but is dwarfed by the
ordset kernel's own WASM penalty; the ratio never approaches 1, never inverts.

**Binding condition — how F-H2 builds masks (this is the whole risk):**

1. **Do NOT build masks naively at load on load-dominated rungs.** Construction
   adds ~18% load_inf and 27%/36% (native/WASM) load_wall. Since search is the
   *minority* of end-to-end on most rungs while load is the majority, a
   load-time build can make F-H2 **net-negative end-to-end** on load-dominated
   rungs (it inflates the dominant term to speed the minor one) — and worse
   under WASM. On this deployment the perceived-latency verdict is CONDITIONAL.
2. **The win is unconditional in two forms:**
   - **Search-dominated regime** (mid dict × large dense grid: g17_50k / g09_full
     class): the counting win (counting = 59-90% of `search_inf`, P-F1 §B) swamps
     the construction tax → clear net win on *both* runtimes even with load-time
     build.
   - **Artifact-borne masks**: fold mask construction into the Phase-3 precomputed
     index artifact (P-F1 §E already ranks this the biggest perceived-latency
     lever). Construction amortizes to ~0 at load → F-H2 is a pure win everywhere.
3. **Sequence**: F-H2 *after* F-H1 (F-H1 cuts the count *volume* — the recount is
   79.5-85.6% of search; F-H2 cheapens each *remaining* count) and folded with
   the Phase-3 artifact decision (which pays the construction once, offline).
4. **Ratchet**: F-H2's win is **wall-only** (bignum = 3 inf/op). `search_inf`
   will drop dramatically and *overstate* the wall gain; construction *adds*
   inference to load. Gate F-H2 on wall medians, not the inference delta alone.

Net: **build the F-H2 counting kernel — the WASM fear is dead — but ship the
masks in the artifact (or scope F-H2 to search-dominated rungs), never as a naive
load-time build on the latency-critical load path.**

## Anomalies / off-brief

1. **The WASM doc's LibBF dismissal is scope-mismatched to F-H2.**
   `swi-vm-wasm-performance.md:89` calls LibBF's penalty "considerably worse
   *rational* arithmetic (irrelevant to arrange's integer arithmetic)". F-H2 is
   *integer* bignum (AND + popcount over ~445-limb numbers). I measured that
   directly — the gate's actual object — and LibBF *does* carry a per-op integer
   penalty (~3-14x GMP, limb-linear), just not enough to matter here. The doc
   reached the right conclusion (LibBF fine) for a reason that did not cover F-H2;
   this probe supplies the missing integer-bignum evidence. Flagged, not silently
   reconciled.
2. **WASM build version banner ≠ commit provenance.** The build reports 10.1.10
   (same string as the pinned native) but is from post-setjmp-fix `aa6289399`.
   The doc warns "our pinned 10.1.10 predates that fix"; the spike correctly
   built post-fix, but nothing in the banner says so. My ~2.0x WASM load penalty
   (vs the doc's feared ~3x post-fix / ~10x pre-fix) corroborates the fix is in.
3. **Inference counts ARE portable native↔WASM here** — kernel ord_inf/op
   (405/4005/40005/96005) and build_inf (314106/1415657/4771904) are bit-identical
   across runtimes; load_inf matches to warmup. The doc (line 68-75) says
   portability "is NOT certified" and asks for a ladder diff; for these
   arithmetic/list workloads it holds. Corroborates the doc's cheap-validation
   suggestion (positive, not a contradiction).
