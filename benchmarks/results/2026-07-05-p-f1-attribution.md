# P-F1 — fill attribution probe (2026-07-05)

Phase 1 instrumentation probe of the fill performance campaign
(`docs/plans/fill-perf-campaign.md`). MEASUREMENT ONLY: everything here lives
on `probe/p-f1-attribution`; `prolog/crosswordsmith/fill.pl` is untouched and
the probe modules (`benchmarks/probe_f1/`) sit beside, never inside, the
measured engine. Machine: same host as Phase 0; SWI-Prolog 10.1.10. Inference
counts deterministic; wall host-specific.

Method summary: (A) the `load_dict/3` pipeline replicated stage-by-stage, each
stage run on the previous stage's output under `call_time/2`, calling the
module-internal originals; (B/C) a counter-instrumented verbatim copy of the
search (`benchmarks/probe_f1/search_instr.pl`) with `statistics(inferences)`
delta-regions around deterministic spans and nb_ counters/trace outside all
regions, run in a coarse (v1) and a leaf (v2) region mode, calibrated
(bookkeeping constants: `d_in=2`, `t_total≈6–7`, `c_cnt=2` inf/entry);
equivalence gated on term-identical filled grids; cross-checked against
`profile/2` exact port counts on g17_50k.

## A. Load decomposition — the 26.6M-inference startup tax, itemized

Stage table (inferences; share of whole `load_dict/3`):

| stage | 10k inf | 10k % | 50k inf | 50k % | 172k inf | 172k % | 172k wall ms |
|-------|--------:|------:|--------:|------:|---------:|-------:|-------------:|
| read_file_lines (read+split) | 30,008 | 1.9 | 150,008 | 1.9 | 518,477 | 1.9 | 43 |
| normalize findall (parse loop) | 1,167,854 | 72.9 | 5,851,486 | 75.4 | 20,229,070 | **76.0** | 1,572 |
| sort/2 dedupe (C) | ~0 | 0.0 | ~0 | 0.0 | ~0 | 0.0 | 34 |
| map_list_to_pairs(length) | 29,999 | 1.9 | 149,999 | 1.9 | 518,468 | 1.9 | 34 |
| keysort LPairs (C) | ~0 | 0.0 | ~0 | 0.0 | ~0 | 0.0 | 40 |
| group_pairs_by_key (len) | 10,020 | 0.6 | 50,023 | 0.6 | 172,847 | 0.6 | 20 |
| list_to_assoc DictByLen | 74 | 0.0 | 88 | 0.0 | 94 | 0.0 | 0 |
| build_index findall triples | 201,405 | 12.6 | 1,008,689 | 13.0 | 3,486,839 | 13.1 | 340 |
| build_index keysort (C) | ~0 | 0.0 | ~0 | 0.0 | ~0 | 0.0 | **1,150** |
| build_index group_pairs | 94,523 | 5.9 | 459,242 | 5.9 | 1,576,325 | 5.9 | 239 |
| build_index ord_set per key | 54,349 | 3.4 | 69,511 | 0.9 | 81,173 | 0.3 | 26 |
| build_index list_to_assoc | 15,036 | 0.9 | 18,098 | 0.2 | 19,764 | 0.1 | 3 |
| SUM(stages) | 1,603,262 | | 7,757,138 | | 26,603,051 | | 3,500 |
| WHOLE load_dict/3 (warm) | 1,603,141 | | 7,757,017 | | 26,602,930 | | 4,717 |
| **residue (whole−sum)** | **−121 (−0.008%)** | | **−121 (−0.002%)** | | **−121 (−0.0005%)** | | +1,217 ms |

- Inference reconciliation is exact to a flat −121 inf at every scale (stage
  runs pay 12 `call_time/2` boundaries; ≤0.008%, well inside the ±1% target).
  The warm wholes reproduce the baseline `load_inf` EXACTLY (26,602,930 /
  7,757,017 / 1,603,141).
- **The pre-registered expectation (build_index > 60% of load_inf) is
  REFUTED**: build_index totals 19.4% at 172k (22.8% at 10k). The parse loop
  (`normalize findall`) is 72.9–76.0% at every scale.
- **Inside the parse loop, the cost is one lambda** (micro-decomposition, full
  ENABLE): the loop is 19.88M inf; `string_upper` 0.35M; `string_chars` +0.17M;
  the `include([C]>>char_type(C,alpha), Cs, L)` filter is **~19.2M inf /
  ~1.4s wall** — the yall lambda is re-meta-called per character (~14
  inf/char over 1.4M chars). Same filter with a named helper: 5.58M (−72%);
  with a direct one-pass recursion: 3.83M (−81%). Output list is identical in
  both variants (order-preserving, same members) — a byte-identity-safe
  product fix candidate, sized below.
- **Inference counts are blind to the C sorts**: `sort/2`/`keysort/2` are ~0
  inf but the build_index `keysort` of 1.57M triples is 1.15s wall at 172k —
  the single largest wall stage after the parse loop. The whole-vs-stages wall
  residue (+1.2s at 172k only) is GC/allocation pressure that materializes when
  the intermediates become garbage mid-run; it belongs mostly to the same
  triple-sort region. Any load verdict argued purely on inferences understates
  build_index's wall share (~35% of wall vs 19% of inf).

### Dictionary shape + per-grid wasted-length fractions

ENABLE 172,823 words: length mode is 8 (28,420 words); lens 3–5 hold only
13,511 words (7.8%). Index keys 5,787 over lens 2–28. (Full per-length tables:
`/tmp` logs + `benchmarks/probe_f1/dict_stats.pl` output; triples = words×len,
the unit of build_index findall/keysort/group work.)

Fraction of index-build work on lengths the grid CANNOT use (full ENABLE):

| grid | slot lens | % words unusable | % triples unusable | % keys unusable |
|------|-----------|-----------------:|-------------------:|----------------:|
| fill_grid_04a | [4] | 97.7 | 99.0 | 98.2 |
| fill_grid_05a | [5] | 95.0 | 97.3 | 97.8 |
| fill_grid_09a | [3,5] | 94.4 | 97.1 | 96.5 |
| fill_grid_11a/13a/15a/17a/21a | [3,4,5] | 92.2 | 96.1 | 94.7 |

(50k and 10k subsets: within ±0.3% of the same fractions.)

**Caveat (decision-grade):** these fixtures are the deliberately short-word
completing ladder. A realistic 15x15 with lights 3–15 wastes only ~2.5% of
words / ~4.6% of triples (lengths ≥16). So the slot-length-filtered index is a
~96% cut of index work on the bench ladder but only a ~5% cut on full-range
real grids; it must be sized against the product's actual grid population, not
the ladder.

### Artifact sizings (fraction of load_inf surviving, measured shares)

| load strategy | surviving stages | surviving inf @172k | @50k | @10k |
|---------------|------------------|--------------------:|-----:|-----:|
| status quo | all | 26.60M (100%) | 7.76M (100%) | 1.60M (100%) |
| in-code parse fix (named helper, no artifact) | all, cheaper parse | ~12.5M (47%, measured kernel) | ~3.6M (~46%, extrapolated) | ~0.8M (~48%, extrapolated) |
| pre-parsed WORD LIST artifact (sorted words) | len-group chain + build_index | 5.87M (22.1%) | 1.76M (22.6%) | 0.40M (25.2%) |
| pre-grouped DictByLen artifact | build_index only | 5.16M (19.4%) | 1.56M (20.0%) | 0.37M (22.8%) |
| precomputed INDEX artifact (full .qlf-style) | ~none of the pipeline | ~0 (+ deserialization, NOT measured here) | ~0 | ~0 |

Wall note: an artifact that skips build_index also deletes the 1.15s C-keysort
+ most of the 1.2s GC residue at 172k — the wall win of the full artifact is
larger than its inference share suggests (~4.7s → ~0.05s + load-artifact time).

## B. Search attribution — counting is everything, materialization is noise

Region shares of the CLEAN baseline `search_inf` (v1 coarse regions verbatim
inside; v2 leaf split, calibrated; instrumented totals alongside; equivalence
= term-identical filled grid on every run):

| component | sq04_50k | g17_50k | g09_full | g21_full |
|-----------|---------:|--------:|---------:|---------:|
| count-path ord_intersection (via select_mrv recount) | 71.85 | 46.07 | 71.48 | 62.08 |
| count-path Bound findall (slot_bucket) | 4.34 | 22.35 | 5.33 | 12.34 |
| count-path shell (length/get_assoc/maplist) | 3.32 | 17.20 | 4.38 | 9.52 |
| winner index_intersection (candidates/4) | 17.81 | 12.68 | 17.56 | 14.56 |
| winner Bound findall | 1.25 | 0.59 | 0.53 | 0.25 |
| winner materialization (maplist nth0_of walk) | 0.50 | 0.25 | 0.24 | 0.63 |
| select_mrv sort+select | 0.44 | 0.62 | 0.25 | 0.44 |
| try loop member/memberchk/unify (residual) | 0.50 | 0.24 | 0.23 | 0.18 |
| SUM | 100.00 | 100.00 | 100.00 | 100.00 |
| clean total (baseline search_inf) | 7,738,070 | 19,637,890 | 34,880,750 | 3,307,580 |
| instrumented total v1 / v2 | 7.94M / 8.22M | 19.88M / 21.80M | 35.31M / 36.64M | 3.32M / 3.51M |

Brief-mapped aggregates:

| brief component | sq04_50k | g17_50k | g09_full | g21_full |
|-----------------|---------:|--------:|---------:|---------:|
| (1) candidates/4 materialization | 0.50 | 0.25 | 0.24 | 0.63 |
| (2) counting (intersections + count shell) | 92.98 | 75.95 | 93.42 | 86.15 |
| (3) Bound-extraction findall (both paths) | 5.59 | 22.94 | 5.85 | 12.59 |
| (4) member + memberchk(Used) + unification | 0.50 | 0.24 | 0.23 | 0.18 |
| (5) everything else (sort+select) | 0.44 | 0.62 | 0.25 | 0.44 |

- **ord_intersection counting carries 59–89% of search_inf** (count path +
  winner path; per count-path intersection: 298–525 inf, tracking index-set
  size with dictionary scale). `select_mrv`'s full per-node recount
  (`r_count_all`) alone is 79.5–85.6% on all four rungs.
- **The plan's ranked gap #1 (candidate materialization) is REFUTED as a
  search cost**: 0.24–0.63%. The winner's `maplist(nth0_of)` walks are tiny
  because MRV picks nearly-decided slots (g17_50k: only 6,660 nth0_of calls
  across 5,579 nodes — ~1.2 candidates materialized per node).
- The used-word scan is noise: memberchk hit 215 of 5,794 tries on g17_50k;
  component (4) ≤0.5% everywhere. F-H4 has nothing to attack.
- V1 vs V2 cross-check: V2's `r_count_all` equals V1's plus the modeled
  nested bookkeeping to within 0.3% on every rung (the model: entries×(t_total
  ≈7) + inter-counter×2), so the leaf split is internally consistent.

## C. Backtrack / fan-out — the inversion explained

| rung | clean inf | nodes | inf/node | tries | placements | unwound | node fails | fail depth mean/max | count calls | counts/node |
|------|----------:|------:|---------:|------:|-----------:|--------:|-----------:|---------------------:|------------:|------------:|
| sq04_full | 529,275 | 100 | 5,293 | 103 | 100 | 92 | 92 | 5.4 / 7 | 280 | 2.8 |
| sq04_50k | 7,738,070 | 4,748 | 1,630 | 4,851 | 4,748 | 4,740 | 4,740 | 4.1 / 7 | 18,646 | 3.9 |
| g17_full | 2,534,409 | 327 | 7,750 | 361 | 327 | 227 | 227 | 26.4 / 99 | 21,755 | 66.5 |
| g17_50k | 19,637,890 | 5,579 | 3,520 | 5,794 | 5,579 | 5,479 | 5,479 | 50.2 / 99 | 277,773 | 49.8 |
| g09_full | 34,880,750 | 9,961 | 3,502 | 10,076 | 9,961 | 9,929 | 9,929 | 20.1 / 31 | 118,983 | 11.9 |
| g21_full | 3,307,580 | 423 | 7,819 | 536 | 423 | 269 | 269 | 102.5 / 153 | 25,796 | 61.0 |

- **The inversion is a node-count effect, confirmed**: same grid, thinner
  dictionary → sq04: 47.5x nodes at 0.31x inf/node = 14.6x total; g17: 17.1x
  nodes at 0.45x inf/node = 7.7x total. Thin-dict nodes are CHEAPER (smaller
  index ordsets to intersect, smaller buckets to length-count); there are just
  vastly more of them.
- **Refinement of the pre-registered hypothesis**: "full-ENABLE trees are
  near-backtrack-free" holds only for loose interlocks (g21_full: 423 nodes,
  64% unwound; sq04_full: 100 nodes). g09_full — full ENABLE — is the hardest
  rung BECAUSE its dense 9x9 interlock churns: 9,961 nodes, 99.7% of
  placements unwound, with 99% of nodes inside a depth-17–24 thrash band
  (deciles 5–6: 9,858 nodes). Dictionary abundance shrinks trees; it does not
  guarantee shallow ones. Cost = (interlock-driven node count) × (dict-scale ×
  slot-count per-node counting price).
- Candidate fan-out over time (BestCount per depth decile): every rung's
  median selected-slot count collapses to 0–1 after the first placements
  (MRV chases forced/dead slots); dead-end discovery is via the RECOUNT, not
  via trying words. E.g. g09_full deciles 5–6: median 0, max 12–26; g17_50k
  decile 0 median 0 (1,830 root-decile nodes re-entered after unwinds).
  Full tables in `/tmp/claude-1000/pf1_*_v1.log` (committed in-doc excerpts
  only; regenerate with `run_search.pl RUNG v1`).

## D. Profiler cross-check (g17_50k) — both metrics needed, each blind alone

`profile/2` exact port counts corroborate the counters EXACTLY: 277,773
slot_candidate_count calls (= trace-derived count_calls); 5,579
select_mrv/candidates; 33,824 index_intersection (28,248 count-path + 5,576
winner, = nodes minus 3 open-slot winners); $memberchk 5,794 calls / 215
succeeds (= tries − placements); nth0_det 6,660 (materialization walk);
$skip_list 561,225 = the modeled length/2 call census (283,352 slot_bucket +
249,525 all-path counts + 28,248 idx counts + 100 emit).

Two misattributions caught, one in each direction:

- **Profiler self-time buries the ordset walk**: ordsets isect2+isect3 show
  9/1,296 ticks (0.7%) self despite 7.3M calls (46% of inferences) — their
  frames vanish into callers under last-call optimization (some re-surface as
  `slot_bucket` self: 176 ticks). Cumulative frames blame callers; treat as
  structure only. The counters are the trustworthy inference attribution.
- **Inference counts bury length/2**: `$skip_list` (C helper of length/2) is
  369/1,296 self ticks = **28.5% of profiled wall** at ~1 inf/call — the
  `Sel == all -> length(Words, Count)` path walks whole length-buckets
  (thousands of cells) 249,525 times on g17_50k. A bucket-length cache (or
  F-H1/F-H3 carrying sizes) buys ~30% of search WALL on this rung while being
  nearly invisible in search_inf.

## E. Verdicts

**(i) Phase 2 ranking** (share attacked, on the measured rungs):

1. **F-H1 — incremental candidate counts**: attacks the per-node full recount
   (`r_count_all` = 79.5–85.6% on all four rungs; the winner-path work,
   12.7–17.8% intersections + materialization, is NOT touched by it).
   Only slots crossing the just-bound cells can change count; with max light
   length 5 that is ≤5 of the ~50–61 slots recounted per node on the big
   grids. Expected order-of-magnitude cut of the dominant share on every rung
   class; also deletes most of the length/2 wall (the hidden 28.5%). Highest
   ceiling, moderate risk (threaded state).
2. **F-H2 — bitset counting** (gated, below): attacks ord_intersection
   58.8–89.7%. Complements F-H1 (cheapens the counts F-H1 still does);
   without F-H1 it caps at ~59–90% alone.
3. **F-H3 — compound-term buckets + arg/3**: as originally scoped
   (materialization) it attacks 0.24–0.63% — near-worthless for search_inf.
   RESCOPE: its real value is O(1) `bucket_size` (kills the length/2 wall) and
   cheaper index-set representation feeding F-H2; build it as infrastructure
   inside F-H1/F-H2, not as its own experiment.
4. **F-H5 — lazy candidate generation**: attacks the same 0.24–0.63% + part
   of (3)'s winner share (~0.5%). DROP unless a later tree regime (deep
   member walks) appears; today tries ≈ placements (words are almost never
   rejected by memberchk, and member yields average ~1.02/node).
5. **F-H4 — Used as a set**: attacks ≤0.5% (215 memberchk hits in 5,794
   tries on the worst rung). DROP.

**(ii) F-H2 gate verdict: YES — run the WASM bignum probe.** ord_intersection
counting carries 58.8% (g17_50k) to 89.7% (sq04_50k) of search_inf — far above
any reasonable gate bar. The 250–600x kernel measurement now has a measured
target share; even a 10x effective speedup on the intersection component is a
~2–8x rung-level win. The WASM/LibBF AND+popcount measurement remains the
build/no-build condition per the plan.

**(iii) Phase 3 sizing + priority.** Phase 0 showed dict_load >50% of
end-to-end latency on 10/11 rungs; this probe shows where inside:

- (a) **slot-length-filtered index**: on the bench ladder saves ~96% of index
  triples (≈ build_index ≈ 19.4% of load_inf + the 1.15s keysort wall); on
  realistic full-range grids only ~5%. Worthwhile ONLY as the filtered variant
  of a rebuilt loader, not as a stand-alone experiment; do not oversell from
  ladder fixtures.
- (b) **pre-parsed word-list artifact**: load_inf → 22.1% of status quo
  (26.6M → 5.87M @172k). Also deletes the parse loop's 1.5s wall.
- (c) **precomputed index artifact**: load_inf → ~0 (+ unmeasured
  deserialization); also deletes the 1.15s keysort + ~1.2s GC wall. Biggest
  single perceived-latency lever; needs a .qlf-load cost measurement first.
- (d) **NEW, found by this probe — the include-lambda parse fix** (named
  helper in `normalize_word`): no artifact, no format, output-identical;
  load_inf 26.6M → ~12.5M (−53%), wall −~1.2s at 172k. Trivial diff, ships
  ahead of any artifact work.

**Priority: Phase 3 leads Phase 2 for perceived latency, starting with (d)
then (c)** — dict_load is 58–84% of end-to-end on 10/11 rungs and its biggest
component is a one-line-fixable lambda; search work (Phase 2) matters mainly
for the search-dominated regime (g17_50k-like: mid-size dicts × large grids)
and for the WASM budget story. Recommended order: (d) → F-H1 → F-H2 (post
WASM gate, folding F-H3's bucket-size infrastructure) → (c)/(b) artifact
decision with a measured .qlf-load probe.

## F. Soundness evidence

- **Equivalence**: every instrumented run (4 rungs × v1+v2, plus sq04_full/
  g17_full v1) printed `equivalence=IDENTICAL` — the instrumented search's
  filled grid is term-identical to the clean engine's on the same fresh slots.
- **Clean reproduction**: the real bench pathway
  (`bench_core:measure` + `fill_subjects:fill_search_sampler`) inside a probe
  process reproduces baseline search_inf EXACTLY (sq04_50k: 7,738,070).
  The probe runner's own clean printout read +1 on fill_attempt call #2 with
  the probe module loaded (a one-time deferred JIT/indexing event; call #3
  is byte-exact) — measurement-context artifact, engine unaffected.
- **Non-perturbation gate**: `check_fill_baseline.pl` on this branch — PASS
  twice: default (7 core rungs) and `--heavy` (all 11 rungs), every rung
  +0.00% on BOTH search_inf and load_inf, exit 0 both times. Probe modules
  live beside the engine; nothing on the measured path changed.
- **Load reconciliation**: flat −121 inf residue (≤0.008%) at all three
  scales; warm wholes equal baseline load_inf exactly.
- Attribution shares are computed against CLEAN totals from region deltas
  whose bookkeeping sits outside the deltas (calibrated d_in=2 subtracted on
  leaf regions); instrumented totals are never compared to baselines.

## G. Anomalies / off-brief observations

1. **Inference-blind wall costs are first-class in fill**: C-implemented
   `keysort` (1.15s at 172k load) and `length/2` (28.5% of g17_50k search
   wall) never show in the gated metrics. The campaign's inference ratchet is
   sound for tree-shape regressions but will not see wins/losses in these C
   walks; wall medians must stay reported alongside.
2. The `Sel == all -> length(Words, Count)` bucket-length walk is O(bucket)
   per open-slot count — 249,525 times on g17_50k. Cheapest independent fix:
   cache per-length bucket sizes at load (an assoc Len→Size), byte-identical
   counts.
3. `get_assoc` on the big Index assoc is a real constant: 345,256 calls
   (g17_50k) ≈ 5% of profiled wall ($btree_find_node+arg self ticks) — F-H3's
   compound-term/arg(3) idea applies to the INDEX lookup path too.
4. The load pipeline's `map_list_to_pairs(length, Words)` re-walks every word
   (0.52M inf @172k) immediately after `normalize_word` already knew each
   word's length — a pre-parsed artifact absorbs this; an in-code fix could
   emit Len-Word pairs directly.
5. GC/allocation pressure adds ~1.2s wall to the un-staged full-ENABLE load
   that no stage exhibits in isolation (whole 4.72s vs stage-sum 3.50s) —
   artifact strategies that skip the triple findall/keysort remove most of it.
