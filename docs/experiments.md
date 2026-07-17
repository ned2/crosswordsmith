# Algorithm experiment log

Append-only record of solver-algorithm experiments: what we tried, what we
measured, and **why** we kept or rejected it. The goal is to never re-run a
dead end and to be able to write up comparable findings over time.

## How to use this log

- **Code is the source of truth for reproduction, not this log.** Every *kept*
  algorithm lives as a selectable strategy in `crossword.pl` (`--strategy`,
  see `select_word/9` and `strategies/1`). To reproduce its numbers, run it;
  do not reconstruct it from prose here.
- **This log is the narrative + index.** Each entry records a hypothesis, the
  change (by strategy id / commit), a results snapshot *with its conditions*,
  and a verdict. Rejected ideas live here as prose (often with no surviving
  code) — that is their whole point: don't try them again.
- **Numbers here are historical snapshots.** Never compare numbers across
  entries directly. For an authoritative comparison, regenerate the whole
  matrix now (`make bench-matrix`) so every cell shares one machine + build.

## Metrics

- **Inferences are the primary metric.** They are deterministic and
  machine-independent, so a logged inference count stays comparable across
  machines as long as code + fixture + config are pinned.
- **Wall time is secondary** ("is it fast enough in practice"). It is
  machine- and load-dependent; do not compare wall time across runs/machines.

## Reproducing a result batch

```sh
make bench-matrix                       # all strategies x all fixtures
swipl -q benchmarks/run_matrix.pl -- baseline mrv_capped   # subset
```

- Per-fixture config (grid, start, iterations, warmup) is read from the
  manifest `benchmarks/fixtures.pl` — the **single source of truth**, mirrored
  for humans in `fixtures/README.md`.
- Raw CSV batches are saved under `benchmarks/results/`, tagged with the SWI
  version and git commit that produced them.

> **Methodology lesson (2026-06-18).** The first benchmarking pass hardcoded a
> 17x17 grid for every fixture and reported `benchmark_20`/`benchmark_26` as
> *unsatisfiable*. They are not — they are designed for 37x37 and 49x49 grids
> and solve fine there. The "unsat" was an artifact of the wrong grid. This is
> exactly why per-fixture config must live in data (the manifest) that travels
> with results, and why the harness reads only the manifest. Don't reintroduce
> a hardcoded grid.

---

## Suite property: start-position bias

Every benchmark pins one start position (`topleft_across` for all manifest
fixtures). A sweep over all four named starts (`start_locs/1`), baseline and
`mrv_inc`, each fixture on its manifest grid, exposes how much that pin matters.
Data: `benchmarks/results/2026-06-20-start-sensitivity.csv` (reproduce with
`swipl -q benchmarks/start_sensitivity.pl`). Inferences per start (TL-A =
topleft_across, TL-D = topleft_down, TR = topright, BL = bottomleft):

| fixture | strategy | TL-A | TL-D | TR | BL |
| --- | --- | --- | --- | --- | --- |
| bundled_17 | baseline | 25.6 k | 16.9 k | 9.3 k | 9.2 k |
| bundled_17 | mrv_inc | 322 k | 305 k | **22 k** | **22 k** |
| benchmark_08 | both | ~10 k / 53 k | ~10 k / 53 k | **unsat** | **unsat** |
| benchmark_14 | both | solves | solves | **timeout** | **timeout** |
| benchmark_16_dense | both | solves (slow) | solves (slow) | **timeout** | **timeout** |
| benchmark_70_mesh | baseline | 293 k | 292 k | **timeout** | **timeout** |
| benchmark_70_mesh | mrv_inc | 5.16 M | 5.14 M | 4.84 M | 4.77 M |
| benchmark_20 | mrv_inc | 79 k | 81 k | 154 k | 158 k |
| benchmark_26 | mrv_inc | 139 k | 143 k | 274 k | 280 k |

Findings:

- **The synthetic fixtures are start-locked by construction.** `benchmark_08`
  (provably, small grid), `benchmark_14`, `benchmark_16_dense` and
  `benchmark_70_mesh` only solve from the two **top-left** starts; from
  top-right / bottom-left they are unsatisfiable or do not finish. They were
  planted/generated with a top-left-anchored solution (the mesh generator seeds
  the witness at cell (0,0)), so no layout exists from the far corners. Pinning
  `topleft_across` is therefore not an arbitrary choice - it is the start most
  of them work from at all. The start is **load-bearing** for the suite.
- **Corner matters, not the first word's direction.** `topleft_across` ≈
  `topleft_down` and `topright` ≈ `bottomleft` everywhere; there are two
  regimes, top-left vs the far corners.
- **On the one start-neutral fixture (`bundled_17`, real words), the start
  matters a lot and the pin is the *worst* choice for the production solver:**
  `mrv_inc` solves it in 22 k inferences from the far corners vs 322 k from
  top-left - a ~15x swing. So "corner start gives nicer/less-cramped layouts"
  is not supported here; the top-left choice is just *fixed*, not better.
- **"Best start" is fixture-dependent:** top-left for the combs (`benchmark_20`/
  `26` cost ~2.3-2.4x more from the far corners), far-corner for `bundled_17`,
  top-left-only for the planted fixtures.

Consequences: (1) benchmark results are all taken at the one start the synthetic
fixtures were built for, which is fine for *comparing strategies* but is not a
representative cross-section of "general crossword solving". (2) **The synthetic
suite cannot evaluate a "let the start fall out of the algorithm" idea** (e.g. a
principled deterministic seed choice) - it is rigged for top-left, so any seed
heuristic that picks elsewhere would make those fixtures look unsolvable. A fair
test needs start-neutral fixtures (real-word sets, or a generator that does not
anchor a corner).

---

## Result batch — 2026-06-23, F007 adopted (mrv_count aggregate_all, on top of F010) — current

Source: `benchmarks/results/2026-06-23-f007-strategy-matrix.csv`. After adopting
F007: `mrv_count/8` counts viable placements with `aggregate_all(count, ...)`
instead of `findall(t, ...)+length`. Stacks on the F010 batch below; same
30-iter / Warmup=3 methodology, so directly comparable. Median inferences.

| fixture (grid) | baseline | mrv (full) | mrv_capped | **mrv_inc** |
| --- | --- | --- | --- | --- |
| benchmark_08 (13) | 9.5 k | 83 k | 50 k | 52 k |
| benchmark_14 (17) | 35 k | 551 k | 265 k | 271 k |
| **benchmark_16_dense (17)** | **431.4 M** | 807 k | 370 k | **378 k** |
| benchmark_20 (37) | 43 k | 319 k | 319 k | **79 k** |
| benchmark_26 (49) | 76 k | 694 k | 695 k | **138 k** |
| bundled_17 (17) | 17 k | 343 k | 315 k | 312 k |
| benchmark_70_mesh (21) | 283 k | 12.62 M | 11.51 M | 5.02 M |

Headlines:
- Marginal and exactly as scoped: vs the F010 batch, the **21 mrv/mrv_capped/
  mrv_inc cells fall −0.09% to −0.37%, 0 worse**; the **7 baseline cells are
  byte-identical** (baseline never calls `mrv_count`), confirming the change is
  perfectly localized.
- The value is readability (`aggregate_all(count)` over `findall(t)+length`) at a
  small inference win — and it **refutes the audit's "~+60% worse on the mrv path"
  prior** (see the audit micro-opt section below).
- Behaviour-preserving: golden byte-identical; 81 plunit tests green.

## Result batch — 2026-06-23, F010 adopted (position/3 → nth1/3) — historical (superseded by the F007 batch above)

Source: `benchmarks/results/2026-06-23-f010-strategy-matrix.csv`. After adopting
F010: `find_intersecting_word/6` now uses the C-builtin `nth1/3` instead of the
hand-rolled `position/3` + `x_position/4` (both removed, including the dead
`:- false` clause). Same 30-iter / Warmup=3 methodology as the post-audit
baseline, so directly comparable. Median inferences (portable metric).

| fixture (grid) | baseline | mrv (full) | mrv_capped | **mrv_inc** |
| --- | --- | --- | --- | --- |
| benchmark_08 (13) | 9.5 k | 83 k | 51 k | 52 k |
| benchmark_14 (17) | 35 k | 552 k | 266 k | 271 k |
| **benchmark_16_dense (17)** | **431.4 M** | 807 k | 371 k | **378 k** |
| benchmark_20 (37) | 43 k | 320 k | 320 k | **79 k** |
| benchmark_26 (49) | 76 k | 696 k | 697 k | **139 k** |
| bundled_17 (17) | 17 k | 344 k | 316 k | 313 k |
| benchmark_70_mesh (21) | 283 k | 12.64 M | 11.53 M | 5.03 M |

Headlines:
- vs the 2026-06-22 post-audit baseline: **26/28 cells fewer inferences, 0 worse**,
  2 unchanged (baseline `benchmark_20`/`26`). Typical −2% to −3%; the headline is
  the one genuinely hard cell, baseline `benchmark_16_dense` **453.09 M → 431.37 M
  (−4.79%, ~21.7 M fewer)**. Strategy story otherwise unchanged.
- `nth1/3` is the C-builtin that the hand-rolled `position/3` reinvented; the win
  is per-lookup cost on `find_intersecting_word/6`, used by every strategy.
- Behaviour-preserving: golden output byte-identical; suite green at **81 plunit
  tests** — the 3 `position/3` unit tests were removed with the predicate. The full
  F006/F007/F010 measurement that led here is in the audit micro-opt section below.

## Result batch — 2026-06-22, post-audit re-verification (SWI 10.0.2) — historical (superseded by the 2026-06-23 F010 batch)

Source: `benchmarks/results/2026-06-22-post-audit-strategy-matrix.csv` (HEAD
`aaf01f4`). Full matrix re-run after the SWI-Prolog code audit
(`docs/prolog-audit-findings.md`) fixes landed: F016 (quality-engine
infinite-loop fix), F026 (benchmark harness throws on a fixture missing
`clues/1`), and F005 (`all_crossword/5` counts via `aggregate_all/3`).

**Result: inferences are byte-identical to the post-I5 batch below across all 28
cells (4 strategies × 7 fixtures), all solved.** No new table — the numbers in
the post-I5 batch carry forward unchanged. Wall ms differs (machine-dependent,
reporting-only, as always).

This is the expected outcome: none of the three fixes touch the `find_crossword`
search the matrix measures. F016 is in the greedy quality engine (`quality.pl`,
reached only via `--quality`); F026 changes the harness's fixture read only on
the *malformed* path (every committed fixture defines `clues/1`, so the measured
path is unchanged); F005 changes the `--all` solution-count path, not the solve.
Recorded as a dated regression checkpoint confirming the audit work is
perf-neutral on the solver. (F005's own large speedup is on the `--all` count
path, which the strategy matrix does not exercise — measured separately at
old/new = N+1 inferences, e.g. 9× at N=8.)

## Audit micro-optimization benchmarks — 2026-06-23 (SWI 10.0.2)

The three "NEEDS-BENCHMARK" items from the SWI audit
(`docs/prolog-audit-findings.md`): F006, F007, F010. Each was a clarity-vs-cost
trade the audit deliberately did **not** auto-recommend, leaving the call to
measurement. Measured here on the real suite; **F010 and F007 were subsequently
adopted** (see the 2026-06-23 result batches above), while F006 was left as-is.

Method: `benchmarks/run_matrix.pl`, all 4 strategies × 7 fixtures. Baseline =
the committed `2026-06-22-post-audit-strategy-matrix.csv`; re-run reproduced it
**byte-identical across all 28 cells** (determinism check on this machine).
Each variant was applied to a fresh `git checkout` of `crossword.pl` in an
isolated worktree and measured via a 1-iteration manifest (inferences are
iteration-independent, so the counts match the 30-iter manifest). Inferences are
the metric, as always.

- **F010 — `position/3`/`x_position/4` → `nth1/3` (+ delete the dead `:- false`
  clause).** Hits `find_intersecting_word/6`, used by *every* strategy. **WIN,
  larger than the audit guessed ("marginal").** 26/28 cells fewer inferences, 0
  worse, 2 unchanged (baseline `benchmark_20`/`26`). Typical −2% to −3%; headline
  is the one genuinely hard cell, baseline `benchmark_16_dense` **453.09 M →
  431.37 M (−4.79%, ~21.7 M fewer)**. Behaviour-preserving (golden byte-identical).
  **Adopted** — see the 2026-06-23 result batch above; on adoption `position/3`/
  `x_position/4` and their 3 unit tests are removed, leaving the suite at 81 tests,
  green. (The pre-adoption measurement was run with F010+F007 applied together —
  84 tests then, since `position/3` still existed.)
- **F007 — `mrv_count/8` `findall(t,…)+length` → `aggregate_all(count,…)`.**
  **The audit's prior ("likely-worse, ~+60% on the mrv/unbounded path") is
  refuted.** On the real suite it is marginally *faster*: all 21 mrv/mrv_capped/
  mrv_inc cells −0.07% to −0.36%, 0 worse; the 7 baseline cells are unchanged
  (baseline never calls `mrv_count`), confirming the change is correctly
  localized. A readability win at worst perf-neutral. (Whatever produced the
  +60% microbench figure did not reflect the predicate in situ — exactly why the
  audit punted to measurement.) **Adopted** — see the 2026-06-23 F007 result batch
  above (mrv cells −0.09% to −0.37% on adoption, baseline byte-identical, 81 tests
  + golden green).
- **F006 — strip-spaces idiom (`delete(Cs, ' ', …)` at 4 sites).** Two rewrites
  measured; **both lose, audit confirmed — keep the four inline `delete/3`s.**
  (a) `exclude(==(' '), …)` (meta-call): **all 28 cells worse**, +0.09% to
  **+4.36%**, total +0.35%; worst on `mrv_inc`, which hits the most strip sites.
  (b) A delete-based DRY helper (`strip_letters/3`): all 28 cells worse but
  tiny, +0.02% to +0.77%, total +0.02% — the per-call overhead buys only
  marginal de-duplication on a hot path. The perf fence holds.

Net: F010 (the material win) and F007 (marginal, mrv-only) are both adopted in
code (2026-06-23 batches above); F006 stays as-is (both rewrites lose).

## Result batch — 2026-06-20, post I5 fix (SWI 10.0.2) — current numbers (re-verified 2026-06-22)

Source: `benchmarks/results/2026-06-20-postI5-strategy-matrix.csv`. After the I5
fix (`no_word_merge/3` maximality check). Median inferences (portable metric);
wall ms in parentheses where notable. Each fixture on its manifest grid.

| fixture (grid) | baseline | mrv (full) | mrv_capped | **mrv_inc** |
| --- | --- | --- | --- | --- |
| benchmark_08 (13) | 9.7 k | 84 k | 51 k | 53 k |
| benchmark_14 (17) | 36 k | 567 k | 273 k | 278 k |
| **benchmark_16_dense (17)** | **453.1 M** (24 s) | 830 k | 380 k | **388 k (20 ms)** |
| benchmark_20 (37) | 43 k | 320 k | 321 k | **79 k (4.2 ms)** |
| benchmark_26 (49) | 76 k | 697 k | 698 k | **139 k (8.2 ms)** |
| bundled_17 (17) | 17 k | 353 k | 324 k | 321 k |
| benchmark_70_mesh (21) | **293 k (14 ms)** | 12.97 M | 11.93 M | 5.16 M (238 ms) |

Headlines:
- Strategy story unchanged: `mrv_inc` keeps the `dense_16` pruning win (~1,150x
  vs baseline) and the large-comb savings over `mrv_capped` (bench_20/26 ~4x
  fewer inferences). Inference counts rose ~5% everywhere vs the pre-I5 batch -
  the per-placement maximality check.
- **`benchmark_70_mesh` flipped:** pre-I5 baseline timed out; post-I5 it solves
  in 14 ms and now *beats* the MRV strategies. Its hardness was largely the
  collinear-overlap bug (baseline thrashing through invalid overlap branches).
  See I5 / I4 (Done). `benchmark_16_dense` remains the only genuine baseline-hard
  fixture - its equal-length full-width words cannot trigger the overlap bug.

## Result batch — 2026-06-20, pre I5 fix (SWI 10.0.2) — historical

Source: `benchmarks/results/2026-06-20-strategy-matrix.csv`. Before the I5 fix;
notable only because `benchmark_70_mesh` shows `baseline = timeout` here (the bug
made it baseline-hard). mrv_inc 161 ms / 3.18 M, mrv 367 ms, mrv_capped 479 ms.
Superseded by the post-I5 batch above.

## Result batch — 2026-06-19 (SWI 10.0.2) — historical (pre-I4 fixture)

Source: `benchmarks/results/2026-06-19-strategy-matrix.csv`. Same six fixtures
as 2026-06-20 minus `benchmark_64_mesh`; numbers match. Superseded by the
2026-06-20 batch.

## Result batch — 2026-06-18 (SWI 10.0.2) — historical (pre-mrv_inc)

Source: `benchmarks/results/2026-06-18-strategy-matrix.csv`. Median wall (ms)
and inferences, each fixture on its manifest grid.

| fixture (grid) | baseline | mrv (full) | mrv_capped |
| --- | --- | --- | --- |
| benchmark_08 (13) | 0.5 ms / 9.0 k | 3.9 ms / 76 k | 2.5 ms / 47 k |
| benchmark_14 (17) | 2.0 ms / 33 k | 25.5 ms / 501 k | 12.9 ms / 247 k |
| **benchmark_16_dense (17)** | **26,256 ms / 450.9 M** | 36.8 ms / 718 k | **18.2 ms / 340 k** |
| benchmark_20 (37) | 2.3 ms / 37 k | 15.6 ms / 272 k | 15.9 ms / 272 k |
| benchmark_26 (49) | 3.8 ms / 65 k | 36.8 ms / 585 k | 37.2 ms / 586 k |
| bundled_17 (17) | 1.0 ms / 16 k | 16.7 ms / 328 k | 15.9 ms / 308 k |

Headline: on the one pathological (dense mesh) fixture, `mrv_capped` is a
**~1,330x reduction in inferences** (~1,450x wall). On every other fixture
(easy in input order) the MRV variants are a net tax — all still < 40 ms. The
comb tax motivated E5 (`mrv_inc`).

---

## Entries

### E1 — baseline (input-order search) — REFERENCE, kept

- **Strategy id:** `baseline` (`crossword.pl`).
- **What:** original algorithm. `assign_words` picks the next word with
  `member/2` in input order; `find_intersecting_word` rejects words that can't
  connect. No search guidance.
- **Result:** fast on every fixture *except* the dense mesh, where it explodes
  to ~26 s / 450 M inferences.
- **Verdict:** keep permanently as the comparison reference. Not the algorithm
  to ship — one pathological input is enough to motivate ordering.

### E2 — MRV, full exact counts — superseded by E3

- **Strategy id:** `mrv` (`crossword.pl`); `mrv_cap(mrv, unbounded)`.
- **Hypothesis:** fail-first (minimum-remaining-values) ordering tames the
  dense blow-up.
- **What:** at each step, count *every* viable placement of each remaining
  word, try words most-constrained-first, **backtrackably** (reorders the same
  tree; completeness preserved).
- **Correctness subtlety (important):** placements must cross an
  already-placed word, so a word's option count *grows* as neighbours land — a
  count of 0 means "not connectable yet", **not** "dead". An earlier prototype
  treated 0 as a dead branch and as grounds to commit deterministically; that
  was *incomplete* and failed fixtures the baseline solved. Fix: only rank
  words with count > 0, keep selection backtrackable, and special-case the
  first (seed) word as a free choice over all words.
- **Result:** dense 450.9 M → 718 k inferences (huge win); but a per-node
  recount tax makes easy fixtures 5–20x slower.
- **Verdict:** proves the win exists, but **strictly dominated by E3** (capped
  is ≤ in every cell). Kept as a documented comparison point only.

### E3 — MRV, count capped at 2 — KEPT (recommended)

- **Strategy id:** `mrv_capped` (`crossword.pl`); `mrv_cap(mrv_capped, 2)`.
- **Hypothesis:** the ordering only needs the buckets 0 / 1 / ≥2, so capping
  the count enumeration at 2 keeps the pruning while cutting the per-node cost.
- **Result vs E2 (inferences):** mesh fixtures roughly halve (dense 718 k →
  340 k; bench_14 501 k → 247 k; bench_08 76 k → 47 k). On the comb fixtures
  (20, 26) the cap is **inert** (272 k/586 k, identical to full MRV) — comb
  words have ≤ 2 placements each, so the cap never triggers.
- **Result vs baseline:** dense ~1,330x fewer inferences; elsewhere a tax of
  ~5–19x that is largest on the big comb grids (per-node recount over many
  words on a large grid, unamortized) but always < 40 ms absolute.
- **Verdict:** strong, but **superseded by E5** (`mrv_inc`), which keeps E3's
  pruning and removes most of its large-grid per-node tax. Kept as a strategy
  for comparison.
- **Follow-up (resolved by E5):** the residual tax is fixed per-node *recompute*
  cost (recounting every word every node), not enumeration depth.

### E4 — static "longest word first" — REJECTED

- **Idea (no surviving code):** approximate MRV for free by sorting words by
  length descending once (longest words have fewer slots), then run the
  baseline `member/2` search on the pre-sorted list. Sketch:
  `order_by_length/2` via `map_list_to_pairs(atom_length, …)`, `keysort`,
  `reverse`, used in `find_crossword` before `assign_words`.
- **Result:** catastrophic. **Timed out (>120 s)** on `benchmark_14`, which
  the baseline solves in ~2 ms, and on the dense fixture. Even `benchmark_08`
  regressed ~7x.
- **Why:** a static order is blind to grid state, which is the entire value of
  MRV. Reordering globally can push the search into a far worse region; "fewer
  slots in the abstract" ≠ "most constrained right now".
- **Verdict:** rejected — do not revisit static/length-based ordering. If
  ordering is wanted, it must be **dynamic and state-aware** (E3).

### E5 — MRV capped + incremental count cache — KEPT (recommended) — was I1

- **Strategy id:** `mrv_inc` (`crossword.pl`); own driver `assign_words_inc/9`.
- **Hypothesis:** E3's residual tax is recomputing every remaining word's count
  at every node; recompute only the counts that can have changed.
- **What:** thread a count cache. After placing word W, only recount words that
  **share a letter** with W; carry the rest forward. Counts capped at 2 (as E3).
- **Correctness invariant:** a placement must cross an already-placed word, so
  placing W can only ADD options to words sharing a letter with W, and only
  REMOVE options from the rest. So a carried-forward count is always **>= the
  true count** — never an under-count, which is what would be unsafe (it could
  prune a placeable word and break completeness). Over-counts are harmless: the
  word is tried, `find_intersecting_word` fails it, search moves on.
- **Bug found & fixed during impl:** first cut built the candidate list with
  `findall(Count-W, …)`, which **copies** terms; the copied `Entry` no longer
  matched the original word in `Words`, so the `==`-based `remove_x` never
  removed it and the search looped forever (even on 2 words). Fix: use
  `map_list_to_pairs`, which keeps the original terms.
- **Result (2026-06-19 batch):** large combs improve sharply — bench_20 272 k →
  72 k inferences (3.8x; 15.7 → 4.1 ms), bench_26 586 k → 126 k (4.6x; 36.7 →
  8.0 ms); comb tax over baseline drops from ~7-9x to ~1.9x. Dense win
  preserved (348 k vs 451 M ≈ 1,300x). Small/mesh fixtures within ~2-3% of E3
  (cache bookkeeping overhead, immaterial).
- **Verdict:** **SHIPPED as the production default** (`default_strategy(mrv_inc)`
  in `crossword.pl`; CLI `--strategy` default and the `crossword/3` +
  `find_crossword/5` wrappers all resolve to it; golden regenerated). Best
  general strategy: E3's pruning with the large-grid tax mostly removed, no
  regressions, valid output, same solvability verdicts as baseline. Residual
  small-fixture tax over baseline (few nodes, so the cache can't help) is
  unchanged but < 16 ms absolute.

---

## Open ideas

- All logged algorithm ideas (I1-I6) are closed. Possible future work: a
  *cheaper* value heuristic that biases `find_intersecting_word/6`'s generation
  order without the full enumerate-and-sort (I2 failed on that cost, not the
  idea); and the infra items deferred by choice (auto-stamp git SHA + a
  fixture-set hash in `run_matrix.pl`; a CI / inference-count regression gate).

## Done / closed ideas

- **I1 — incremental counts.** DONE → E5 (`mrv_inc`). Shipped; 3.8-4.6x fewer
  inferences on the large combs, dense win preserved.
- **I2 — value ordering.** TRIED, REJECTED (reverted from code, like E4). Idea:
  among the chosen word's candidate placements, try them ordered by how many
  crossings each makes with already-placed words, instead of
  `find_intersecting_word/6`'s native order. Sketch: in the placement step,
  `findall(Count-(S-D), distinct(S-D, (find_intersecting_word(...,S,D),`
  `crossing_count(Letters,S,D,GridLen,GIn,Count))), Pairs)`, `sort/4` by Count,
  then `member` best-first; `crossing_count/6` counts word cells already holding
  the matching letter.
  - **Result (post-I5 solver):** both directions are *worse* than `mrv_inc`.
    Most-crossings-first: `benchmark_16_dense` 398 k -> 488 k inferences,
    `benchmark_14` 289 k -> 358 k (~20-24% more); `benchmark_20` ~3% more.
    Fewest-crossings-first (the least-constraining-value intuition): essentially
    identical (dense 491 k, bench_14 360 k).
  - **Why:** the ~20% penalty is the per-node cost of enumerating and sorting
    *all* of a word's placements - inherent to value ordering. Both heuristic
    directions lose by the same margin, so the loss is the overhead, not the
    ranking. On top of `mrv_inc`'s strong variable ordering the search already
    backtracks little, so there is not ~20% of pruning available to offset the
    cost. (Value ordering reorders placements *within* a word; it cannot help
    baseline's `dense_16` blow-up, which is a *variable*-ordering problem.)
  - **Verdict:** not worth it for this suite. Revisit only with a value
    heuristic that does not require full per-node placement enumeration, and/or
    a fixture where `mrv_inc` itself backtracks heavily (none found).
- **I6 — `mrv_inc_deg` (MRV + degree tie-break).** TRIED, REJECTED (reverted
  from code, like E4/I2). Idea: among the words tied at the same capped placement
  count, place the highest-DEGREE one first - degree = the number of other
  remaining words sharing a letter (the count of future crossing constraints) -
  the classic MRV + degree fail-first pairing. Implemented as a `mrv_inc` variant
  reusing the incremental count cache, ordering by key `Count - (-Degree)` so
  `keysort` gives count ascending then degree descending; `mrv_inc` itself left
  byte-identical (golden + 66 tests unchanged).
  - **Result:** worse than `mrv_inc` on EVERY fixture (median inferences):
    bundled_17 321k->511k, bench_08 53k->58k, bench_14 278k->317k,
    `benchmark_16_dense` 388k->447k (+15%), bench_20 79k->252k (+220%),
    bench_26 139k->583k (+320%), `benchmark_70_mesh` 5.16M->14.57M (+182%).
    Full batch: `benchmarks/results/2026-06-21-mrv_inc_deg.csv`.
  - **Why (node counts settle it):** the degree ordering gives NO search benefit.
    Instrumented node counts are IDENTICAL where it matters - `dense_16` 16 = 16,
    bench_26 26 = 26 - so the entire inference increase there is pure per-node
    cost: an O(n^2) degree scan (each ranked word vs every remaining word) at
    every node, which balloons with the word-set size (hence +220-320% on the
    larger sets, whose search is near-linear with no pruning to gain). Worse, on
    `benchmark_70_mesh` the degree ordering actively STEERED INTO A WORSE TREE:
    `mrv_inc` solves in a clean 70-node forward pass; `mrv_inc_deg` blows up to
    775 nodes (11x). So degree loses three ways: no node reduction anywhere (a
    tie at best), an 11x tree blow-up on one fixture, plus O(n^2) per-node cost.
  - **Verdict:** not worth it, and - unlike I2 - NOT worth revisiting with a
    cheaper degree computation: the ordering produced zero tree reduction even
    when notionally free, so cost is not the only problem. The capped-count
    buckets (0 / 1 / >=2) already place forced words first; among the large `>=2`
    bucket, "most-connected first" is not a useful discriminator for this suite
    (and on one fixture actively misleads). `mrv_inc` stays the default.
- **I5 — solver soundness: collinear "word inside a word".** DONE (fixed in
  `crossword.pl`). The solver could place a word and a collinear word whose span
  contains it (e.g. `DFEE` inside `DFEEE`, same start cell, same direction) -
  `check_prev_cell`/`check_next_cell` validate a word's boundaries only at its
  own placement time, and the `X==L` branch of `assign_letters/7` lets a later
  longer word pass straight through and fill a shorter word's end boundary. The
  result is two same-direction answers at one start cell, which
  `assign_clue_numbers/2` cannot number, so the CLI backtracked into the search
  and hung (the benchmark harness uses `find_crossword` with a cut, so it never
  hit numbering - "solves" != "full pipeline succeeds").
  - **Fix:** `no_word_merge/3` in `assign_word/10` - a newly placed word's cells
    must not occupy the boundary cell (just before start / after end, in its
    direction) of any already-placed word. With the existing checks (which guard
    the new word's own boundaries) this maintains "every word is a maximal run"
    inductively. Forbids the prefix, suffix and interior-containment cases.
    Regression tests added (`rejects_collinear_prefix_overlap`,
    `allows_separated_collinear_word`). Validated end-to-end: a prefix-laden mesh
    that hung the old CLI now emits valid JSON. Golden unchanged (no existing
    fixture's first solution contained an overlap).
  - **Surprise / important:** the fix made `benchmark_70_mesh` (I4) **baseline-
    easy** - baseline went from timeout (>120 s) to 14 ms. Its baseline-hardness
    was largely the bug: baseline was thrashing through invalid overlap branches.
    Probing under the corrected solver found **no robust baseline-hard /
    MRV-tractable mesh window** (grid 21: N<=70 baseline-easy, N>=72 hard for
    all). So `benchmark_16_dense` remains the suite's only genuine baseline-hard
    fixture - and its hardness is real (equal-length full-width words cannot form
    a collinear containment, so the bug never applied to it). The generator's
    prefix-pair rejection is now redundant for correctness (solver enforces it);
    kept only to keep witnesses tidy and the fixture reproducible.
- **I4 — hard large-grid fixture.** DONE → shipped `fixtures/benchmark_70_mesh_words.pl`
  (`benchmarks/gen_mesh_fixture.py 21 70 8 3 5 1`, deterministic; grid 21, 70
  short words, alphabet 8, seed 1). Generator builds a witness layout by
  simulating the solver's own placement rules (short perpendicular words
  crossing placed ones, seed at (0,0) across), so it is **satisfiable +
  reachable by construction**; output order is shuffled so baseline cannot
  replay the build.
  - **Result (as shipped, PRE-I5):** baseline timed out (>120 s), while `mrv`
    (367 ms), `mrv_capped` (479 ms) and `mrv_inc` (161 ms) all solved — looked
    like the first suite fixture where baseline is hard but MRV wins on a 21x21
    grid. Added a 60 s per-cell timeout to `run_matrix.pl` so such fixtures
    cannot hang the matrix.
  - **SUPERSEDED by I5:** that baseline-hardness was largely the collinear-overlap
    bug. Once I5 fixed it, baseline solves this fixture in 14 ms and now beats
    the MRV strategies; no robust baseline-hard mesh window survives (see I5).
    The fixture is kept, relabelled as a large dense mesh; `benchmark_16_dense`
    is the genuine baseline-hard fixture.
  - **Key finding (the bigger lesson):** MRV's advantage is **structure-specific,
    not generic.** On *sparse/flexible* random meshes baseline wins (its cheaper
    per-node cost beats MRV, which has nothing to prune); only *near saturation*
    does a search hard enough to favour pruning appear, and there the cliff is
    **instance-dependent** — density sweeps at grid 21 found neighbouring
    instances favouring opposite strategies (e.g. ~62 words baseline-easy/
    mrv_inc-timeout; ~70 words baseline-timeout/mrv_inc-fast; >=72 both timeout).
    `benchmark_70_mesh` is the chosen baseline-hard instance.
  - **Honest caveat:** the fixture sits near that cliff, so it is a point-in-the-
    design-space discriminator, not a smooth "MRV wins by Nx"; the committed file
    is deterministic so the *numbers* are reproducible.
  - **Latent solver bug surfaced (→ I5):** an earlier candidate exposed that the
    solver can place a word and its **prefix** collinearly at the same start cell
    (e.g. `DFEE` and `DFEEE` both across from one cell) — `check_next_cell/4`
    validates at placement time but a later, longer word retroactively fills a
    shorter word's "next" cell. `assign_clue_numbers/2` then fails (a start cell
    with two same-direction words), so the CLI backtracks into the hard search
    and hangs. Worked around in the generator (reject prefix pairs); the real fix
    belongs in the solver — see I5.
- **I3 — hard large-grid fixture.** ATTEMPTED, no fixture shipped. A lattice
  generator (full-width words on even rows/cols of a deterministic letter grid;
  satisfiable + adjacency-valid by construction, alphabet/stride as difficulty
  knobs) was built and swept. **Finding:** full-width lattice words are
  *pathologically* hard regardless of parameters — even grid 17 with a near-
  unique alphabet times out under `mrv_capped` (>30 s), because a full-width
  word can spuriously "cross" a perpendicular word at any matching letter and
  only fail deep. Difficulty is not tunable this way; the generator was removed.
  The pruning-vs-overhead discrimination I3 wanted is instead obtained **across
  the existing suite**: combs (large grid, easy search) isolate per-node
  overhead, `dense_16` (small grid, hard search) isolates pruning — which is
  what let E5 be evaluated. A better fixture is tracked as I4.
  *(Op note: those sweeps balloon Prolog stacks into the GBs before any timeout
  fires — run heavy solves under `ulimit -v` to avoid an OOM kill.)*

---

## Crosswordsmith arrange campaign — 2026-07 (post-pivot engine)

Everything above this section predates the pivot to the crosswordsmith
arrange engine (the rebuilt strict construct/rescore/emit path in
`prolog/crosswordsmith/arrange.pl` over `core.pl`'s mrv_inc). The kept
learnings (mrv_inc, capped counts, the I5 merge guard) are folded into the
current engine; the REJECTED verdicts above (I2 value ordering, I6 degree,
E4 static ordering) were adjudicated on the old engine + old fixture suite
and are treated as advisory, not binding. The quantitative reference for this
campaign is `benchmarks/baseline.json` over the 9/15/21 cost ladder
(`benchmarks/workloads.pl`); never compare against pre-pivot numbers.
Background research distilled for this campaign lives in `docs/research/`.

### E-H1 — strict-mode transpose-pair corner dedup — KEPT (shipped)

- **Where:** `arrange.pl` `arrange_corners/1` (strict construct path only);
  commit 30628c7.
- **Hypothesis:** the 4 start corners in `start_locs/1`
  ([topleft_across, topleft_down, topright, bottomleft]) form two
  TRANSPOSE-PAIRS — {topleft_across, topleft_down} and {topright, bottomleft}.
  Transposing the grid (row<->col) maps across<->down while preserving letter
  order, so a pair's mrv_inc search trees are exact transposes: same first
  solution modulo transposition, and `layout_reward` is transpose-invariant
  (it counts checked cells), so a pair's first-solution rewards are identical.
  Strict mode therefore did ~2x redundant work sweeping all four.
- **What:** scope strict `arrange_best_layout/6` to a new
  `arrange_corners([topleft_across, topright])` — one representative per pair,
  the same two the reward-tie break already prefers (stable `sort(1,@>=)`), so
  the emitted best is byte-unchanged. `start_locs/1` stays 4-corner for the
  greedy/best-effort/candidates paths (transposed layouts are legitimately
  distinct candidates there) and the fragment/enumerate paths.
- **Verification (before changing anything):** per-corner first solution over
  every ladder rung (grids 9/15/21, 8..80 words). Within each pair, rewards are
  identical AND the partner placement is the literal grid-transpose of its
  mate, on all 9 rungs incl. the two heaviest. The two pairs are NOT equal to
  each other (e.g. 15x15_32w: pair-A 347 vs pair-B 349), so both
  representatives are retained.
- **Result (search inferences, `make bench-check --heavy` vs the pre-change
  baseline):** -49.6% to -50.1% on every rung (09x09_08w 98,058 -> 49,035;
  15x15_36w 182.6M -> 91.8M). No regressions. `make test`: 186 plunit pass,
  goldens byte-identical. Baseline re-recorded at the new counts
  (`make bench-record BENCH_ARGS=--heavy`).
- **Verdict:** SHIPPED. A ~2x strict-search win at zero output change, since
  the dropped corners were exact transpose-twins of the kept ones. Guard:
  relies on a square canvas + transpose-symmetric branch step/reward — do not
  route the candidate/fragment/enumerate paths through `arrange_corners/1`.

### E-H2 — var-cell grid term replaces the assoc grid — KEPT (shipped)

- **Where:** `core.pl` (init_grid, assign_letters, check_prev/next_cell,
  adj_is_free), `arrange.pl` (crossing_count, clashing_cell,
  word_best_placement); merge 56bad96 (branch commit 79bca05).
- **Hypothesis:** the AVL-assoc grid (cell -> `empty`|letter) was the top
  measured constant-factor cost: every cell read O(log N^2) via
  `$btree_find_node`, every write an allocated rebalanced path, all inside
  mrv_count's innermost loop (~15% self-time on a mid rung). A single
  compound `grid(C1..C(N^2))` — unbound var = empty, bound atom = letter —
  makes reads arg/3 + var/nonvar (O(1), allocation-free) and writes a
  unification undone by the trail, with an identical search tree.
- **What:** pure representation swap; GIn/GOut collapse to one term. Two
  forced adjustments: (1) assign_word rejects Start < 1 up front (arg/3
  THROWS on negative indices where get_assoc merely failed — the old code
  relied on that failure); (2) greedy word_best_placement scores BEFORE
  assign_word, since the in-place bind would otherwise count the word's own
  letters as crossings. White-box init_grid tests rewritten to the new rep.
- **Result (vs the post-E-H1 baseline, `make bench-check --heavy`):** -14.8%
  to -23.9% search inferences on every rung (15x15_36w 91.75M -> 72.79M);
  wall medians down ~13-18% on the heavy tail; RSS flat; goldens
  byte-identical; 186 plunit pass. Cumulative with E-H1: every rung -55% to
  -62% vs the campaign start (36w rung 182.6M -> 72.8M).
- **Verdict:** SHIPPED. Same search order, same output, materially fewer
  inferences and less GC pressure (doubly valuable under WASM, whose
  setjmp-heavy backtracking path is the expensive one — see
  docs/research/swi-vm-wasm-performance.md).
- **Caveat for later:** the greedy-path scoring reorder means placement_key
  now runs for every crossing candidate, not only legal ones (the old
  ordering's micro-optimization). The greedy path is unbenchmarked; if
  best-effort/candidates rungs are ever added, revisit that ordering.

### E-H5 — saturating placement counter for the capped mrv_count path — KEPT (small)

- **Where:** `mrv_count/8` (`core.pl`) Cap=2 path (mrv_capped + production
  mrv_inc); the mrv/unbounded path stays on aggregate_all. Merge d59f0d4
  (branch commit b766f4d).
- **What:** `count_upto2/2` — the SWI-manual nb_setarg failure-driven counter
  with early exit at the 2nd solution, wrapped in `\+` so Start/Dir bindings
  are discarded like findall; fresh counter(0) holder per call (reentrant).
  Candidates measured on real mid-search goal tuples: aggregate_all+limit
  (ref); nb_setarg loop -0.15% (winner); staged call_nth +24%; \+ \+ probes
  +49% — count=1 is the PLURALITY case at search nodes, so re-running
  variants lose badly (useful fact for future counting code).
- **Result:** -0.04% to -0.24% per rung, 0 regressions, goldens
  byte-identical; new count_upto2 plunit unit locks 0/1/2 + no-residual-
  binding semantics. R1 research (docs/research/swi-vm-wasm-performance.md)
  had already corrected the original hypothesis: the profiled has_type tax
  belonged to list_to_set/2, not aggregate_all — hence the small size.
- **Verdict:** KEPT (marginal but monotonic and behaviour-locked). Ratcheted
  together with E-H3.

### E-H3 — tabled crossing/letter memos on the mrv_inc hot path — KEPT (shipped)

- **Where:** `core.pl` — `:- table pair_crossings/3` (the ordered x(PPos,Pos)
  crossing sequence per ordered letter-list pair, replacing per-call
  intersection/3 + list_to_set/2 + two nth1/3 backtracking passes in
  find_intersecting_word/6) and `:- table answer_letters/2` (the atom_chars +
  delete/3 letter footprint behind entry_letters/2). Merge 365c7a6 (branch
  commit 43c16b3).
- **Order preservation (hard requirement, held):** the tabled findall lifts
  the former inline conjunction verbatim, so candidate order and every golden
  are byte-identical.
- **Determinism invariant (load-bearing):** tables are abolished at the
  mrv_inc find_crossword/6 boundary. A search amortizes ~|Words|^2 pair
  entries across millions of nodes but starts cold, so search_inf is
  self-contained and run-order independent (core-only == --heavy counts,
  verified). Do NOT remove the abolish for a warm-table micro-win.
- **Result (standalone, vs pre-campaign baseline):** -8% to -30.6% per rung,
  heavy tail gains most (36w 182.6M -> 126.7M standalone). Composed with
  E-H1/E-H2/E-H5 and ratcheted: 36w now 44.8M (-75.5% cumulative).
- **Cost:** whole-process RSS rises with table size (21x21_80w 14.2 ->
  24.9 MiB, +75%; mid rungs +13-20%), bounded per search, freed at the
  boundary. Tabling is available under the WASM build (single-threaded, so
  shared-table caveats moot) — but the footprint matters for the ~300MB
  Android ceiling; revisit if WASM memory becomes tight.
- **Verdict:** SHIPPED. The biggest single constant-factor win of the
  campaign after the corner dedup.

### P1 — backtrack-distance / failure-clustering probe — MEASUREMENT (no merge)

- **What:** nb_setval place/unplace/wipeout instrumentation on
  assign_words_inc/select_inc, run on the dense rungs from both production
  corners. Artifacts on probe branch (worktree commit a74ffcf+0b32369):
  benchmarks/probe_backtrack.pl + scratch/probe_p1_output.txt.
- **Findings:**
  1. CBJ headroom: NONE. Retreats are 84-86% distance 1-2, the distribution
     dies by 5, max observed 7. Conflict-directed backjumping / dynamic
     backtracking is NOT worth building for this search. (Closes the
     question docs/research/arrange-search-algorithms.md #2 left open.)
  2. wdeg signal: STRONG. Wipeouts concentrate on ~8-10 words per rung
     (top-5 = 75-91%) — an adaptive weighted tie-break has real signal.
     Addressable surface: the thrashing corner on 34w/36w only.
  3. 21x21_80w places all 80 words with ZERO backtracking from both
     corners: its ~10M-inference cost is pure per-node counting work.
     Only constant-factor cuts help it; no ordering change can.
  4. Churn is corner-asymmetric (36w: topleft_across 4 unplaces vs topright
     3,418) and broad-but-shallow (thousands of 1-2-deep retreats) — the
     signature of a locally-suboptimal ordering, favouring wdeg over any
     backjumping scheme.

### E-H6 — failure-driven word weights (dom/wdeg tie-break) — TRIED, REJECTED (reverted from code)

- **Where (measured, then reverted):** `core.pl` — a non-backtrackable weight
  store (`:- dynamic word_weight/2` keyed by answer atom, `weights_active/0`
  guard), reset per search via `reset_search_weights/0` at both mrv_inc entries
  (`find_crossword(mrv_inc,...)` and `construct_from_seed/5`), bumped on WIPEOUT
  (select_inc filters Placeable to [] → blame the just-placed word, the classic
  dom/wdeg analogue), and read as a SECONDARY sort key inside `select_inc`
  (new `order_placeable/3` + `weight_key/2`: keysort on `k(Count,-Weight)` →
  count ascending, then weight descending, input order last). No surviving code.
- **Hypothesis (P1 #2):** wipeouts concentrate on ~8-10 words per dense rung
  (top-5 = 75-91%), so failure-driven weights used ONLY as an equal-count
  tie-break should steer the dense 34w/36w search away from its shallow thrash
  for a large inference cut, without regressing the backtrack-free rungs.
- **Design (held to the guidance):** weights SURVIVE backtracking (the point —
  a threaded arg would undo the bump exactly when it matters), so a
  non-backtrackable dynamic store; reset per corner (self-contained like the
  E-H3 table abolish, so counts stay run-order independent); tie-break within
  equal capped-count buckets ONLY (fail-first preserved); fully deterministic
  (no RNG on the default path); seeded path left byte-identical. Per-node cost
  is O(1)-ish: a first-arg-indexed `word_weight/2` read per candidate, and — key
  design win — a `weights_active/0` guard so a search with no recorded failure
  yet (and any backtrack-free rung) takes the ORIGINAL weight-free ordering,
  byte-identical in both output and cost.
- **The signal fired exactly as P1 predicted, and still bought nothing.**
  Per-corner (the thrashing `topright`): 34w 1,756 bumps on 7 words
  (FCC 583, EFD 345, ...), 36w 2,095 bumps on 12 words (DEAA 396, CBF 338, ...);
  `21x21_80w` recorded ZERO bumps from both corners (guard kept it inert). The
  blame concentrated on the right culprits — but reordering equal-count
  candidates by that blame did not shrink the tree: 34w topright ≈ 44.2M→44.0M,
  and the full ladder moved within noise.
- **Result (`make bench-check --heavy`, warm, vs baseline.json):**

  | rung | baseline | measured | delta | note |
  | --- | --- | --- | --- | --- |
  | 09x09_08w | 33,615 | (win) | ~-5% (cold probe) | light-rung fluke |
  | 15x15_12w | 126,842 | 117,329 | **-7.50%** | light-rung fluke |
  | 21x21_25w | 421,556 | 421,734 | +0.04% | inert |
  | 15x15_28w | 557,743 | 557,866 | +0.02% | inert |
  | 15x15_32w | 1,458,240 | 1,458,804 | +0.04% | inert |
  | **09x09_16w** | 868,655 | 906,771 | **+4.39%** | **REGRESSION** |
  | 21x21_80w | 10,005,936 | 10,006,258 | +0.00% | byte-identical (guard) |
  | **15x15_34w** | 18,722,171 | 18,589,349 | **-0.71%** | target rung, within noise |
  | **15x15_36w** | 44,804,532 | 44,605,072 | **-0.45%** | target rung, within noise |

  Ratchet verdict: **FAIL** (1 regression). Reward (arrange_best_layout) held
  on every rung EXCEPT **15x15_36w, which dropped 395 → 392** — a strictly WORSE
  layout for a sub-noise "win". All other rewards unchanged.
- **Why it fails (three ways, echoing I6):** (1) NO tree reduction on the target
  — on top of mrv_inc's strong MRV variable ordering the low-count buckets are
  mostly forced (0/1 placements), so there are few equal-count ties at the
  critical shallow-retreat nodes for a tie-break to act on; floating a blamed
  word within a tie just picks a DIFFERENT equally-doomed branch (the retreats
  are structural near saturation). (2) It STEERED `09x09_16w` INTO A WORSE TREE
  (+4.39%, reward unchanged so purely more work) — the same "reordering pushes
  the search into a worse region" failure I6 hit on the mesh, and E4 hit
  globally. (3) The only "wins" are light rungs (12w -7.5%, 08w) where a lucky
  reorder found a shorter first path — not robust and not the target. A broader
  bump event (event (c), bump on EVERY unplace via a backtracking-safe
  `(true;bump,fail)` hook) was also measured: it recorded MORE bumps (34w 3,072,
  36w 3,380 — matching P1's unplace counts) but gave the SAME sub-noise tree
  (34w topright 17.5M, 36w topright 44.0M) while adding a per-node choicepoint on
  the hot path, so it was discarded in favour of the cleaner wipeout event.
- **Validation:** `make test` — 201 plunit pass, and ALL goldens byte-identical
  (the shipped golden fixtures — bundled_17 fixed/max/fragment/candidates — never
  wipe out on their winning path, so their first solution is unchanged; the only
  output change was the 36w LADDER reward, not a golden). `bundled_17` sanity:
  size 15 stays infeasible; size 17 reward 60/placed 6 unchanged.
- **Verdict:** REJECTED — reverted from code (like I2/I6/E4). The P1 wdeg signal
  is real and lands on the right words, but an equal-count tie-break cannot
  convert that concentration into a search cut here: the addressable target
  (34w/36w) moves within noise while a lighter dense rung (09x09_16w) regresses
  +4.39% and 36w's layout reward drops. This closes the
  docs/research/arrange-search-algorithms.md #1 (wdeg) bet as a NEGATIVE result.
  Do not revisit as a within-bucket tie-break; any future attempt must change the
  PRIMARY ordering (which I6/E4 already show this search punishes) or attack the
  per-node counting cost that dominates the backtrack-free rungs instead.

### E-H7 — placed-word record: dict -> pw/8 compound — KEPT (shipped)

- **Where:** all placed-word consumers — core.pl (assign_word builds
  pw(Answer, Letters, Cells, Dir, Len, Start, End, Num); accessors
  pw_answer/2..pw_num/2 exported; placed_boundary_cell reads the
  precomputed End instead of last/2 per candidate; add_clue_word binds the
  fresh Num var via pure rebuild), metrics.pl, arrange.pl, lint.pl,
  stockgrid.pl, fill.pl (the metric layer is duck-typed and shared by the
  solver AND validator paths, so all moved in lockstep; lint/stockgrid/fill
  build records with unused fields unbound). export.pl reads only emitted
  JSON — untouched. Merge 43a4729 (branch commit c2e3ebd).
- **Hypothesis:** placed-word get_dict field reads (~1.5x compound access;
  hot in find_intersecting_word, no_word_merge, the scorers) plus the
  per-candidate last(Cells) recomputation are a residual constant factor.
- **Result (composed on Wave 1, vs the ratcheted baseline):** -1.9% to
  -9.7% per rung, 0 regressions; goldens byte-identical; 201 plunit pass.
  Composed value EXCEEDS the standalone-vs-old-tree measurement (-0.95%..
  -7.0%): Wave 1 shrank the surrounding costs, so the dict reads were a
  larger share of what remained — complementary, not overlapping. Gradient:
  broad forward searches win most (21x21_80w -9.7%, many placed-word field
  reads per candidate), the trail-dominated 36w rung least (-1.9%).
- **Verdict:** SHIPPED. And per the honest diminishing-returns read: the
  representation seam is now MINED OUT — Letters-as-packed-atom would lose
  (the crossing machinery needs lists), cells-as-compound buys nothing.
  The one remaining constant-factor idea with substance is algorithmic:
  a per-word letter-presence bitmask to short-circuit non-crossing pairs
  before the pair_crossings lookup and inside inc_count's shares_letter
  (tracked as E-H8, the campaign's final perf experiment).

### E-H8 — per-word letter-bitmask prefilter — REJECTED (no surviving code)

- **Idea:** an integer letter-presence mask per answer (bit = 1 << (char
  code /\ 63); exact on real crossword alphabets, conservative under
  collision) to short-circuit "do these words share a letter?" in O(1) AND:
  gating find_intersecting_word's per-placed-word crossing enumeration and
  replacing inc_count's shares_letter list scan.
- **Execution caveat:** the experiment ran against the PRE-Wave-1 tree
  (worktree branched at abbc770), so its gate target was the old inline
  intersection/3 rather than the current tabled pair_crossings lookup. The
  rejection stands anyway because the decisive measurement is
  base-independent (below). Rejected diff inspectable on the experiment
  branch (commit bd1c064); nothing merged.
- **Result:** +1.3% to +2.9% inferences on 8/9 ladder rungs (0 wins);
  bundled_17 (real words) +0.65%. Root cause: the prefilter almost never
  fires — measured AND=0 skip rates 2.6-13.7% (find_intersecting_word) and
  10-17% (inc_count). The brief's alphabet intuition INVERTED: real answers
  are letter-DENSE (esp. multi-word answers), so real-word masks are dense
  and pairs nearly always share a bit — bundled_17 had the WORST yield
  (2.6%). Per-word letter density, not alphabet size, governs disjointness.
- **Verdict:** REJECTED, do not revisit on this solver. With <=14% skippable
  pairs, no gate-target cost (inline or tabled) leaves room for the mask
  constant paid on the 86-97% of pairs that pass through. This closes the
  constant-factor idea queue: together with E-H7's "representation seam
  mined out" verdict, the campaign's per-node cost reduction is done.

### P2 — density-envelope close-out probe — MEASUREMENT (rungs adopted)

- **What:** density sweeps past the pre-campaign cliffs at all three grid
  sizes, strict search under the SHIPPED 5e8 default budget on the
  post-campaign engine, with multi-seed discipline at each edge (the I4
  lesson: cliffs are instance-dependent). Full data:
  benchmarks/results/2026-07-04-p2-envelope-probe.md.
- **New envelope (robust = places across seeds):** 9x9 16w -> 17w (18w noisy
  edge; 19w near geometric packing saturation - the generator itself
  struggles to plant a witness); 15x15 36w -> 42w (44w noisy edge); 21x21
  80w -> 82w robust, 84w noisy - the pre-campaign ">=82 hard for all" line
  is broken.
- **Cliff character:** instance-noisy AND search-bound. The seed-11 first
  failures stay not_proven at 4x budget (2e9) while sibling seeds at the
  same N place marginally at ~5e8 - the classic heavy-tailed backtracking
  distribution. RAISING THE BUDGET IS NOT THE LEVER; the right future play
  for the near-cliff regime is opt-in randomized restarts on not_proven
  (docs/research/arrange-search-algorithms.md #3), which is a product
  feature, not a search-speed change.
- **Adopted:** three envelope-guard rungs (heavy tier) in workloads.pl +
  baseline.json: ladder_09x09_17w (~46.7M), ladder_15x15_40w (~13.5M),
  ladder_21x21_82w (~10.8M). Near-cliff marginal instances deliberately NOT
  promoted (they'd flip on any solver change).
- **Also notable:** near the envelope, cost is instance-dominated, not
  word-count-dominated (15x15: 37w costs 23M, 42w only 4.2M under one seed).

---

## Campaign wrap-up — 2026-07-04

Seven experiments (5 shipped, 2 rejected), two research sweeps, two
measurement probes. Search inferences on the original ladder: **-67% to
-76% per rung** (15x15_36w: 182.6M -> 43.9M, ~9.5s -> ~2.3s wall), all
output byte-identical throughout. Shipped: E-H1 corner dedup (~2x), E-H2
var-cell grid, E-H3 tabled crossing/letter memos, E-H5 count_upto2, E-H7
pw/8 records. Rejected with recorded reasoning: E-H6 wdeg tie-break, E-H8
letter bitmask. Probes: P1 (killed CBJ before it was built; found the 80w
rung backtrack-free), P2 (envelope above).

State of the idea space at pause (why the campaign stopped here):
- Tree-size: CLOSED for the deterministic search (P1: no backjump
  headroom; E-H6: failure-concentration signal exists but is not
  convertible by tie-breaking; pre-pivot I2/E4/I6 cover value ordering and
  static orders).
- Constant-factor: MINED OUT (E-H7 verdict on representation; E-H8 verdict
  on prefilters; R1 checklist of confirmed non-hot-spots).
- Remaining levers are PRODUCT features, not search speed: (1) randomized
  restarts on not_proven for the heavy-tailed near-cliff regime; (2) WASM
  build guidance + one-time ladder inference-parity validation under
  swipl-wasm (docs/research/swi-vm-wasm-performance.md); (3) greedy-path
  (best-effort/candidates) benchmarking + incrementality, currently
  unmeasured (E-H2 left a known greedy micro-regression: placement_key now
  runs for every crossing candidate).

### E-H9 — check-only counting legality + O(1) boundary grid — SHIPPED (large)

- **Where:** core.pl (mrv_count_goal -> check_word_fits/5 + check_letters/6;
  assign_word/10 -> gs/2 bundle + no_word_merge_bg/2 + mark_boundary/5;
  init_gs/2 exported; no_word_merge/3 and placed_boundary_cell/3 deleted),
  arrange.pl (seed_from_fragment / greedy_construct use init_gs;
  fragment_conflict and word_best_placement destructure the bundle).
  Merge of experiment/e-h9-boundary-grid (e91ef6e).
- **What (two composed changes):** (1) the Cap=2 counting goal calls
  check_word_fits/5, a non-binding twin of assign_word (same checks, same
  order: Start>=1, check_prev_cell, per-cell letter-match-or-adj_is_free,
  check_next_cell, boundary test) that binds no grid var and builds no
  Cells/pw record - the parts assign_word materialized only for count_upto2
  to discard. Sound because a run's cells are distinct and adj_is_free reads
  only perpendicular off-line cells, so nothing a check would bind is
  re-read by the same check. (2) the no_word_merge placed-words scan
  (O(|Placed|) per placement attempt AND per counted candidate) replaced by
  a second grid-sized var term BoundaryGrid, bundled as
  gs(LetterGrid, BoundaryGrid): assign_word marks a placed word's 0-2
  on-grid boundary cells (atom b, trail-undone), the test is arg/3+nonvar
  per cell. Boundary cells stay UNBOUND in the letter grid (they must read
  as empty to the adjacency checks - the sentinel-in-letter-grid trap). The
  counting chain never inspects the grid, so the bundle threads with zero
  signature churn; fragment pinning and the greedy constructor share the
  same assign_word, so maintenance is uniform.
- **Attribution (measured separately):** Part 1 alone -0.04% to -0.81% per
  rung (monotonic, 12/12). Part 2 is the win. A Part-1 variant replacing
  no_word_merge's memberchk with arithmetic run-membership REGRESSED
  (+0.7..+2.4%) - building the short Cells list and memberchk-ing it beats
  per-placed-word arithmetic; kept the list.
- **Result (composed, --heavy, all 12 rungs):** -12.9% to -50.3%
  inferences, 0 regressions (21x21_80w 9.03M -> 4.49M, -50.3%; 15x15_36w
  43.9M -> 38.3M; envelope rungs 17w -17.5%, 40w -22.7%, 82w -36.8%). Wall
  -12..-35% on the heavy tail; RSS flat. On-model with P1: the
  backtrack-free, counting-dominated rungs gain most.
- **Identity verification (beyond byte-identical goldens, 204 plunit):**
  arrange output byte-identical on all 12 rungs vs pre-change code;
  FULL-TREE all_crossword solution counts identical
  (baseline/mrv_capped/mrv_inc x 4 corners); --seed and --best-effort
  outputs identical. Three new plunit tests lock init_gs shape,
  boundary-mark maintenance/edge omission, and check_word_fits/assign_word
  agreement + probe purity.
- **Caveats:** assign_word (exported) now requires the gs/2 bundle and
  ignores its PlacedWords arg; the boundary invariant rests on "every
  placement goes through assign_word" (true for search/fragment/greedy).
  Greedy findall snapshots copy two grids per kept candidate (path
  unbenchmarked, but it also gains the O(1) merge test).
- **Verdict:** SHIPPED, ratcheted. **Correction to the wrap-up above:** the
  "constant-factor mined out" claim was scoped too broadly - it was true of
  the REPRESENTATION seam (E-H7) and pair prefilters (E-H8), but this O(P)
  -> O(1) complexity cut on the per-node counting work sat outside both.
  Lesson recorded: distinguish "this seam is mined out" from "per-node cost
  is minimal"; P1's finding that counting dominates the backtrack-free
  rungs was the pointer, and it was still pointing after E-H8.

### Bench-tool fix (found during E-H9 acceptance): --record dropped new rungs

`check_baseline.pl --record` mapped over EXISTING baseline pairs only, so a
measured rung absent from baseline.json was silently not persisted while
the report printed "new rung, recorded". The three P2 envelope guards were
therefore never actually ratcheted (they lived only in history.jsonl).
Fixed: run_arrange rows now carry tier/warmup/budget/words; the recorder
builds a complete spec for first-seen rungs. All 12 rungs now gate.

### P3 — post-E-H9 envelope re-check — MEASUREMENT (envelope unchanged)

- **What:** re-ran all 15 of P2's marginal/failed near-cliff instances at
  the shipped 5e8 budget on the post-E-H9 engine (fixtures regenerated,
  byte-verified against P2's). Full data:
  benchmarks/results/2026-07-04-p3-envelope-recheck.md.
- **Result:** NOTHING flipped. No hard failure placed; the budget-marginal
  placers reproduced their counts to within ~10 inferences, still pinned at
  the ceiling. Only already-comfortable placements got cheaper (-18..-20%,
  in E-H9's ladder range). Envelope statement identical to P2.
- **The decision-relevant finding:** node cost and envelope are DECOUPLED.
  P2 showed more budget is not the lever (fails persist at 4x); P3 shows
  cheaper nodes are not the lever either - the campaign's biggest single
  win moved the near-cliff envelope by zero words. Mechanism: a
  search-bound corner consumes its full budget slice regardless of
  per-node cost, so only within-budget solutions take the discount.
  Constant-factor wins buy LATENCY (and WASM headroom); envelope gains
  need a DIFFERENT TREE per attempt -> randomized restarts on not_proven
  remain the right (product-level) play.

### E-H10 — watched-witness recount skipping — REJECTED (unsound as briefed; sound residue fails ratchet)

- **Idea:** cache each word's first up-to-2 witness placements; on the
  shares-letter recount trigger, revalidate witnesses via check_word_fits
  (O(len)) and keep the previous bucket if they survive, skipping the full
  rescan. Branch experiment/e-h10-watched-witnesses (92179c8), not merged.
- **The core finding - the brief's soundness argument was WRONG:** the
  orchestrator's monotonicity premise ("a word's placement set only shrinks
  within a descent") is false for exactly the triggered words: sharing a
  letter with the just-placed word is when a count can RISE (new crossings
  against that word - the documented reason the trigger IS letter-sharing).
  A bucket-1 word with a surviving witness can have true count 2; keeping
  bucket 1 under-counts and reorders the search. Empirically: the bucket>=1
  variant kept the full-tree solution SET identical but changed the emitted
  first layout on 11/12 rungs. Only Count==2 is sound (the cap saturates:
  both witnesses valid => still >=2).
- **Sound Count==2 version (byte-identical everywhere, 204 plunit):**
  bench FAIL - 3 regressions / 6 wins. 21x21_80w -18.5%, 21x21_82w -6.0%,
  15x15_32w -3.4%, but light rungs +0.8..+1.4% from the flat tax (wc-unwrap
  in the ordering loop + capture in every recount). Witness-survival
  instrumentation explains the split exactly: the shortcut fires on 35% of
  triggers at 80w but ~5-6% on light rungs (hit rate tracks fan-out). The
  sound scope also excludes the volume: count<=1 is the plurality bucket
  (P1) and a bucket-2 rescan already early-exits at the 2nd hit.
- **Verdict:** REJECTED as shipped. One focused variant remains (E-H10b):
  witnesses in a SPARSE side-map holding only cap-saturated words so the
  count map stays bare ints (no ordering-loop unwrap) and capture cost
  moves off the light path - aiming to keep the 21x21 win with zero
  regressions. If that still fails the ratchet, this avenue closes.

### E-H10b — sparse witness side-map — REJECTED (watched-witness avenue CLOSED)

- **What:** E-H10's sound Count==2 shortcut rebuilt so the light path is
  byte-for-byte master: CountMap stays bare ints, witnesses in a
  backtrackable global (b_setval; the briefed state-threading variant
  regressed +1.1..+3.1% from per-word wrap/unwrap and was dropped), capture
  confined to the recount fallback behind a previously-saturated gate.
  Branch experiment/e-h10b-sparse-witnesses (0ae0ec9), not merged.
- **Result:** did its narrow job - E-H10's structural taxes gone (36w/17w
  now pass), full 21x21 win retained (80w -17.6%, 82w -5.6%), identity
  verified to the E-H9 standard (byte-identical on all 12 rungs, full-tree
  counts identical across 4 strategies x 4 corners). But 3 light rungs
  still regress (+0.6..+1.5%) past the 0.5% ratchet.
- **The decisive diagnostic:** with shortcut AND capture disabled, the
  bare admission probe (one assoc lookup per shares-letter word) alone
  costs +0.60..+0.74% on the four smallest rungs - the eligibility test
  exceeds tolerance before the shortcut does any work. On light rungs the
  hit rate is ~5-6% and a bucket-2 rescan already early-exits at the 2nd
  hit, so there is nothing to bury the probe under. No size-independent
  signal can skip the probe only where it cannot pay off; input-size
  switches are off-charter.
- **Verdict:** REJECTED; the watched-witness avenue is CLOSED as
  engineering. What remains is a POLICY question, not an experiment: the
  blocked trade is ~+400 absolute inferences on a 28k rung vs -790k on a
  4.5M rung. If the ratchet ever adopts an absolute-floor tolerance
  (e.g. regressions under N absolute inferences do not gate), E-H10b is
  shippable as-is from its branch.
- **Policy decision (2026-07-05, Ned):** ratchet stays strict (relative
  0.5%, no absolute floor). Rationale: no product evidence yet that dense
  21x21-class arranges will be common, and the win is only perceptible in
  that scenario (under WASM); absent evidence, ratchet discipline is worth
  more than an unproven latency gain. REVISIT TRIGGER: if browser telemetry
  or product direction shows dense large-grid arranges are a real usage
  class, adopt an absolute-floor tolerance (~10k inf) and merge
  experiment/e-h10b-sparse-witnesses (0ae0ec9) as-is - it was verified
  byte-identical and its win scales with word count.

## Fill campaign — 2026-07

Plan: docs/plans/fill-perf-campaign.md (red-teamed GO 2026-07-05). Bench:
benchmarks/run_fill.pl; ratchet: check_fill_baseline.pl + fill_baseline.json
(search_inf gates at 0.5% relative; load_inf reported, non-gating until
Phase 3 decides); equivalence oracle: check_fill_identity.sh (CLI
byte-identity on every rung).

### Phase 0 (fill) — benchmark build-out — COMPLETE (merged eeb2611)

- **Where:** benchmarks/{run_fill,fill_subjects,fill_workloads,
  check_fill_baseline}.pl, check_fill_identity.sh, gen_fill_dict.py,
  fixtures/dict/ (ENABLE 172,823 words, public domain, + seeded 10k/25k/50k
  subsets, all byte-frozen), fixtures/fill_grid_{04,05,09,11,13,15,17,21}a
  .json + fill_seed_11a.json, fill_baseline.json + fill_history.jsonl +
  fill_identity.sha256, tests/golden/fill_15_bench.json,
  benchmarks/results/2026-07-05-fill-phase0.md. Branch campaign/fill-phase0
  (base 9751190), fast-forwarded to master. NO product-code change.
- **What:** the fill product bench (four attribution buckets — command /
  dict_load / grid / search; fresh slots rebuilt OUTSIDE the timed goal per
  sample, since fill_attempt destructively unifies the shared cell vars), an
  11-rung completing ladder, a forked ratchet, a CLI byte-identity oracle
  (the Phase 2 gate), and a 15x15-vs-full-ENABLE scale golden.
- **Ladder:** 11 rungs, search 0.53M -> 34.9M inf, all `filled`
  deterministically, min headroom 57x under the 2e9 bench budget. Baseline
  verified IN-FILE after --record (all 11 rungs, complete specs — the
  arrange silent-drop bug did not recur). search_inf and load_inf reproduce
  +0.00% across 3 full runs; cold/warm search delta is a fixed 319 inf
  (arrange's was ~7.6k).
- **Headline finding (decision-grade):** dict_load exceeds 50% of
  end-to-end CLI latency on 10 of 11 rungs (75-84% on full ENABLE; a flat
  ~3.25s / 26,602,930 inf regardless of grid; 1.6M inf at 10k words). Only
  g17_50k is search-dominated. The plan's Phase 3 conditional FIRES:
  startup/dict-load outranks per-node search work for perceived latency.
  P-F1 is therefore weighted toward the load layer (where inside
  load_dict/build_index do the 26.6M go?) alongside the search-internals
  attribution.
- **Surprise (P-F1 must explain):** search cost is INVERSELY related to
  dictionary size in the completing regime — sq04 costs 0.53M at 172k words
  but 7.7M at 50k (14.6x); 9x9 x full ENABLE (34.9M) is the ladder's
  hardest rung while 21x21 x full is 3.3M. Hypothesis: a bigger dictionary
  makes early candidates succeed (shallow tree); a thinner one forces deep
  backtracking. If confirmed, per-node candidate-cost wins (F-H3/F-H5) pay
  most exactly where trees are deep — the thin-dictionary/hard regime.
- **Envelope:** curated blocked_15a x full ENABLE -> not_proven at a 2e8
  budget (~6x the hardest completing rung): full-width 15-cell lights
  over-constrain the interlock. The bench-authored masks (shorter lights)
  complete; the ladder ceiling (~35M) is honest, not forced. The plan's
  ~100M ladder target sits inside the non-completing regime on real grids.
- **Fixture trap found (latent product bug, out of campaign scope):** a
  seed pin whose answer is ALSO reachable from the dictionary can crash the
  emit — seeded slots are exempt from the no-duplicate rule (fill.pl:196),
  so the search may re-place the pinned word elsewhere and
  assign_clue_numbers throws unique_key_pairs (CLI exit 1 with an internal
  error, e.g. pinning AAH against ENABLE). The seed fixture avoids the
  overlap. Deserves a product fix ticket (clean report, not a throw);
  output-semantics territory, so not this campaign's to change.
  **RESOLVED (2026-07-06, fix/fill-seed-pin-crash;
  docs/plans/fill-seed-pin-crash-fix.md):** the throw site was
  answer_meta_assoc/2 (core.pl), not assign_clue_numbers. Seed answers now
  pre-seed the search's Used set (fill.pl seed_used/3 — a correctness fix:
  the search takes the non-duplicate alternative), duplicate seed answers
  are rejected up front (fill_seed_duplicate, hooked), and fill's emit
  boundary re-runs check_unique_answers/1 as defense in depth. The
  g11_full_seed baseline was re-recorded (correctness-driven).
- **Orchestrator decision (2026-07-05):** the plan's pre-campaign
  "--budget CLI hook or reconcile 800M/500M" fix is DEFERRED: the bench
  drives the budget-explicit fill_attempt/8 in-process, so no CLI surface
  is needed for measurement, and moving the shipped 800M default can flip
  not_proven outcomes (product semantics). REVISIT TRIGGER: envelope/
  restart product work, or any Phase 3 change that touches CLI startup.
- **Process note:** the original Phase 0 build agent stalled twice
  (600s dead stream during bench invocations); recovered by a fresh agent
  over the same worktree with commit-first ordering and swipl hygiene
  (`</dev/null` on stdin so a load error cannot park swipl at the
  interactive toplevel; timeout prefixes; one rung per invocation). The
  inherited harness itself was sound — it reproduced the orchestrator's
  independently measured load counts exactly.

### P-F1 — fill attribution probe (load stages, search internals, backtrack shape) — MEASUREMENT

- **Where:** branch probe/p-f1-attribution (base a178923, commit 0329b2f;
  never merges): benchmarks/probe_f1/*.pl +
  benchmarks/results/2026-07-05-p-f1-attribution.md. NO product-code
  change; check_fill_baseline PASS +0.00% on all 11 rungs (core and
  --heavy) with the probe modules loaded beside the engine.
- **Load:** the normalize parse loop is 72.9-76.0% of load_inf at
  10k/50k/172k (the pre-registered "build_index >60%" REFUTED at
  19.4-22.8%); stage sums reconcile to the whole at a flat -121 inf
  (<=0.008%). Inside the loop, the include/3 yall lambda
  ([C]>>char_type(C,alpha), fill.pl:131) costs ~19.2M of ENABLE's 26.6M —
  a named-helper rewrite is output-identical, ~-53% load_inf, ~-1.2s wall.
  Inference-BLIND C walls found: build_index keysort 1.15s at ~0 inf;
  ~1.2s GC residue visible only in the un-staged whole.
- **Waste:** the ladder grids (lights 3-5) cannot use 92-98% of ENABLE
  words / 96-99% of index triples, but realistic 3-15-light grids waste
  only ~2.5%/4.6% (lengths >=16) — slot-length filtering is
  ladder-inflated; size it against real grids, not the bench masks.
- **Search:** ord_intersection counting is 59-90% of search_inf on all
  four probe rungs; select_mrv's full per-node recount is 79.5-85.6%;
  candidate MATERIALIZATION is 0.24-0.63% — the plan's ranked gap #1
  REFUTED (its rank came from a wall-clock micro-bench measuring a real
  cost in the wrong currency for where it sits in the tree);
  member/memberchk/unify <=0.5%. Equivalence: term-identical output on
  every instrumented run; profiler port counts match the counters exactly;
  two misattribution classes caught (LCO hides the ordset walk from
  self-time; inference counts hide length/2 = 28.5% of g17_50k profiled
  wall). Wall medians stay FIRST-CLASS beside the inference gate.
- **Backtrack:** the Phase-0 inversion is a node-count effect: thin
  dictionaries multiply nodes faster than they shrink per-node cost
  (sq04 47.5x nodes x 0.31x inf/node = 14.6x; g17 17.1x x 0.45x = 7.7x).
  Refinement of the hypothesis: abundant dictionaries shrink trees but do
  NOT make them backtrack-free — g09_full churns 9,961 nodes with 99.7%
  of placements unwound in a depth-17-24 thrash band. Dead ends are
  discovered by the RECOUNT (candidate medians collapse to 0-1 after
  early depths), not by trying words — further concentrating value in the
  counting path.
- **Orchestrator decisions (2026-07-05, at adjudication):**
  1. Campaign order RESTRUCTURED: Phase 3 leads (dict_load is 58-84% of
     end-to-end on 10/11 rungs). First experiment F-P3-L1: named-helper
     rewrite of the normalize lambda (dispatched).
  2. Phase 2 order: F-H1 (incremental counts; attacks the 79.5-85.6%
     recount) -> F-H2 (bitset counting; attacks the 59-90% intersection;
     WASM bignum gate probe dispatched in parallel) -> F-H3 RESCOPED to
     bucket/index access infrastructure only.
  3. F-H4 and F-H5 CLOSED on mechanism: their targets (used-list scans,
     candidate materialization) measure <=0.6% of search_inf — no
     experiment can pay there. Revisit only if a later change makes
     materialization hot again.
  4. load_inf gating: decide promotion from reported to gated when the
     first Phase 3 win is recorded.
- **Process notes:** the probe's first worktree was provisioned 10 commits
  stale (pre-merge) and then deleted mid-session — the brief's base check
  caught it before any measurement (the arrange campaign's stale-worktree
  lesson paying off); recovered by cutting probe/p-f1-attribution from
  a178923. A transient untracked file (benchmarks/probe_p_f1.pl, a
  foreign alternative probe implementation) appeared in the new worktree
  mid-session and was gone by adjudication; never committed or used.
  The probe briefly switched the MAIN checkout's branch (restored, refs
  identical, zero content change) — future briefs say: never touch the
  main checkout.

### Correction + adoption of the P-F1 entry above — split-brain incident (2026-07-05)

- **The entry above (36432d2) was written by a second, concurrent
  orchestrator instance**, not the active one. At the 17:07 session
  restart, a background agent named "Check if agent is still running" was
  spawned as a fork-resume of the orchestrator's own transcript with auto
  permissions; having inherited full campaign context, it went beyond its
  status check and re-performed the orchestrator role: adjudicated P-F1,
  committed 36432d2 to master, and began creating an F-L1 experiment
  branch in the main checkout — concurrently with the active orchestrator
  doing the same adjudication from the same evidence. Ned attached to it
  via the background-agent view and paused it.
- **Adoption:** the active orchestrator ADOPTS 36432d2 after independent
  verification (commit scope re-checked: 7 bench-side files, no product
  code; results-doc numbers spot-verified; the +0.00% gate additionally
  witnessed by a second, independently implemented instrument — the
  re-dispatched P-F1 agent — whose shares matched to rounding before it
  discovered the collision and stood down). Both orchestrator instances
  reached the same verdicts independently, which is corroboration, not
  double-counting: both read the same probe.
- **Corrections to the entry above:** (1) its two "(dispatched)" claims
  were aspirational — no F-L1 or WASM-probe agent, branch, or worktree
  existed; the real dispatches follow this entry. (2) Its transient
  branch experiment/f-p3-l1-parse-lambda carried no unique work and is
  deleted; the real experiment runs as experiment/f-l1-normalize. (3) Its
  "foreign alternative probe implementation" process note describes the
  legitimate re-dispatched P-F1 agent's instrument, not an intruder.
- **Process rule (standing):** exactly ONE active orchestrator; ledger
  writes and dispatches are its alone. Dispatched agents return draft
  entries and never commit to master or touch the main checkout. Never
  fork-resume an orchestrator transcript into an autonomous background
  agent — a status-check prompt does not survive contact with inherited
  campaign context and auto permissions.
- Naming note: the campaign's load-layer experiments are F-L1 (normalize
  parse fix), F-L2 (precomputed wordlist/index artifact), F-L3
  (slot-length-filtered index build; deprioritized — P-F1 measured ~2.5-5%
  waste on realistic grids vs 96%+ only on the short-word bench masks).

### F-H2 gate probe — bignum AND+popcount under WASM — MEASUREMENT (accepted)

- **Branch/commit:** probe/f-h2-wasm-bignum @ 62dc2b4 (instruments stay
  branch-only, per P-F1 precedent); results doc copied to master:
  benchmarks/results/2026-07-05-f-h2-gate-probe.md. No product code, no
  baseline writes — verified at adjudication.
- **Verdict: GO on the counting kernel.** The plan's central fear — "the
  one change most likely to be a native win and a WASM loss" — is
  REFUTED. Backends confirmed: native = GMP 6, WASM = LibBF (swipl-wasm
  under node, post-setjmp-fix build). ord_intersection+length vs
  `A /\ B` + popcount: 438-5152x same-N under WASM (47-3881x native);
  on REALISTIC pairings from P-F1's measured set-size distribution the
  worst case is 9x (median-size set × full-ENABLE bucket, WASM) — never
  near 1, never inverted. LibBF's ~3-14x per-op penalty vs GMP is offset
  by the ordset VM path's own ~3.6-4.0x WASM slowdown. Red-team's
  250-600x native anchor reproduced (494x @ N=1000).
- **The binding condition (the probe's real finding):** bitset mask
  CONSTRUCTION costs ~18% of load_inf and 27% (native) / 36% (WASM) of
  load WALL at 172k — not "small" as pre-registered, and worse exactly
  where latency hurts most. Since load dominates end-to-end (P-F1:
  58-84% on 10/11 rungs), a naive load-time mask build spends the
  dominant term to speed the minor one. F-H2 is therefore GO only as:
  (i) masks shipped inside the Phase-3 precomputed-index artifact (F-L2)
  — construction amortized to ~0, pure win both runtimes; or (ii) scoped
  to search-dominated workloads (g17_50k/g09_full class). Never a naive
  load-time build.
- **Ratchet implication (orchestrator adjudication):** the bignum kernel
  is 3 inf/op flat vs 405-96005 — an inference-blind, wall-only win
  (same class as P-F1's keysort/length findings). When F-H2 is built,
  search_inf will collapse for the wrong-looking reason while the
  construction tax lands in load. The charter's gate stays inference-
  based: load_inf MUST be promoted to gating before or with F-H2
  (promotion decision already scheduled at first Phase 3 acceptance),
  with wall medians recorded as corroborating evidence. The probe's
  "gate on wall" suggestion is adopted in substance, not in mechanism.
- **Sequencing confirmed:** F-H2 after F-H1 (F-H1 cuts count volume;
  F-H2 cheapens each remaining count), design folded into the F-L2
  artifact decision.
- **Bonus findings:** inference counts verified bit-identical native↔WASM
  for these kernels (first direct evidence for the ladder's portability
  premise); docs/research/swi-vm-wasm-performance.md's LibBF dismissal
  reaches the right conclusion via a non-covering reason (rational vs
  integer arithmetic scope) — this probe supplies the integer-bignum
  evidence directly.
- **Process notes:** the probe's worktree was again cut from stale
  9751190 (session-start HEAD — now a known harness behavior, third
  occurrence); the agent self-corrected onto 397615d before the
  orchestrator's fix arrived and re-verified base. Separately, the
  resumed F-L1 runner's shell cwd reset to the MAIN checkout and its
  branch creation briefly switched the main checkout off master (caught
  with zero content damage; agent relocated to its worktree; master
  restored). Standing rule extended: resume/fix-up messages must repeat
  the absolute-path + never-main-checkout discipline, not assume the
  original brief's context survives a session resume.

### F-L1 — normalize parse loop: yall lambda → first-order filter — ACCEPTED (merged 16452de)

- **Where:** branch experiment/f-l1-normalize @ 1b3220e (base ed5d556),
  merged to master 16452de. prolog/crosswordsmith/fill.pl:
  normalize_word/2 now calls new alpha_chars/2; results doc
  benchmarks/results/2026-07-05-f-l1-normalize.md.
- **What + mechanism (verified independently by the runner):** replaced
  `include([C]>>char_type(C,alpha), Cs, L)` with a first-order
  alpha_chars/2 recursion making the identical char_type(C, alpha)
  decision per character (string_upper path untouched). fill.pl never
  imports library(yall), so the lambda was runtime-meta-called with a
  copy per character: 12.33 inf/char, 19.36M of ENABLE's 26.6M load_inf
  — confirming P-F1's ~19.2M attribution. Control experiment: WITH yall
  imported the lambda compile-expands and costs helper-price; the
  premium was purely the runtime meta-call. The per-char cost model
  predicted the product's load_inf to 1 inference (12,468,069 predicted
  vs 12,468,070 measured for V1; 10,724,706 vs 10,724,707 for V2).
- **Result (composed == standalone; master had only gained docs over the
  branch base):** load_inf 26,602,930 → 10,724,707 (−59.69%) on the 9
  full-ENABLE rungs; 7,757,017 → 3,164,117 (−59.21%) on the 50k rungs.
  search_inf +0.00% on all 11 rungs (bit-identical, verified in the
  recorded baseline diff). Command wall −12% to −30% per rung (warm
  in-process load_dict 3.47s → 1.96s); RSS ±0.5%. Identity oracle 11/11
  OK; full test suite green on branch and on merged master (orchestrator
  re-ran all three gates independently). Baseline re-recorded at
  16452de; history appended.
- **Variants:** V1 (named helper inside include/3 — the P-F1 −53%
  instrument variant) measured 12,468,070 @172k: real but strictly
  slower than V2's first-order loop (10,724,707, also drops include/3's
  per-element call/3). V1 rejected. No V3: avoiding the intermediate
  char list would leave the frozen string_chars path — semantic risk
  for a non-dominant payoff.
- **Provenance:** V2 originated as the twin orchestrator's uncommitted
  14-line draft (preserved at its worktree teardown). Treated as a
  hypothesis per the never-partially-trust rule: the runner re-derived
  equivalence by fuzz (42 edge strings incl. ß→SS, combining marks,
  zero-width/NBSP, non-Latin digits, ligatures; 20k random char lists
  raw+uppercased; every line of all four frozen dicts, elementwise ==
  against the frozen include/lambda reference) plus a determinism check
  before shipping.
- **Expectation vs outcome:** pre-registered −53% was V1's number; the
  shipped V2 beat it at −59.7% for the pre-stated reason (include/3
  overhead also deleted). g17_50k's smaller wall gain (−11.9%) is the
  expected signature of the one search-dominated rung.
- **Follow-on signal:** remaining load wall is now build_index-bound
  (~1.15s C keysort + GC residue, inference-blind) — strengthens F-L2
  (precomputed index artifact) as the next Phase 3 experiment and
  the F-H2 mask-shipping vehicle.
- **Ratchet policy (scheduled decision, taken now):** load_inf is
  PROMOTED from reported-only to GATING at the same 0.5% relative
  tolerance, effective with the baseline recorded at 16452de. Grounds:
  load is the dominant end-to-end term and now carries real wins worth
  defending; F-H2's measured construction-tax risk lands exactly in
  load; determinism is proven (runner's 3× md5-identical count tables +
  orchestrator reproduction). Implemented as a separate infra commit.

### F-H1 — incremental candidate counts for MRV — ACCEPTED (merged 9b1d34c)

- **Where:** branch experiment/f-h1-inc-counts @ 1724d8f (base d38a649),
  merged to master 9b1d34c. prolog/crosswordsmith/fill.pl: fill_search/4
  now drives fill_search_inc/4 (carried cnt(Count,Slot) list +
  select_min_count/remove_slot/newly_bound_cells/recount_crossing);
  select_mrv/6 retained unchanged as the reference selector (white-box
  test pins R6/P13 call it directly; it documents the selection rule the
  incremental drive reproduces). Results doc
  benchmarks/results/2026-07-05-f-h1-inc-counts.md.
- **What + mechanism:** per-slot candidate counts carried as backtrack-
  restored threaded state; after each placement, recount EXACTLY the
  slots crossing a newly-bound cell (measured crossings/slot: min 3 /
  mean ~3.8 / max 5, vs slot lists up to 154) — replacing select_mrv's
  per-node full recount (79.5-85.6% of search_inf per P-F1). Tree
  preserved node-for-node: counts exact by construction (candidate_count
  is a pure function of cell binding state and never reads Used —
  try-time only), identical min c(Count,Start,Dir) selection (Start+Dir
  unique ⇒ total order, no genuine ties), same winner materialization,
  completed 0-count slots retained so dead branches fail exactly as
  before. Affected detection by ground cell NUMBER identity — no shared
  crossing variable inspected or copied (the findall trap structurally
  avoided).
- **Result (composed == standalone; master did not move between base and
  merge):** search_inf down on ALL 11 rungs, −21.66% (sq04_50k) to
  −57.37% (g09_full, the top rung 34,880,750 → 14,870,446); deep/dense
  rungs cut hardest (g21_full −43.96%, g11_full_seed −40.01%). load_inf
  +0.00% on all 11 (gated since F-L1 — first experiment to pass under
  the promoted gate). cmd_wall down on every rung; g17_50k −40.65% wall
  vs −29.86% search_inf — the predicted extra inference-blind win from
  deleting C-native length/2 bucket walks. Identity 11/11 OK, 204 tests
  green, 3× reproduced (runner) + independently re-run by the
  orchestrator on the branch (all three gates) and on merged master.
  Baseline re-recorded at 9b1d34c.
- **Expectation vs outcome:** pre-registered −30..−70% on deep rungs;
  g09_full −57.4% and g21_full −44.0% in-band, g17_50k −29.9% at the
  edge. Win correctly bounded below the naive ~85% by the untouched
  winner-path candidates/4 (12.7-17.8% per P-F1) + the O(|Rest|)
  affected-detection scan.
- **Anomaly (accepted, ledgered):** cmd_rss +82% on g21_full ONLY
  (278 → ~507 MiB, stable across runs, reproduced by the orchestrator).
  Correlates with backtrack DEPTH (g21_full depth 153; g09_full depth 31
  does not move): per-placement rebuild of the counted-list spine raises
  global-stack high-water under deep backtracking. Bounded, not a leak;
  RSS is reporting-only by Phase 0 design. Standing mitigation if the
  Phase 3/WASM posture needs footprint: (Start,Dir)-keyed assoc carry
  (O(log n) path-copy) — a redesign requiring full re-validation; do NOT
  fold it into an unrelated change.
- **Variant:** cnt/2 term-reuse for unaffected slots — search_inf
  byte-identical, RSS unmoved; reverted for the simpler diff.
- **Phase 2 status:** F-H1 done. F-H2 (bitset counting, gate GO) remains
  scoped to the F-L2 artifact design per the gate probe's binding
  condition. F-H3 rescoped to infra; F-H4/F-H5 closed on mechanism.

### F-L2 — precomputed index artifact — ACCEPTED (merged 540a2ce)

- **Where:** branch experiment/f-l2-index-artifact @ c23ec1c (3 commits,
  base f9e5a83), merged to master 540a2ce. prolog/crosswordsmith/fill.pl
  (fill_save_index/2, fill_load_index/4, fill_solve_index/5; fill_solve
  refactored into shared fill_prepare + fill_place_and_emit so both
  modes run the identical engine); CLI `fill --save-index FILE` (build)
  and `fill --index FILE` (consume), mutually exclusive, raw path
  unchanged; new artifact-mode oracle
  benchmarks/check_fill_identity_artifact.sh; +8 plunit tests (212
  total). Results doc benchmarks/results/2026-07-05-f-l2-index-artifact.md.
- **What:** the dictionary parse + index build is a pure function of the
  frozen dict file; serialize the EXACT runtime terms (DictByLen + assoc
  Index — AVL shape preserved, loaded == built proven at 10k/50k/172k)
  once via `fill --save-index`, and load them back in ~0.1s instead of
  recomputing ~2s per invocation. Artifact term fill_index(Version,
  Meta, DictByLen, Index); Meta carries dict_sha256 / swi_version /
  provenance and is the F-H2 extension point (masks = new Meta key under
  a Version bump). Loader REFUSES (clear fill_index_* errors, no silent
  rebuild) on missing/unreadable/malformed/version/SWI/hash mismatch.
- **Format bake-off (measured at 172k, warm):** fastrw 0.099s / 4 inf /
  14.6MB beat .qlf-of-fact 0.269s / 710 inf / 15.5MB; build cost fastrw
  0.2s vs qcompile 9.3s; no qlf clause-size limit at 172k. SHIPPED
  fastrw; .qlf documented as the WASM-aligned alternative (localized
  swap, schema unchanged).
- **Result:** warm artifact load ~20x the post-F-L1 raw load at 172k
  (1.96s → 0.099s). End-to-end CLI (raw → --index): sq04_full 2530→230ms
  (11x), g11_full_seed 2520→240ms, g17_full 2620→410ms, g21_full
  2850→500ms, g09_full 3220→870ms (3.7x — search-dominated, least gain,
  as expected). Orchestrator smoke: g11a full-ENABLE 2.70s → 0.40s,
  output byte-identical (cmp).
- **Gates:** raw path +0.00% on BOTH gated layers, all 11 rungs
  (baseline intentionally UNMOVED — nothing to ratchet; post-merge
  history row logged at 540a2ce with counts identical to baseline).
  Identity 11/11 in BOTH raw and artifact modes against the same pinned
  fill_identity.sha256 (the cross-mode equivalence gate). Tests green on
  branch and merged master. All independently re-run by the orchestrator.
- **Honest-metric note:** artifact-mode load_inf ≈ 4 (C-level fastrw) —
  WALL is the metric for artifact mode; the inference collapse is a
  measurement blind spot, not the win. The 11 raw rungs remain the
  inference-gated floor.
- **Risks:** fastrw's binary format is SWI-version-bound → the artifact
  is a build-time cache, not a distribution format; embedded swi_version
  + refusal handles it. Browser milestone must verify fastrw under
  swipl-wasm or emit .qlf instead (schema unchanged either way).
- **Deferred proposals (results doc §h, not adopted):** artifact-mode
  bench-manifest rungs and an artifact file-size drift tripwire — revisit
  if/when F-H2 masks land in the artifact (schema v2).

### F-H2 — bitset counting via artifact v2, masks opt-in — ACCEPTED (merged 2cb7c0f)

- **Where:** branch experiment/f-h2-bitset-count @ a2248bd (build) +
  0aa13ca (adjudicated follow-up: masks opt-in), base 7181337, merged
  2cb7c0f. fill.pl: candidate_count/5 dispatch (none → ordset kernel
  unchanged; masks(_) → bignum `/\` chain + popcount), Masks threaded
  through fill_search/fill_search_inc/recount_crossing as pure state;
  fill_attempt/8 kept as the pristine ordset seam (bench/test-pinned);
  artifact schema v2 with masks as an OPTIONAL Meta key; CLI `--masks`
  on --save-index (default emits NO masks). Phase A instruments in
  benchmarks/probe_fh2/ (merged — they document the attribution method);
  results doc benchmarks/results/2026-07-05-f-h2-bitset-count.md.
- **Phase A gate (pre-registered 20% rule):** on the post-F-H1 engine the
  counting kernel is 37.8% (g21_full) / 44.1% (g17_50k) / 47.3%
  (g09_full) of fill-phase wall — decisively past the threshold; build
  warranted. (P-F1's 59-90% share was correctly treated as stale.)
- **Result:** masks mode cuts fill-phase wall −24..−33% on the hard
  rungs (3× reproduced; counts byte-exact); end-to-end −11% (g09_full) /
  −14% (g17_50k). search_inf collapse is only 1.6-1.8x and is NOT the
  win (inference-blind kernel — wall is the metric, per the standing
  note). Default path: v2-no-masks artifact BYTE-IDENTICAL in size to
  v1, CLI wall equal to ≤1ms interleaved-paired on g21_full/g11_full.
  Raw path +0.00% both gated layers all 11 rungs; identity 11/11 in all
  THREE modes (raw / artifact-default / artifact-masks); 214 tests
  green. All independently re-run by the orchestrator; baseline
  intentionally unmoved (post-merge history row at 2cb7c0f).
- **Scoped correction to the F-H2 gate probe verdict above:** "masks
  shipped in the artifact → construction amortized to ~0 → pure win
  everywhere" is HALF right: construction amortizes as predicted, but
  the mask SIZE (+32% artifact) does not — every load pays it (~+0.09s
  at 172k native, measured in-process fast_read +0.031s vs fill −0.021s
  on g21). Hence masks are opt-in for search-heavy use, never
  default-on. The gate's GO stands; its deployment clause is corrected.
- **Measurement integrity (runner self-caught, on the record):** (1) a
  post-commit git-stash A/B cycle silently measured version-refusal
  errors as "v1 timings" — discarded, replaced with a detached
  measurement worktree + per-sample exit/output checks; (2) the
  originally reported "+20% g21 end-to-end regression" was ~5x inflated
  by ~100ms CLI-timer quantization; contemporaneous interleaved pairing
  puts the true masks-mode cost at +2-3% on g21. Direction unchanged,
  magnitude corrected in the results doc; the opt-in decision predates
  and survives the correction.
- **WASM note:** the gate probe measured the kernel ratio ≥9x under
  WASM/LibBF at realistic sizes; masks-on is the recommended default for
  the browser milestone IF fastrw artifacts verify under swipl-wasm
  (else emit .qlf; schema unchanged).

### Fill campaign — CLOSE-OUT (2026-07-05)

- **Executed to plan:** Phase 0 (bench substrate, 11-rung ladder,
  byte-identity oracle, two-layer inference ratchet), Phase 1 (P-F1
  attribution), both do-not-skip gates (P-F1 before experiments; WASM
  bignum probe before F-H2), and four accepted experiments — F-L1,
  F-H1, F-L2, F-H2 — every one output-byte-identical, independently
  verified at adjudication, and ratcheted. F-L3 deprioritized on P-F1
  evidence (~2.5-5% at stake); F-H3 rescoped and absorbed; F-H4/F-H5
  closed on mechanism (≤0.6% targets).
- **Composed end state on the hardest rung (g09_full, 15x15-class
  search-bound, full ENABLE 172k):** search 34,880,750 → 14,870,446 inf
  (−57.4%); end-to-end CLI 5.18s (campaign start) → 3.23s raw → 0.87s
  artifact → ~0.77s artifact+masks. Load-dominated rungs: ~2.5-4s →
  0.23-0.50s (up to 11x). Warm dictionary load at 172k: 3.47s →
  1.96s (F-L1) → 0.099s (F-L2 artifact). Both ratchet layers gated at
  0.5% defend the floors.
- **Standing residue for future work:** (1) g21_full RSS depth anomaly
  (+82%, F-H1 entry; assoc-carry mitigation specified); (2) browser
  milestone: verify fastrw under swipl-wasm or emit .qlf, run the
  ladder-parity check, then consider masks default-on there; (3)
  artifact-mode bench rungs + size tripwire (F-L2 §h proposals);
  (4) latent product bug outside campaign scope: seed+dict overlap can
  crash unique_key_pairs (found during Phase 0; needs its own ticket).
- **Process record:** one split-brain incident (duplicate orchestrator
  from a fork-resumed transcript — contained, adopted-with-corrections,
  standing one-orchestrator rule); a recurring harness bug (isolation
  worktrees cut at session-start HEAD) neutralized by base-check-first
  briefs after catching it three times; the F-H2 runner's self-caught
  A/B contamination validates the measure-against-a-pinned-worktree
  rule. Command hygiene that held: `</dev/null` + timeout on every
  swipl call, absolute paths, agents never touch master or the main
  checkout, verify state not success messages.


---

## Appendix: full experiment-time writeups (landed 2026-07-06 from record branches)

The close-out entries above summarize E-H8/E-H10/E-H10b; the sections below are
the full writeups as committed on the (now-deleted) experiment branches at the
time the work was done. The rejected implementations themselves are preserved as
patches under `docs/experiments/` (e-h8-bitmask-prefilter.patch,
e-h10-watched-witnesses.patch, e-h10b-sparse-witnesses.patch,
p1-backtrack-instrumentation.patch).

### E-H8 — per-word letter bitmask prefilter — REJECTED

- **Idea (code present in the experiment worktree, do-not-merge):** short-circuit
  the two "do these words share a letter?" hot paths with an O(1) integer AND on
  a per-word letter-presence bitmask, instead of a list scan.
  - `find_intersecting_word/6`: compute the candidate word's mask once, then gate
    each placed word by `CandMask /\ PMask =\= 0` before the `intersection/3` +
    `list_to_set` scan (skip pairs that share no letter, hence no crossing cell).
  - `inc_count/10` (the `mrv_inc` recount gate): replace the exact
    `shares_letter/2` member scan with `EMask /\ LastMask =\= 0`.
- **Bit mapping:** bit index = `CharCode /\ 63` (`1 << (Code /\ 63)`). a-z→33..58,
  A-Z→1..26, 0-9→48..57 are all DISTINCT, so on every real crossword alphabet the
  mask is *exact* (no collisions). The test is conservative by construction: a
  shared letter is the same char → same bit → survives the AND, so AND=0 PROVES
  disjoint (no false negatives → `mrv_inc`'s no-under-count invariant is safe).
  A hypothetical mod-64 collision yields only a false POSITIVE (fall through to
  the exact path / an unnecessary-but-exact recount) — safe. Masks stored as a
  **tabled `answer_mask/2`** keyed on the answer atom (a pure function, scanned
  once per answer per process); every placed dict already carries `answer`, so no
  dict-shape change and no golden-leak risk. The `inc_count` change also drops the
  per-node `entry_letters/2` (atom_chars) call, replaced by the tabled lookup.
- **Result (goldens byte-identical, 66 plunit + goldens green):** a REGRESSION
  everywhere. Ladder (`make bench-check --heavy`) rose on 8 of 9 rungs, +1.3% to
  +2.9% inferences (9x9/08w +2.31%, 15x15/28w +1.91%, 15x15/32w +2.06%, 9x9/16w
  +2.56%, 21x21/80w +1.34%, 15x15/34w +2.87%, 15x15/36w +2.77%; only 15x15/12w
  ~flat at -0.48%). `bundled_17` (real 26-letter words, in-process strict arrange)
  also regressed, 673,333 → 677,707 inferences (+0.65%).
- **Why (measured pre-filter skip rates):** the prefilter almost never fires, and
  what it gates is already cheap. Fraction of pair tests with AND=0:
  `find_intersecting_word` 2.6% (`bundled_17`), 6.1% (9x9), 11.6% (15x15/32w),
  13.7% (21x21/80w); `inc_count` 10-17%. Crucially the REAL-word fixture had the
  LOWEST find-intersecting skip rate (2.6%): crossword answers (esp. multi-word,
  e.g. `NARRATIVE FALLACY`) have many distinct letters → dense masks → pairs
  almost always share a bit. So "26-letter alphabet ⇒ more disjoint pairs" is
  wrong — per-WORD letter density, not alphabet size, decides, and it is high.
  Every call pays the mask machinery (candidate mask rebuild + a tabled
  `answer_mask` lookup + the AND) on the hottest path; the saved work is one
  `intersection/3` over short lists, only 3-14% of the time. Constant cost >
  constant saving → the net rise.
- **Verdict:** REJECTED, do not revisit on this solver. Unlike the tabled
  `pair_crossings` design the hypothesis assumed (skip an expensive hash+copy
  table lookup), this codebase's `find_intersecting_word/6` does a cheap inline
  `intersection/3`; there is no expensive lookup to amortize away, so an O(1)
  prefilter with its own non-trivial constant cannot pay off at these skip rates.
  (A bitmask would only help a variant that made the *gated* operation expensive.)

### E-H10 — watched witnesses for the incremental recount — REJECTED (branch only)

- **Idea (watched-literals analogue):** cache, per word, the WITNESSES for its
  last count - the first up-to-2 legal placements, each a re-checkable
  `w(Start,Dir)`. At a recount trigger (a word shares a letter with the just-placed
  word), first re-validate the cached witnesses with `check_word_fits/5`
  (O(word length) each); if they still hold, KEEP the previous bucket and skip the
  full `find_intersecting_word x pair_crossings x check_word_fits` rescan.
- **Where (measured, branch `experiment/e-h10-watched-witnesses`, not merged):**
  core.pl only - `count_upto2_wit/5` (count_upto2 that snapshots ground
  `w(Start,Dir)` copies into the mutable holder via nb_setarg, order- and
  count-preserving), `mrv_count_wit/8`, the cache value widened `Count -> wc(Count,
  Witnesses)` threaded in the existing `state(CountMap, LastLetters)` (pure
  backtrackable arg, so it unwinds on backtrack - the one hard requirement), and
  `recount_or_reuse/9` on the inc_count shares-letter branch. `count_of/3` unwraps
  the wc on the hot ordering loop.
- **Soundness — the headline mechanism is UNSOUND (verified in code, not on
  faith).** The brief's monotonicity premise ("within one descent a word's legal
  set only SHRINKS") is FALSE for exactly the words this optimisation touches. The
  recount trigger is "shares a letter with the last placement", which is precisely
  the case where placing that word can ADD a new crossing placement (the mrv_inc
  invariant, core.pl §"Correctness rests on this invariant"): a shares-letter
  word's count can RISE, not only fall. So a bucket-1 word whose single witness
  still validates can nonetheless have TRUE count 2 (a fresh placement crossing the
  just-placed word), and keeping bucket 1 UNDER-counts. Empirically confirmed: the
  brief's bucket>=1 policy left `all_crossword` full-tree solution counts identical
  (the solution SET is unchanged - completeness holds) but changed the emitted
  first layout on 11 of 12 ladder rungs (only 21x21_25w survived) - a byte-identity
  violation. Only the SATURATED bucket is sound: `Count == 2` (== cap) cannot rise
  past the cap, so both-witnesses-valid => still >= 2 => same bucket. The shipped
  experiment therefore short-circuits ONLY `Count == 2`; buckets 0/1 always recount.
- **The sound scope excludes the volume.** P1 found count<=1 is the plurality
  bucket at nodes, and a full recount of a bucket-2 word already EARLY-EXITS at the
  2nd hit (count_upto2). So the sound (bucket-2-only) short-circuit skips
  candidate-regeneration on the recounts that were already the CHEAPEST, while a
  flat witness-capture + `wc` unwrap tax is paid on every count/every ordering pass.
- **Witness-survival instrumentation (first-solution descent, topleft_across;
  temporary counters, removed before commit):** of shares-letter recounts, the
  fraction that are bucket-2 and the fraction of those whose witnesses survive:
  21x21_80w 53.8% bucket-2 x 65.7% survival (fires on 35.3% of triggers);
  21x21_82w 36.2% x 45.3% (16.4%); 15x15_36w 30.5% x 45.7% (14.0%);
  09x09_08w 29.7% x 18.2% (5.4%); 15x15_12w 29.8% x 21.4% (6.4%). The hit rate
  tracks fan-out: high on the dense 21x21 rungs, low on the light/trail rungs.
- **Result (`make bench-check --heavy`, warm, vs baseline.json; sound bucket-2
  version, all 12 rungs byte-identical):**

  | rung | baseline | measured | delta | status |
  | --- | --- | --- | --- | --- |
  | 09x09_08w | 28,212 | 28,616 | **+1.43%** | **REGRESSION** |
  | 15x15_12w | 101,588 | 102,808 | **+1.20%** | **REGRESSION** |
  | 21x21_25w | 302,257 | 304,600 | **+0.78%** | **REGRESSION** |
  | 15x15_28w | 374,521 | 369,794 | -1.26% | win |
  | 15x15_32w | 969,070 | 936,588 | -3.35% | win |
  | 09x09_16w | 674,961 | 656,651 | -2.71% | win |
  | 21x21_80w | 4,487,855 | 3,658,866 | **-18.47%** | win |
  | 15x15_34w | 13,650,369 | 13,556,053 | -0.69% | win |
  | 15x15_36w | 38,275,505 | 38,375,768 | +0.26% | ok |
  | 09x09_17w | 38,533,195 | 38,600,393 | +0.17% | ok |
  | 15x15_40w | 10,411,252 | 10,359,444 | -0.50% | ok |
  | 21x21_82w | 6,822,207 | 6,410,131 | -6.04% | win |

  Ratchet verdict: **FAIL** (3 regressions, 6 wins). The gradient is exactly P1's
  prediction - the counting-dominated, high-fan-out 21x21 rungs gain most
  (80w -18.5%, 82w -6.0%), the trail/low-fan-out rungs least (they REGRESS on the
  flat tax).
- **Identity:** 204 plunit pass; all goldens byte-identical; `all_crossword`
  full-tree solution counts identical across baseline/mrv_capped/mrv_inc x 4 corners
  on 2 fixtures (09x09_08w, bundled_17); arrange output byte-identical on all 12
  ladder rungs. (The sound version only; the bucket>=1 version fails identity - see
  Soundness above.)
- **Verdict:** REJECTED - branch only, nothing merged, not bench-recorded. The
  headline win the hypothesis targeted (making the PLURALITY count<=1 recounts
  O(word length)) is unreachable: that bucket is exactly where a shares-letter
  recount can gain a placement, so the witness check cannot certify it without a
  rescan. The sound residue (bucket-2) helps only the high-fan-out counting rungs
  and taxes the rest past the 0.5% gate. This closes the "cheaper per-node counting
  via cached placements" idea in its stated form. The ONE avenue that could still
  ship the 21x21 win: keep witnesses in a SPARSE side-assoc (only ever-bucket-2
  words) so `count_of` stays a bare-int read and the flat ordering-loop tax
  disappears - untried here; it is a different experiment, not this hypothesis.

### E-H10b — sparse-side-map watched witnesses — REJECTED (light-path structural tax removed; the shortcut still cannot pay its own admission test on the smallest rungs)

- **Idea:** rebuild E-H10's sound Count==2 witness shortcut so the light path
  is untouched. Keep CountMap values BARE INTS (no ordering-loop unwrap); hold
  witnesses in a SPARSE side-map keyed only by cap-saturated words; confine
  witness capture to recounts. Aim: keep the 21x21 win with ZERO regressions.
  Branch experiment/e-h10b-sparse-witnesses, not merged.
- **Where:** core.pl only. New: mrv_count_wit/8 + count_upto2_wit/5 +
  wcounter_witnesses/2 (the witness-capturing count, reused near-verbatim from
  92179c8 but reachable ONLY from the recount fallback); recount_or_reuse/9 +
  witness_valid/4 (the shortcut). Changed: inc_count/11's shares-letter branch
  gains a `get_assoc(A, PrevMap, 2)` gate that routes ONLY previously-saturated
  words into the shortcut; inc_counts(none,...) resets the witness map. Every
  other line of the mrv_inc counting path (count_of, select_inc, the state/2
  term, full_count, the carry-forward branch, inc_counts(state,...)) is
  BYTE-FOR-BYTE master.
- **Design decisions (deliverable b):**
  - *Side-map location.* Threading the map through the state term via the foldl
    accumulator (a cw(CountAcc,WitAcc) pair) was tried FIRST and REGRESSED the
    light rungs +1.1..+3.1% - the per-word compound wrap/unwrap taxes every
    carry-forward word. Instead the map lives in a BACKTRACKABLE GLOBAL
    (b_setval/b_getval `cw_witnesses`): itself backtrack-restored state, but the
    carry-forward path and the ordering loop never touch it, so they stay
    byte-identical. Reset at inc_counts(none,...) (the per-subtree entry both
    the search and greedy/fragment drivers pass through); sibling seeds start
    clean via ordinary backtracking.
  - *Capture confinement.* Witnesses are captured only in the recount FALLBACK
    (mrv_count_wit), never in the initial full pass and never on carry-forward.
    The `PrevMap==2` gate further restricts capture to words that were saturated
    last node - the alternative (gate on sparse-map membership, capturing on
    every shares-letter recount) REGRESSED +2.4..+3.2% and is recorded as the
    losing variant.
  - *Stale-entry policy.* A recount OVERWRITES the entry on re-saturation and
    LEAVES a stale entry when the count drops (the zero-assoc-op branch). SOUND
    REGARDLESS of staleness: two currently-valid distinct witnesses are two real
    placements, so validation succeeding proves true count >= 2; if the word
    truly dropped, at most one witness validates and we fall through. A blocked
    witness cell stays blocked deeper on a descent, so a dead entry never
    mis-fires - it wastes at most two check_word_fits.
- **Result (--heavy, all 12 rungs, inferences vs baseline.json):** the 21x21
  win is FULLY retained - 21x21_80w 4,487,855 -> 3,698,716 (-17.6%), 21x21_82w
  6,822,207 -> 6,440,060 (-5.6%); also 15x15_32w -2.6%, 09x09_16w -2.7%,
  15x15_28w -0.8%. But THREE light rungs still regress past the 0.5% ratchet:
  09x09_08w +1.5%, 15x15_12w +1.4%, 21x21_25w +0.6%. The E-H10 structural taxes
  are GONE (36w +1.07% -> +0.26%; 17w +1.77% -> +0.29% - both now pass), but the
  smallest rungs do not clear the bar.
- **The decisive finding (why this closes the avenue):** a diagnostic that keeps
  ONLY the gate probe (`get_assoc(A, PrevMap, 2)` on every shares-letter word)
  with the shortcut and capture DISABLED still costs +0.60..+0.74% on the four
  smallest rungs (08w +0.74%, 12w +0.71%, 25w +0.60%, 28w +0.60%) and +0.2..+0.4%
  everywhere else. That is the IRREDUCIBLE admission test the sound shortcut
  needs: to decide a word is shortcut-eligible you must look up whether it was
  cap-saturated, and that single bare-int assoc probe per shares-letter word is
  ALREADY over the ratchet on the smallest rungs, BEFORE the shortcut does any
  useful work. On the high-fan-out rungs the -18%/-6% savings bury it; on the
  light rungs (hit rate ~5-6%, and a bucket-2 rescan already early-exits at the
  2nd hit) there is nothing to bury it with. No size-independent signal skips the
  probe only where it cannot pay off, and an input-size switch / magic constant
  is off-limits by charter.
- **Identity verification:** 204 plunit green; goldens BYTE-IDENTICAL (not
  regenerated); arrange output BYTE-IDENTICAL on all 12 ladder rungs vs master;
  FULL-TREE all_crossword counts identical across all four strategies
  (baseline/mrv/mrv_capped/mrv_inc) x 4 corners on 2 non-zero fixtures
  (bundled_17: 16315/16315/4238/4238; ladder_09x09_08w: 36553/36553/39040/39040).
  The shortcut is SOUND - it only changes speed, never the tree.
- **Verdict:** REJECTED. E-H10b did its narrow job - it removed E-H10's
  bookkeeping tax (bare-int count map, byte-identical light path, capture off the
  light path) and kept the entire 21x21 win - and by doing so ISOLATED the real
  obstacle: not the bookkeeping, but the per-shares-word eligibility lookup the
  sound Count==2 scope inherently requires. That lookup is a flat ~0.6-0.7% tax
  on the smallest rungs that the shortcut's light-rung savings cannot offset.
  **This closes the watched-witness avenue** (E-H10 + E-H10b). The E-H9 lesson
  restated: a per-node cost that is minimal on the heavy rungs can still be the
  whole budget on the light ones - and here it is the shortcut's own admission
  test, which no amount of bookkeeping cleverness can remove.

### A2-KS — Fenwick order-statistic replay for seeded permutation — KEPT

- **Idea:** replace `seeded_permutation/2`'s O(n²) survivor-list selection
  with an invocation-local Fenwick tree over live original positions. Keep
  every PRNG draw, `V mod remaining-count`, selected live rank, and deletion
  exactly unchanged.
- **Where:** `prolog/crosswordsmith/core.pl`; fill's default pinned
  `prng_shuffle/2` is untouched. Built on `experiment/a2-seeded-fenwick`,
  commit `42abb6a`.
- **Soundness:** an independent test-only implementation of the old walk
  matched the Fenwick implementation for seeds 0, 1, 7, 42, and 2^64-1 at
  every length 0..256. All 1,285 permutations and the following PRNG draw
  matched. The published seed-42 answer and empty-list no-draw, determinism,
  failure, and error contracts remain locked.
- **Result:** the trivial 3×3 `cwl50` seeded CLI reproducer fell from 171.14s
  to 9.57s, a **17.88×** speedup, with exit 0 in both cases. Peak RSS moved
  from 785,912 to 819,404 KiB (+4.26%). A 172,823-word plain seeded load/fill
  completed in 4.91s at 601,148 KiB.
- **Identity and ratchets:** `make test` passed with every golden unchanged;
  all 11 fill identity digests matched. Heavy fill search and load inference
  counts were exactly unchanged on all 11 rungs (+0.00% throughout). Fuzz
  passed 66/66 and `make test-wasm` confirmed native/LibBF seeded parity. No
  baseline, identity, or golden data was regenerated.
- **Verdict:** **KEEP.** The historical seeded sequence is replayable in
  O(n log n); A2-KS clears the <=15s and >=10× gates with zero output or
  inference change. Residual cost is about 33 MiB additional peak RSS on the
  measured `cwl50` run. The re-contract arm is unnecessary and was not taken.

### B0-I — shipped MAC instrumentation — MEASURED, all Track-B arms survive

- **Method:** benchmark-only exact-replay twin at commits `42c7eac` and
  budget-boundary correction `9a88f36`. Setup and unchanged helpers call
  `crosswordsmith_fill` internals; local twins only instrument restart,
  search, placement, propagation, revision, and bump. The correction moved
  MAC setup inside the same inference limit as root/search, matching product
  `fill_attempt/8`; the first CWL stop numbers were discarded and rerun.
- **Soundness:** product/twin outcome, exact layout/words, and final-attempt
  nodes matched on four easy/seeded/scored rows, both multi-attempt authority
  rows (`919/919` @30, `228/228` @1), and all five quality masks. Duplicate
  authority/CWL counter runs matched every non-wall field. Counter overhead
  was +5.39%, below the pre-registered 15% ceiling.
- **Reference result:** STW `blocked_13a @30/@1` filled in 7/2 attempts and
  11,318/729 nodes. Queue length <=2 was only 0.22%/1.04%, while length >=4
  was 99.51%/98.08%. Fruitful revisions were 27.89%/28.93%; all DWOs came
  from propagation, not placement.
- **Learning signal:** 10,945/599 edge bumps; the top edge quartile held
  83.54%/91.49% of bumps and Gini was 0.737/0.799. F1 queue ordering, F2/F4
  credit variants, and F3/F5 aging/probing therefore all **SURVIVE** their
  numeric kill-tests. This is structural headroom, not a win verdict.
- **Cost mechanism:** propagation/support consumed 99.54%/99.49% of search
  wall. Track-B latency arms must win by reducing bignum propagation work;
  at `@1`, setup still accounts for about half of end-to-end wall.
- **Matrix:** all 11 ladder rows filled; five-mask assignments matched the
  product with min/mean 50; corrected `cwl50 @50` stopped `not_proven` after
  charging MAC setup and search to the shipped 800M budget (27,342 nodes,
  236.1s end to end, before the 240s cap). Product tests, all identity rungs,
  and both inference ratchets passed at +0.00%; no artifact was re-recorded.
- **Verdict:** **B0 complete; build B1, then B2/B3 in the registered order.**
  Runnable evidence:
  `benchmarks/results/2026-07-16-fill-b0-mac-instrumentation.md`.

### B1-O — MAC revision ordering — MEASURED, no arm graduates

- **Method:** rejected probe branch `experiment/b1-revision-order`, commit
  `76147fc`. Its baseline call graph reproduces corrected B0 authority/CWL
  counters exactly; seam checks lock stable ties, full-unfilled-neighbour
  weighted degree, inactive-edge retention, and shared-cell bind validity.
  Eight authority duplicate pairs matched fills, limited-goal inferences,
  counters, attempts, and weights exactly.
- **`q_dwd`:** authority inferences -70.4%/-49.1% and revisions
  -86.2%/-80.0%, but `@30` mean fell 45.0 → 43.9. The ladder also exposed
  `g17_50k` +4,381% inferences/+6,969% nodes and `g09_full` +31%/+352%.
  **REJECT.** Conditional `q_dwd_edge` was therefore not built.
- **`q_wdeg`:** authority inferences -66.7%/-27.3%, revisions
  -82.4%/-62.8%, means 46.5/40.0. It still failed the registered ladder
  guard: `g17_full` +6.0% inferences and `g17_50k` +184.8% inferences/+298.0%
  nodes. **REJECT as a standalone/global policy.**
- **`edge_desc`:** exhausted the `@30` 800M budget where baseline fills and
  regressed `@1` +15.6% inferences/+7.6% revisions; `g17_50k` was +572.5%.
  **REJECT.** High learned weight is not a useful local ASAP edge signal here.
- **Matrix:** every ladder variant still filled, all five quality masks stayed
  mean/min 50, and no arm completed `cwl50 @50`. Product tests, identity, and
  both inference ratchets remained unchanged because all code was probe-only.
- **Mechanism/verdict:** B0's premise was right: ordering can remove most
  expensive wide-mask revisions. The failure is generality: order changes
  conflict learning and later restart trajectories. No B1 arm enters Track D.
  One distinct B2 sequel is warranted by the thesis's interaction result:
  pair `q_wdeg` with `alldel`/fully-assigned order-robust credit and require it
  to remove the two-rung regression before it can graduate.

### B2-C — MAC conflict credit — MEASURED, no arm graduates

- **Method/soundness:** rejected probe branch `experiment/b2-credit`, commit
  `80e3edf`. Executable seams pin target-specific H1/H2/H3 formulas,
  immediate/final `alldel`, whole-episode distinct fully-assigned credit,
  fresh backtrack-safe episode lifetime, persistent learned weights, and
  B1-exact q_wdeg ties/denominator. Baseline reproduces B0 exactly; all 16
  authority duplicate pairs match every non-wall field.
- **Completion:** H1, H3, `alldel`, fully-assigned, `alldel_qw`, and
  `fully_qw` all lose STW `@30` completion at 800M. Most also produce severe
  `@1` and ladder regressions.
- **H2:** the only two-row completer; `@30` inferences/revisions
  -43.4%/-49.5%, but quality mean 44.4 < 45.0, `@1` only -1.2%/-12.6%, and
  `g17_50k` becomes not-proven at 2B. **REJECT.** Deleted-count magnitude
  creates highly concentrated, search-shape-specific guidance.
- **Broad credit mechanism:** `alldel`/fully-assigned flatten learned-weight
  Gini from 0.737/0.799 to 0.257-0.328 standalone. More deletion evidence is
  not better conflict guidance; dom/wdeg loses discrimination and restart
  trajectories deteriorate. The q_wdeg combinations do not stabilize B1.
- **Matrix/verdict:** all five quality masks stay mean/min 50; no arm
  completes CWL. No B2 arm enters Track D, and alternate-credit/order
  combinations are closed. B3 remains independently warranted using baseline
  credit only. Product tests, 11 identity rungs, and heavy ratchets pass
  unchanged because all variants remain probe-only.

### B3-R — MAC aging and probing init — MEASURED, no arm graduates

- **Method/soundness:** rejected probe branch `experiment/b3-aging-probing`,
  commit `9782772`. Seams pin normal-failure-only age counting, cap/exception
  bypass, exact x0.5 aging, per-attempt reset, same-root cap-125 probes,
  discarded probe fills, and a separate pinned splitmix64 stream that cannot
  shift production draws. All eight authority duplicate pairs match exactly.
- **`age20`:** `@1` inferences/revisions -30.8%/-38.0% and no ladder rung
  over +5%, but `@30` becomes not-proven at 800M after 871 halvings.
  Persistent concentrated weights are necessary on the wide row. **REJECT.**
- **`probe4`/`probe8`:** genuine `@30` guidance (-50.9%/-90.2% inferences),
  but `@1` regresses +8.6%/+38.6%, `@30` mean is 44.6 < 45.0, and both lose
  `g17_50k` at 2B. Warm cost fails directly: 504/940 warm nodes save only
  156/490 real `@1` nodes. **REJECT.** Conditional combination not built.
- **Closure:** all five quality masks remain mean/min 50 and no arm completes
  CWL. A smaller probe dose cannot satisfy the global ladder node guard because
  the easy baseline rungs take only 8–107 nodes before any warm tax. Track B
  closes with no Track D candidate. C1 remains the independent envelope bet.

### C1-P — exact projection set branching — LOST

- **Method/soundness:** rejected probe branch `experiment/c1-projection`,
  commit `896ea1c`. Class masks are lazy checked-letter bignum intersections;
  complete score/bit-order leaf matching carries seed/global all-different and
  backtracks projection search; shared-cell binding validates every fill.
  Five executable seams cover propagation equivalence, used-representative
  alternate, Hall backtracking, exact high-bit coverage, and class top-3.
- **Target result:** **0/5** unclosed rows complete. Across target states,
  77,495 classes represent 77,513 words (1.00023×); only `blocked_13b @1`
  sees non-singletons, just 18 size-2 classes among 22,388. The assumption
  that unchecked cells create large exact projection classes is falsified.
- **Ladder mechanism:** every class is singleton, so sound deferred
  all-different buys no quotienting and postpones duplicate rejection to
  complete assignments. Of 526,651 matching calls, 526,643 fail;
  `g17_full`, `g21_full`, and `g17_50k` lose completion at 2B.
- **Gates/verdict:** authority completion and quality pass; five-mask min/mean
  remains 50. No envelope completion plus three ladder losses makes C1
  **LOST**. Approximate clustering would invalidate whole-class refutation
  and was not attempted. DP-6 pins stand; no Track D dossier input.

### Arrange P0 — strict/greedy ratchets and probe substrate — INFRASTRUCTURE, PASS
- **Scope:** composed P0.1-P0.4 on `campaign/arrange-p0-integration` from
  `aa4d3c8`, with strict identity/promotion, a separate four-layer greedy
  benchmark, counter-free product authority, exact-replay instrumentation,
  fixed cliff fixtures, and the attempt JSONL schema. Product search code and
  existing goldens are unchanged.
- **Strict gate:** core and heavy inference ratchets reproduced exactly. All 15
  fresh-process identity rows passed twice; the 500M latency incumbent retained
  reward 768 and digest
  `90289af7db529b0132bae8bb910a18e90daa6f7b1c19de7cb97508883a56c81b`.
  Promotion read-back tests cover retained heavy rows, new rungs, and
  latency-only rows; integration also makes baseline replacement atomic.
- **Greedy gate:** two independent heavy runs reproduced all 21 construction,
  sweep, and postprocess inference cells at +0.00%. The ordered raw-pool and
  selected-output identity passed twice with digest
  `ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29`.
- **Probe gate:** authority matches direct product at tiny and 500M budgets;
  easy, hard, seeded, and second-representative twins match outcome, reward,
  layout signature, and lean/full decisions. On the serialized adjudication
  rerun, lean wall overhead was +1.36% and -4.11% (inferences +3.65% and
  +1.38%), below the 15% steering ceiling. Full wall overhead was +124.2% and
  +104.2% (inferences +132.61% and +22.50%), so full remains sampled-only.
- **Frozen inputs/schema:** all 12 requested-count cliff fixtures and all 64
  SplitMix64 seeds regenerated identically twice. Python and product seed
  streams agree; all 10 schema tests pass, including mixed-rig rejection and
  configured-cutoff/termination separation.
- **Verification/verdict:** composed `make test` passed 411 plunit assertions,
  every golden comparison, and all CLI contracts. Strict documentation says
  two non-transpose representatives under one shared operation budget; greedy
  remains four-corner. **Phase 0 PASS; A-G1 may begin.** Detailed probe evidence:
  `benchmarks/results/2026-07-16-p0-arrange-probes.md`.

### A-G1 premise — legality before score — MEASURED, PRODUCT EXPERIMENT GATED IN
- **Method/soundness:** benchmark-only exact replay at base `3ce0b8d`, probe
  commit `ea65c1a`. Product and replay matched setup outcomes, layouts, rewards,
  and dropped order for every corner x seed attempt on all seven greedy rungs
  before and after counter instrumentation.
- **Result:** legality rejected 51,802/55,080 scored candidates (94.05%) on
  15x15/32w and 273,312/280,386 (97.48%) on 21x21/80w. Legality-first would
  avoid exactly 388,884/412,496 (94.28%) and 2,000,112/2,050,244 (97.55%)
  current score-cell visits. Both targets clear the pre-registered 25% gate.
- **Correction:** `find_intersecting_word/6` produced 6,621 and 25,601
  `Start < 1` descriptors on the dense targets, but `check_word_fits/5` guards
  `Start >= 1` before any `arg/3`; the product comment claiming score-first is
  needed for safety is stale. Legality-first preserves surviving candidate
  order and first-tie-wins because both operations are pure.
- **Verification/verdict:** full greedy identity remained
  `ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29` and
  focused tests passed 7/7. **Gate A-G1 into a separate product experiment;**
  no product, baseline, identity, ratchet, or golden changed.

### A-G1 — legality before score — KEPT
- **Change/soundness:** integration commit `ef2e0b0` moves pure
  `check_word_fits/5` before `placement_key/8` in the greedy candidate
  generator. Legality guards `Start >= 1` before any grid read and binds no
  cell; scoring is likewise pure. All seven exact replays preserve setup
  outcomes, placed terms, rewards, dropped order, raw pools, selected layouts,
  first-tie-wins, and CLI bytes.
- **Inference result:** construction/sweep changed by -25.81%/-23.18%
  (bundled17), -1.04%/-0.47% (bundled11), -17.07%/-17.76% (benchmark08),
  -47.15%/-48.53% (real13), -54.99%/-59.26% (real15), -60.06%/-61.88%
  (15x15/32w), and -51.87%/-67.65% (21x21/80w). Postprocess was exactly
  unchanged on every rung: 13 gated wins and zero regressions.
- **Mechanism/wall:** scored candidates fell 55,080 -> 3,278 and
  280,386 -> 7,074 on the dense targets, eliminating 388,884 and 2,000,112
  score-cell visits. Two warmup plus 21 serialized baseline/candidate pairs
  gave sweep-wall median ratios 0.401 (p10-p90 0.380-0.417) and 0.344
  (0.338-0.353), or 59.95% and 65.59% median improvements.
- **Verification/adoption:** composed `make test` passed 412 assertions and all
  goldens/CLI contracts. Strict identity and all strict core+heavy ratchets
  stayed exact; greedy identity remained
  `ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29`.
  `bench-greedy-promote --heavy` persisted and read back all seven improved
  rows and appended history. **KEEP; A-G2 now starts from this ratcheted base.**

### A-G2 premise — greedy direct transpose pairs — MEASURED, PRODUCT EXPERIMENT GATED IN
- **Method/soundness:** measurement-only benchmark replay at base `c584962`,
  probe commits `7d82081`/`d7b5dfb`. It enumerated the unchanged four
  corner-major then seed-major blocks, retaining setup failures, strict
  eligibility, complete `pw/8` records, and original dropped terms. Direct
  pairs were checked with an independent row/column transpose, involution,
  normalized direction, raw-pool order, and fresh clue variables.
- **Result:** all 108 manifest slots (54 pairs) matched across all seven greedy
  rungs with zero setup, geometry, placed-order, dropped-order, reward,
  eligibility, normalized-assoc, or raw-pool mismatches. Attempt counts were
  `20,20,20,20,16,8,4`; 12 setup failures paired symmetrically. Three
  additional generated square sets checked 28 slots spanning setup failure,
  dropped words, and all-fit behavior, also with zero mismatches.
- **Implication/verdict:** direct attempts can fall exactly 108 -> 54 (50%) if
  the product searches `topleft_across`/`topright` only and inserts freshly
  rebuilt transpose partners into the omitted historical slots. Dropped terms
  must retain original terms/order and each rebuilt `pw/8` must receive a fresh
  non-aliased clue variable. **Gate A-G2 product experiment in.**

### A-G2 — derive greedy transpose partners — KEPT
- **Change/soundness:** integration commit `7a14082` directly searches the TLA
  and TR seed blocks under the shared memo reset, rebuilds literal square-grid
  TLD/BL transposes with fresh clue variables and copied original dropped
  terms, and concatenates the four blocks in historical TLA/TLD/TR/BL order.
  Failed source setup omits both symmetric entries; strict eligibility remains
  after construction.
- **Identity:** the independent four-corner observer matched all 108 manifest
  slots and 28 generated slots against synthesized pools with zero mismatches.
  Ordered raw pools, ties, dropped order, rewards, selected layouts, normalized
  assocs, candidate counts/distances, and CLI bytes are unchanged; greedy
  identity remains `ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29`.
- **Result:** direct attempts halve exactly by rung:
  `20->10,20->10,20->10,20->10,16->8,8->4,4->2`. Sweep inferences improve
  44.05%, 35.45%, 44.06%, 44.15%, 44.82%, 46.83%, and 47.55%, with zero
  construction/postprocess or strict regressions. Dense serialized wall ratios
  are 0.5436 and 0.5500 (paired p10-p90 0.5308-0.5583 and 0.5355-0.5666),
  approximately 1.84x and 1.82x.
- **Verification/adoption:** composed `make test` passed 418 assertions and all
  goldens/CLI contracts; direct-vs-derived and greedy identity pass. Branch
  verification covered all 15 strict identities and every strict core+heavy
  inference gate exactly. `bench-greedy-promote --heavy` persisted/read back
  all seven wins and appended history. **KEEP.**

### P-D0 — support delta and proof multiplicity — MEASURED, A-D1/A-D2 NOMINATED
- **Rig/soundness:** benchmark-only exact replay at shared Phase 2 base
  `1bccf47`, integrated as `9ee66ed`. Pure threaded residues and independent
  proof scans cannot influence MRV, grids, tables, or budget semantics. Previous
  buckets 0/1 are classified from newest-source proofs plus the sole old
  residue; bucket 2 always falls back. Every claim was checked against a full
  exact Cap=2 recount.
- **Result:** classified 186/246 light refreshes (75.61%) and 2,985/6,232 dense
  refreshes (47.90%). Hypothetical candidate checks changed `2,191 -> 954`
  light and `144,647 -> 67,572` dense, with reductions on every light corner.
  All 3,171 classifications were exact.
- **Multiplicity/representation:** 31 recounts had more legal proofs than
  unique geometries (`153/151` light, `8,122/8,093` dense), killing geometry
  dedup. Dense controls also performed 7,544 assoc writes, 8,636 assoc reads,
  and sorted 5,798 items over 322 interior nodes, warranting isolated A-D1.
- **Verdict:** nominate A-D1 with zero light-rung regressions mandatory; if it
  passes, run proof-preserving A-D2 for previous buckets 0/1 only. Do not reopen
  geometry-only or bucket-2 residue designs without a complete proof model.

### P-C0 — duplicate failed work — MEASURED, GLOBAL PREMISE PASSES; A-C1 KILLED
- **Rig/soundness:** observe-only exact replay at `1bccf47`, integrated as
  `55a0fde`. Absolute canonical keys contain grid size, sorted remaining stable
  IDs, and sorted `(ID,start,dir)` placements. Dual fingerprints only select a
  bucket; every hit is full-key verified. States become dead only after normal
  exhaustive failure, never cutoff or interruption.
- **Result:** the 15x15/34w and 15x15/36w hard `topright` controls revisited
  proved-dead states at 2,763/3,107 (88.93%) and 3,019/3,454 (87.41%) recursive
  entries. Each repeated subtree was a one-node wipeout, so node share is not
  yet an inference/wall claim.
- **Parent-local result:** same-parent legal duplicate children and post-failure
  captures were zero on both hard controls and all sentinels despite
  4,006/8,347 duplicate crossing proofs. A-C1 misses its 80% gate at 0% and is
  closed.
- **Verdict:** admit one separately gated global negative-state-cache candidate,
  registered as A-C2, subject to zero light regressions, <=32 MiB RSS, exact
  identity, and WASM checks. Do not build parent-local A-C1.

### P-R0 — fixed-instance seeded distributions — MEASURED, TRACK R GATED IN
- **Protocol:** unchanged counter-free standalone authority at both
  non-transpose corners for four fixed cliff instances over frozen pilot seed
  indices `[0,16)`, plus one completing control per grid size. Integrated runner,
  131-row JSONL, and report are in `c73a049`; the schema CLI fix is separate at
  `910ddd2`.
- **Result:** cliff rows produced 60/128 placements and 68/128 right-censored
  `not_proven`, with no infeasible or interrupted rows. Either-corner completion
  ranged from 10/16 to 15/16 and every fixture gained +1 to +4 seeds over its
  better individual corner. The three controls completed.
- **Budget/scope:** consumed 1,887.630 child CPU-seconds (0.524 hours), below the
  1.5-hour pilot cap. The early stop did not apply. Tuning `[16,32)` and held-out
  `[32,64)` remain untouched; standalone complementarity is not a live
  operation-wide controller claim.
- **Verdict:** Track R passes its premise gate. A-R1 remains parked at the
  explicit output-changing product-policy checkpoint before any tuning or
  tournament is run.
