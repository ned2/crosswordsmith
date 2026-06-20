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

## Result batch — 2026-06-20, post I5 fix (SWI 10.0.2) — current

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

- All logged algorithm ideas (I1-I5) are closed. Possible future work: a
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
