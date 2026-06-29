# `arrange` engine — implementation plan

*Build spec for Flavour A's unified layout engine. The **design rationale** for every decision encoded here lives in [`cryptic-setter-research.md`](./cryptic-setter-research.md) → "Addendum — Two flavours, one CLI: the Flavour-A engine spec"; this document is the forward-looking implementation plan.*

## Goal

A **deterministic MRV-first layout engine with admissible-bound pruning and within-budget improvement** that places words into an N×N grid maximizing a capped interlock objective, with `--strict`/`--best-effort`, `--size-mode fixed|max`, fragment seeding, and diverse `--candidates N` — replacing the path-dependent `pack` (`quality.pl`) / `solve` (`crossword.pl`) split with one `arrange` verb. *(Earlier drafts called this "a branch-and-bound optimizer"; the design-exploration pass showed that label oversells what ships — see "Empirical reframe" below.)*

**Key framing:** `arrange` is a new *search driver + objective* on top of placement machinery that **already exists and is correct**. The hard part of crossword placement — legal intersection, no-merge, adjacency, the grid model — is all reusable. What is new is *what the search optimizes* and *when it stops*.

**Empirical reframe (design-exploration pass, 2026-06-30).** Probes on the real fixtures established: (a) the full solution space is explosive — 6 short words → **41,106** layouts, and 8/14-word enumeration does not finish in 30 s — so the admissible bound is expected to **prune little** at our sizes; and (b) the per-word cap is **often inert** on realistic sparse inputs (toc_demo: **1/16** words reach `ceil(L/2)`; `--min-half` makes it infeasible), so the objective **degenerates toward plain total-crossings** there (it stays active on dense meshes — quality_22: 19/22). Net: output quality is governed primarily by the branch-ordering / construction heuristic plus from-scratch re-scoring, with the bound an *opportunistic* extra — not by the search proving an optimum. **A Phase 1.5 gate measures whether the search layer earns its place at all before Phase 2 builds it.** A second structural consequence (below): best-effort and candidates ride the **greedy constructor**, not a drop-branch on the strict DFS, because the DFS cannot even seed the dense meshes where dropping matters.

## Reuse map

| Capability | Status | Predicate(s) |
|---|---|---|
| Square grid model (assoc `cell→empty\|letter`, row-major numbering) | reuse as-is | `init_grid/2`, `cell_coord/3`, `next_cell`/`prev_cell`/`calc_num`/`calc_start`/`fits_on_grid`/`is_start_cell`/`is_end_cell` |
| Placement legality (the hard part) | reuse as-is | `assign_word/10`, `assign_letters/7`, `adj_is_free/4`, `no_word_merge/3`, `check_prev_cell/4`, `check_next_cell/4` |
| Candidate placements crossing existing words | reuse as-is | `find_intersecting_word/6` |
| MRV branch ordering (strong-incumbent heuristic) | reuse | `mrv_count/8`, `select_inc/10`, the `state(CountAssoc, LastLetters)` cache |
| Crossing count of a placement | reuse (from `quality.pl`) | `crossing_count/6` |
| Bounding box (for `max` cropping) | reuse (from `quality.pl`) | `placed_bbox/4` |
| Half-check threshold = objective cap | reuse (from `quality.pl`) | `word_meets_half/2` → `ceil(L/2)` |
| **Greedy constructor (best-effort + candidates path)** | **reuse (from `quality.pl`), re-keyed** | `greedy_construct`/`greedy_loop`/`next_move`/`word_best_placement` — swap `placement_key` for the capped-objective delta |
| Clue numbering, JSON emit, input parsing, dup-answer check | reuse as-is | `assign_clue_numbers/2`, `emit_json/3`, `load_clues/2`, `check_unique_answers/1` |
| `--enumerate` (count/emit all solutions) | reuse as-is | `all_crossword/5`, `count_solutions/4` |

**Explicitly out of scope** (decisions already made): `(W,H)` rectangular plumbing (square only); the `free`/size-search *loop* (`grid_candidates` retires — single `--size` stays); priority scores; barred grids. *(Correction from the design-exploration pass: `quality.pl`'s greedy **constructor** is **not** retired — it is promoted to the best-effort/candidates path, since the strict DFS cannot seed dense meshes. Only the `grid_candidates` size-search loop retires.)*

## The one new core: branch-and-bound optimizing search

*(**RESOLVED by Phase 1.5 (2026-06-30): the gate descoped this.** `prunes=0` and `best==first` across the fixtures, so the full B&B superstructure below is **not built** — the core collapses to "construct + rescore + emit" reusing the same state/objective. The B&B description is retained below for the record and in case dense-mesh work ever reopens it. See "Phase 1.5 — RESULT.")*

Today `assign_words/9` returns the *first* complete placement and strategies only reorder the tree. The new driver keeps the DFS skeleton but optimizes:

- **State:** `(RemainingWords, PlacedWords, Grid, Reward, CheckedMap)` — the current tree plus a running integer `Reward` and a per-word checked-cell map.
- **Branch:** unchanged — `select_inc` (MRV ordering, for a strong incumbent fast) → `find_intersecting_word` → `assign_word`. The legality core is untouched.
- **Objective (integer, incremental):** `Reward = WCap·Σ min(checked(w), ceil(L/2)) + WTail·Σ checked(w)`, with `WCap ≫ WTail` (default `5,1` → ε = WTail/WCap = 0.2). Integer arithmetic for deterministic, stable tiebreaks (mirrors `placement_key`'s `Crossings*10000 − Growth` idiom). Balance comes from the cap; the tail breaks ties toward crossier layouts.
- **Admissible bound (pruning):** `UB = Reward + maxRemainingReward(unplaced + under-capped placed)`; prune any branch with `UB ≤ incumbent`. Valid and cheap by construction (additive/capped objective); tightening is the optimization frontier.
- **Incumbent + anytime safety valve:** keep the best complete placement; MRV-first finds a good incumbent early; a **node/inference-count budget returns best-so-far** if the optimum is not proven in budget. *Essential* — optimization is strictly harder than first-solution search, so the contract is "best within budget," not "proven optimum." **The budget MUST be a node/inference count, not wall-clock** — wall-clock is non-reproducible and would break determinism. A wall-clock cap MAY exist only as a last-resort safety valve that, if it fires first, marks the output `truncated`. *(Golden test: byte-identical output under a fixed node budget.)*

## Phased milestones

Each phase is independently testable; `crossword.pl`/`quality.pl` keep working until the Phase 7 cutover.

| Phase | Deliverable | Touches | Exit criterion |
|---|---|---|---|
| **1. Scoring infra** | per-word `checked(w)` + capped integer reward over a *complete* placement | new `arrange.pl` | hand-scored fixtures match |
| **1.5. Search-value gate** ✅ | min. incumbent-recording B&B + greedy seed, with node/prune counters; measured whether any search beats the first incumbent | `arrange.pl` | **DONE (2026-06-30): DESCOPE** — `best==first` on 4/5, `prunes=0` on 5/5, cap inert. See "Phase 1.5 — RESULT" below. |
| **2. Strict layout (fixed N)** | *construct + rescore + emit* (B&B search descoped): first MRV-inc complete placement, rescored by `layout_reward/4`; place-all-or-fail; budget-aware 3-outcome semantics | `arrange.pl` driver | legal full placement (re-validated), rescored; reward ≥ old `solve` fill |
| **3. Sizing + emit framing** | `--size N` (default 15); `fixed` (pad N×N) vs `max` (crop `placed_bbox`), default `max` | `arrange.pl` + emit wrapper | golden output both framings |
| **4. Best-effort (drop)** | served by the **greedy constructor path** (drops naturally), *not* a drop-branch on the strict DFS; lexicographic *most-placed → reward* across greedy seeds; report dropped | `arrange.pl` + greedy | strict-fail cases degrade gracefully; dropped set correct |
| **5. Fragment seeding** | parse emit-schema fragment, reconcile by answer, pre-place + validate, search remainder (words-only v1) | `arrange.pl` + loader | emit → re-ingest → identical; conflicts reported up front |
| **6. Candidates** | distinctness from **constructor breadth (multiple deterministic seeds) + greedy diversity**, then τ-filtered (a single deterministic B&B yields near-duplicate leaves) | `arrange.pl` + greedy | `--candidates 3` → 3 layouts ≥ τ apart |
| **7. CLI + migration** | subcommand dispatch (`arrange`/`lint`/`fill`/`export`); `--enumerate`; "did you mean `arrange`?" shim; README/`run_tests.sh`/golden updates | `crossword.pl` `opts_spec`/`run` | old flag CLI removed; tests green |

## File / module layout

- **New `arrange.pl`** — driver, objective, bound, candidates, fragment loader. Consults `crossword.pl` for primitives (same pattern `quality.pl` uses).
- **`crossword.pl`** — keep all primitives; the only refactor is lifting `crossing_count`/`placed_bbox`/`word_meets_half` (currently in `quality.pl`) into the shared layer so both `quality.pl` and `arrange.pl` can call them. CLI restructure deferred to Phase 7 so nothing breaks mid-build.

## Testing & calibration

- **Determinism:** do *not* carry over `solve_shuffled`'s randomization; break ties by a stable order (lowest start cell) so golden files are stable.
- **Legality as a free property test:** re-validate every emitted layout with the existing legality predicates.
- **Reachability calibration is a *required pre-weighting step*, not an afterthought.** Before locking `WCap:WTail` (ε) and the cap target, measure per-fixture how often `ceil(L/2)` is actually reachable and **report when the objective degenerates to plain total-crossings** (probes: toc_demo 1/16 words reach it, quality_22 19/22). Make `--check-target` a real tunable; lower it below `ceil(L/2)` where half is unreachable; calibrate τ alongside. **Reachability calibration fixtures: `toc_demo` (sparse, cap inert) and `quality_22_mesh` (dense, cap active)**; the full set (`benchmark_{08,14,16,20,26}`, `quality_{22,61}`) for timing.

## Risks & mitigations

1. **Tractability (the big one):** optimal placement ≫ first-solution. → MRV-first incumbent + admissible bound + **node/inference-count budget (not wall-clock) with anytime best-so-far**. Probes show the bound prunes little at our sizes, so do not rely on it for speed — the Phase 1.5 gate decides whether the search layer is even worth building.
2. **Objective correctness:** the entire objective core (`layout_reward`, `placement_delta`, `cell_occupancy`, `upper_bound`, the `nb_setval` incumbent) is new, untested code, and the capped reward is not a plain per-placement scalar (crossings bump *crossed* words' caps). → the incremental-delta-vs-from-scratch-oracle property test (Phase 1 is the oracle) **must run on every fixture in CI** — a wrong delta silently corrupts the incumbent comparison.
3. **Best-effort is harder** (subset selection) → Phase 4, after strict is solid.
4. **Bound tightness** governs pruning power → ship the cheap valid bound first; tighten only if budgets bite.

---

## Detailed scoping — Phases 1 & 2 (the engine core)

Concrete signatures to start coding. These are a scaffold; argument shapes may shift in implementation. Placed words are the existing dict `word{answer, letters, cells, dir, len, start}` (+ `num` after numbering); a placement is a list of these.

### Phase 1 — scoring infra (the correctness oracle)

```prolog
% ceil(L/2): the per-word checking cap (== word_meets_half's threshold).
check_target(Len, Target) :- Target is (Len + 1) // 2.

% Cell occupancy across all placed words: assoc cell -> count (max 2 in a
% legal grid: one across + one down). A cell is "checked" iff count >= 2.
cell_occupancy(PlacedWords, Occ) :-
    empty_assoc(A0),
    foldl(add_word_to_occ, PlacedWords, A0, Occ).

add_word_to_occ(PW, AIn, AOut) :-
    get_dict(cells, PW, Cells),
    foldl(bump_cell, Cells, AIn, AOut).
bump_cell(Cell, AIn, AOut) :-
    ( get_assoc(Cell, AIn, N) -> N1 is N + 1 ; N1 = 1 ),
    put_assoc(Cell, AIn, N1, AOut).

% checked(w): how many of W's cells are shared (occupancy >= 2).
word_checked(PW, Occ, Checked) :-
    get_dict(cells, PW, Cells),
    aggregate_all(count, ( member(C, Cells), get_assoc(C, Occ, N), N >= 2 ), Checked).

% Per-word integer reward contribution: WCap*min(checked,target) + WTail*checked.
word_reward(WCap, WTail, PW, Occ, R) :-
    get_dict(len, PW, L), check_target(L, T),
    word_checked(PW, Occ, C),
    Capped is min(C, T),
    R is WCap * Capped + WTail * C.

% Total objective over a complete placement (the from-scratch oracle).
layout_reward(WCap, WTail, PlacedWords, Reward) :-
    cell_occupancy(PlacedWords, Occ),
    foldl([PW,Acc,Acc1]>>(word_reward(WCap,WTail,PW,Occ,R), Acc1 is Acc+R),
          PlacedWords, 0, Reward).

% default weights (ε = WTail/WCap = 0.2); make these tunable later.
arrange_weights(5, 1).
```

**Exit test:** `layout_reward/4` on a few hand-placed fixtures equals a paper calculation; it becomes the oracle that Phase 2's incremental delta is checked against.

### Phase 1.5 — search-value gate (decides whether Phase 2+ get built)

**The question:** at our densities, does *any* optimizing search produce a measurably higher capped-objective reward than the **first MRV-feasible incumbent** (or the greedy-construct seed), within a realistic node budget — or has the objective degenerated to total-crossings where one constructor pass is already near-optimal?

Build only: (1) the Phase-1 `layout_reward/4` oracle; (2) a minimal `arrange_strict` that records the first MRV incumbent (`select_inc` → `find_intersecting_word` → `assign_word`), then *continues* the DFS recording better complete placements, with the cheap loose `upper_bound` and **both a node counter and a prune counter**; (3) greedy-construct re-keyed onto the capped-objective delta (clone `quality.pl`'s loop, swap `placement_key`).

Run on `benchmark_{08,14,20,26}` and `toc_demo` (puzzle-shaped — DFS seeds in 5–462 ms; *not* the dense `quality_*` meshes, whose bottleneck is feasibility, a separate question). Record per fixture: reward of the first incumbent, reward of the best layout at a fixed node budget (e.g. 1e7 inferences), nodes explored, branches pruned by the bound, count of words where the cap binds (`checked(w) ≥ ceil(L/2)`), and the greedy-construct reward.

**Decision rule:**
- `best == first` **and** `prune_count ≈ 0` across fixtures → the search superstructure is decoration: **descope Phase 2 to "construct + rescore + emit"** (drop the B&B search loop and the deferred LNS entirely).
- `best > first` by a meaningful margin on some fixture → the B&B layer earns its place; **Phase 2 proceeds as written.**

This gate is cheap (reuses the legality core verbatim) and converts the central architectural uncertainty into a measurement instead of a bet.

### Phase 1.5 — RESULT (2026-06-30): **DESCOPE the search superstructure**

Implemented in `arrange.pl` (the `layout_reward/4` oracle + a minimal incumbent-recording B&B with node/prune counters + a greedy-seed comparison) and run on the five puzzle-shaped fixtures under an 8×10⁸-inference budget, weights 5:1:

| Fixture | grid | first | best | Δ | nodes | prunes | cap-binding | greedy (reward/placed/dropped) |
|---|---|---|---|---|---|---|---|---|
| benchmark_08 | 13 | 192 | 192 | 0 | 156,794 | **0** | 0/8 | 192 / 8 / 0 |
| benchmark_14 | 17 | 588 | 588 | 0 | 49,541 | **0** | 0/14 | 588 / 14 / 0 |
| benchmark_20 | 37 | 228 | 228 | 0 | 50,803 | **0** | 1/20 | 228 / 20 / 0 |
| benchmark_26 | 49 | 300 | 300 | 0 | 28,864 | **0** | 1/26 | 300 / 26 / 0 |
| toc_demo | 25 | 180 | 192 | **+12** | 88,214 | **0** | 1/16 | 156 / 14 / 2 |

**Findings:**
- **The bound never pruned — `prunes = 0` on all five fixtures** across ~30k–157k nodes. The admissible bound is decoration at our sizes, exactly as predicted.
- **`best == first` on 4/5.** Continuing the DFS past the first MRV incumbent found *nothing better* on the four benchmark fixtures. The only gain was toc_demo, **+12 (≈6.7%)**, unproven, after 88k nodes of bound-free search — a terrible trade.
- **The cap is inert** (`cap-binding` 0–1 of 8–26 words), so the objective is degenerating to plain total-crossings: reward ≈ 6·Σchecked. Confirms the reachability caveat on realistic inputs.
- **The first MRV incumbent equals greedy-construct** on all four all-placed benchmarks (192/588/228/300). On toc_demo the strict DFS placed all 16 (192) where greedy on the fixed grid dropped 2 (156) — i.e. greedy's edge needs *its own* grid choice, not a fixed N.

**Decision (per the gate's rule — `best == first` ∧ `prune ≈ 0`): descope Phase 2 to "construct + rescore + emit."** Drop the continue-past-first B&B search loop, the admissible `upper_bound`, the incremental `placement_delta`, and the deferred LNS. The engine becomes: **construct one complete layout (reuse the MRV-first `assign_words_inc` path for strict; greedy for best-effort/candidates) → score it with the `layout_reward/4` oracle → emit.** The oracle (Phase 1) is the keeper; the search superstructure is not.

**Consequences for the phases below:**
- **Risk 2 (delta correctness) largely evaporates** — with no incremental search there is no `placement_delta`; the oracle scores the final layout from scratch, so AC-ARR-9 reduces to "the oracle is correct," not "the delta matches the oracle."
- The marginal toc_demo gain, if ever wanted, is far more cheaply captured by **construction/seed diversity + rescore-and-pick-best** (the candidates path, Phase 6) than by exhaustive search.
- **Scope caveat:** this verdict is for the puzzle-shaped Flavour-A domain (the intended one). Dense `quality_*` meshes are a separate *feasibility* question, deliberately out of the gate's scope.
- *(The `gate_*` predicates in `arrange.pl` are the throwaway measurement harness; `layout_reward/4`, `word_reward/5`, `check_target/2`, `cap_binding_count/2` are the production scoring layer that stays.)*

### Phase 2 — strict layout (fixed N) — *descoped to construct + rescore + emit*

*(Original B&B design retained below for the record; per the Phase-1.5 result the continue-past-first search loop, `upper_bound`, and `placement_delta` are NOT built. Strict mode = the first MRV-inc complete placement, rescored by `layout_reward/4`; place-all-or-fail with the budget-aware 3-outcome semantics.)*

**BUILT (2026-06-30).** `arrange.pl` now implements Phases 2–3: `arrange_strict_solve/3` constructs the best of the four start-corner MRV-inc placements (`construct_one/7` under an inference budget), rescores with `layout_reward/4`, and emits via `emit_arrange/4`. Outcomes: `placed` | `infeasible` (names words with no possible crossing via `unplaceable_words/2`, else "grid too small") | `not_proven` (budget hit). Phase-3 framing: `fixed` → `emit_json/3` (exact N×N); `max` → `emit_cropped/6` (tight enclosing **square**, side `max(H,W)`, anchored at the bbox top-left — keeps the single-`gridLength` square schema, design-spec §6.1). Validated: valid+self-consistent JSON, all words placed, reward ≥ single first-solution, byte-identical across runs, both failure outcomes correct, existing 81/81 + golden still green. *(Formal plunit/golden tests for `arrange` are deferred to Phase 7, which owns `run_tests.sh`/golden updates.)*

Standard Prolog branch-and-bound: a mutable incumbent via `nb_setval`, exhaustive DFS that *fails* after recording each complete solution (to keep searching for better), pruned by the bound. The branch step reuses `select_inc`/`find_intersecting_word`/`assign_word` unchanged.

```prolog
% Entry: best-scoring complete placement of all Words on an N x N grid.
arrange_strict(Words, GridLen, BestPlaced, BestReward) :-
    arrange_weights(WCap, WTail),
    init_grid(GridLen, G0),
    nb_setval(arrange_best, none- (-1)),         % incumbent: Placed-Reward
    ignore(search_strict(Words, [], 0, none, empty_checked,
                         GridLen, WCap, WTail, G0)),
    nb_getval(arrange_best, BestPlaced-BestReward),
    BestPlaced \== none.                         % fail if nothing placeable

% Complete placement: record if it beats the incumbent, then FAIL to continue.
search_strict([], Placed, Reward, _St, _Chk, _GL, _WC, _WT, _G) :-
    !,
    update_incumbent(Placed, Reward),
    fail.
% Partial: prune by bound, else branch (MRV) and recurse with incremental reward.
search_strict(Words, Placed, Reward, St, Chk, GridLen, WCap, WTail, G) :-
    upper_bound(Words, Placed, Reward, Chk, WCap, WTail, UB),
    nb_getval(arrange_best, _-BestR),
    UB > BestR,                                  % else prune this branch
    select_inc(Words, Placed, St, GridLen, Start, Dir, G, Entry, RemWords, St1),
    Entry = [Word|_],
    atom_chars(Word, L0), delete(L0, ' ', Letters), length(Letters, WLen),
    find_intersecting_word(Letters, WLen, Placed, GridLen, Start, Dir),
    assign_word(Word, Letters, WLen, Start, Dir, GridLen, Placed, G, PW, G1),
    placement_delta(PW, Placed, Chk, WCap, WTail, Delta, Chk1),
    Reward1 is Reward + Delta,
    search_strict(RemWords, [PW|Placed], Reward1, St1, Chk1,
                  GridLen, WCap, WTail, G1).

update_incumbent(Placed, Reward) :-
    nb_getval(arrange_best, _-BestR),
    ( Reward > BestR -> nb_setval(arrange_best, Placed-Reward) ; true ).
```

**Strict outcome semantics (budget-aware).** A truncated search cannot prove unplaceability, so `--strict` distinguishes three outcomes: **(a)** all words placed (best within budget) → emit; **(b)** search *completed* without placing some word ⇒ that word is genuinely infeasible → fail naming it; **(c)** budget exhausted before any complete placement ⇒ fail with *"not proven; search did not complete within budget"*, **not** "word X unplaceable". The "fail naming the unplaceable word" contract (spec AC-ARR-1) applies only to case (b).

**Incremental reward delta** (the perf-critical bit; oracle = `layout_reward/4`). Placing `PW` with `k` crossings: `PW` gains `k` checked cells, and each of the `k` perpendicular partner words gains exactly one newly-checked cell (a legal cell tops out at occupancy 2). `Chk` is an assoc `answer → checked-count`.

```prolog
% Delta = WCap*min(k, T_PW)                              % PW's capped contribution
%       + WCap * (#partners still below their cap)       % each bumped partner
%       + WTail * 2k                                     % PW(k) + partners(k) raw checked
% Chk1 = Chk with PW set to k and each partner incremented by 1.
placement_delta(PW, Placed, Chk, WCap, WTail, Delta, Chk1) :- ...
% (k via crossing_count/6; partners via the crossed cells -> owning placed word.)
```

**Admissible upper bound** (cheap, loose, valid):

```prolog
% Max extra reward still achievable from here, optimistically:
%   unplaced word u : WCap*T_u + WTail*L_u
%   placed P below cap (checked c<T): WCap*(T_P - c) + WTail*(L_P - c)
% UB_total = Reward + that sum. Prune when UB_total <= incumbent.
upper_bound(Unplaced, Placed, Reward, Chk, WCap, WTail, UB) :- ...
```

**Seeding note:** absolute position of the first word is immaterial under `max` (cropped) and, for `fixed`, the search + `fits_on_grid` handle fit — so seed the first word once at a canonical interior position (reusing `start_loc` is fine) rather than branching over all four corners as `solve` does. Trying multiple seeds is an optional breadth knob, not required for correctness.

**Exit test:** for each benchmark fixture, `arrange_strict/4` returns a placement that (a) re-validates as legal, (b) places all words, and (c) has `layout_reward` ≥ the reward of the old `solve`'s arbitrary fill — and the incremental `Reward` it carries equals `layout_reward/4` recomputed from scratch (delta correctness).

### After Phase 2

Phases 3–7 hang off this core: framing is an emit-time wrapper (`max` crops via `placed_bbox`); **best-effort and candidates ride the greedy constructor path** (Phases 4/6 above) — best-effort uses the lexicographic *(count, reward)* incumbent across greedy seeds, candidates τ-filter the seed breadth — *not* a drop-branch or top-K buffer on the strict DFS; fragments pre-seed `Placed`/`Grid` and shrink `Words`.

**Deferred enhancement — LNS polish.** A destroy-and-repair (large-neighbourhood) pass over the incumbent is a *deferred* option, built **only if** the Phase-1.5 gate shows measurable reward headroom over the first feasible / greedy layout. The cap-inertness data predicts little headroom, so this is not on the critical path.
