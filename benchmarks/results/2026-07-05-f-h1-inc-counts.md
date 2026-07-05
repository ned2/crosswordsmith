# F-H1 — incremental candidate counts for MRV (2026-07-05)

Phase 2 lead experiment of the fill performance campaign
(`docs/plans/fill-perf-campaign.md`). Branch `experiment/f-h1-inc-counts`
(base d38a649, the post-F-L1 baseline). `prolog/crosswordsmith/fill.pl` is the
only engine file changed; goldens, fixtures, `fill_baseline.json`,
`fill_identity.sha256` untouched. Machine: host `lana x86_64`, SWI-Prolog
10.1.10 (same as the baseline). Inference counts deterministic; wall/rss
host-specific.

## Hypothesis + mechanism

`select_mrv/6` recounted EVERY unfilled slot at EVERY search node
(`maplist(slot_candidate_count, Slots, ...)`), measured by P-F1 at 79.5-85.6%
of `search_inf` on all four probed rungs. But a slot's candidate count is a
pure function of its cells' binding state, and placing a word binds ONLY the
previously-free cells of the placed slot. So carry the per-slot counts as
backtrack-restored threaded state and, after each placement, recount EXACTLY
the slots crossing a NEWLY-BOUND cell; every other count carries over unchanged.

## Results — all 11 rungs (3× reproduced, count tables byte-identical)

`search_inf` GATED, `load_inf` GATED (both at 0.5% relative). RESULT: PASS
(11 improvements, 0 regressions) on both the core pass and the `--heavy` pass.

| rung | search_inf before | search_inf after | Δ% | load_inf Δ% | cmd_wall ms before→after |
|------|------------------:|-----------------:|----:|:-----------:|:------------------------:|
| sq04_full     |    529,275 |    371,343 | **−29.84%** | +0.00% | 2660 → 2620 (−1.50%) |
| g11_full      |    671,484 |    514,725 | **−23.35%** | +0.00% | 2700 → 2660 (−1.48%) |
| g11_full_seed |    699,377 |    419,582 | **−40.01%** | +0.00% | 2600 → 2570 (−1.15%) |
| sq05_full     |  1,461,408 |    876,182 | **−40.05%** | +0.00% | 2750 → 2710 (−1.45%) |
| g17_full      |  2,534,409 |  1,822,010 | **−28.11%** | +0.00% | 2880 → 2660 (−7.64%) |
| g21_full      |  3,307,580 |  1,853,583 | **−43.96%** | +0.00% | 2970 → 2860 (−3.70%) |
| g13_full      |  3,754,832 |  2,934,984 | **−21.83%** | +0.00% | 2910 → 2660 (−8.59%) |
| sq04_50k      |  7,738,070 |  6,062,379 | **−21.66%** | +0.00% | 1000 →  940 (−6.00%) |
| g15_full      | 10,734,114 |  7,967,168 | **−25.78%** | +0.00% | 3115 → 2870 (−7.87%) |
| g17_50k       | 19,637,890 | 13,773,804 | **−29.86%** | +0.00% | 2005 → 1255 (−37.41%) |
| g09_full      | 34,880,750 | 14,870,446 | **−57.37%** | +0.00% | 4110 → 3305 (−19.59%) |

- **`load_inf` +0.00% on all 11 rungs** (the load path is untouched; gated).
- **The top/deep rungs win biggest**: g09_full (the hardest rung, 9,961 nodes)
  −57.37%; g21_full (154 slots, depth 153) −43.96% — matching the
  pre-registered −30..−70% band. g17_50k −29.86% sits at the band edge.
- **Wall improves MORE than `search_inf` on the search-dominated rung**:
  g17_50k `search_inf` −29.86% but `cmd_wall` −37.41% — the predicted signature
  of also deleting the inference-blind `length/2` bucket walks (P-F1 §D: 28.5%
  of that rung's search WALL, ~0 `search_inf`). On dict-load-dominated rungs
  (g09_full: load = 68% of wall) the wall gain (−19.59%) is smaller than the
  `search_inf` gain, as expected.

## Design

Carried structure: `Counted` = a list of `cnt(Count, Slot)` in slot-list order,
each `Slot` the ORIGINAL `slot(Start,Dir,Cells,Vars)` term (shared crossing
variables intact). Built once at the root by `maplist(slot_with_count, ...)`
(the one full recount), then threaded down the recursion.

Per node (`fill_search_inc/4`):
1. `select_min_count/3` — the `cnt/2` with the minimum `c(Count,Start,Dir)` in
   standard order (a fold, not a sort), and `Rest` = the others.
2. `newly_bound_cells/3` — the placed slot's cells whose var is still unbound
   (PRE-placement): exactly the cells this placement will newly bind. A plain
   list of ground cell NUMBERS.
3. `candidates/4` on the winner (unchanged) → try each `Word` in dict order.
4. After `BestVars = Word`: `recount_crossing/5` walks `Rest` and recounts (via
   the unchanged `candidate_count/4`) EXACTLY the slots whose `Cells` intersect
   the newly-bound cell numbers; every other `cnt/2` carries its count forward.
5. Recurse on the updated `Rest1`.

Affected-slot detection is by CELL NUMBER (ground integer) identity, not
variable identity — a slot crosses the placement iff it shares a newly-bound
cell number. This sidesteps the findall/copy trap entirely: no cell variable is
ever inspected for identity or copied. `remove_slot/4` drops the winner by its
ground `(Start,Dir)` key compared with `==` (no binding), preserving list order
and term-sharing.

Backtracking restores counts for free: `Rest`/`Rest1` are fresh list spines per
level; unwinding discards them and returns to the parent's `Counted` with its
counts intact. No `nb_*`/`assert`/`setarg` — pure backtrack-restored state.

`select_mrv/6` and `slot_candidate_count/4` are RETAINED unchanged as the
reference selector (the white-box tests R6/P13 call `select_mrv/6` directly);
they are simply no longer on the engine's hot path.

### Affected-slot-count distribution (why the win is bounded, not ~90%)

Static crossings-per-slot (each slot's perpendicular crossings — the upper
bound on affected-per-placement, since affected ≤ newly-bound cells ≤ slot
length ≤ 5):

| grid | slots | crossings/slot min/max/mean |
|------|------:|:---------------------------:|
| 04a | 8 | 4 / 4 / 4.00 |
| 05a | 10 | 5 / 5 / 5.00 |
| 09a | 32 | 3 / 5 / 3.75 |
| 11a | 44 | 3 / 5 / 3.73 |
| 13a | 56 | 3 / 5 / 4.18 |
| 15a | 74 | 3 / 5 / 3.92 |
| 17a | 100 | 3 / 5 / 3.82 |
| 21a | 154 | 3 / 5 / 3.86 |

So each placement recounts ≤5 of up to 154 slots (vs the old full recount of
all remaining). The win is not ~95% because two costs are NOT removed: the
winner's own `candidates/4` materialization (P-F1: 12.7-17.8% of `search_inf`)
and the new O(|Rest|) `recount_crossing` scan that replaces the counting. These
set the floor; the deeper/denser rungs (more redundant recounts per node in the
old code) realize the largest cuts.

### Variants tried

- **`recount_crossing` term reuse** (reuse the input `cnt/2` for unaffected
  slots instead of rebuilding it): `search_inf` byte-identical (1,853,583 on
  g21_full — SWI counts the unification identically), and it did NOT move the
  g21_full RSS. No measurable benefit, so REVERTED to the plain rebuild for the
  simpler diff (the version validated by the identity gate).

## Equivalence argument + evidence

The tree is preserved node-for-node:
1. **Counts exact by construction.** Root counts are a full `candidate_count`
   per slot; thereafter a slot untouched by a placement keeps its exact count
   (candidate_count is a pure function of the cells' binding state) and a
   crossing slot is recounted by the same `candidate_count/4`. Counts NEVER read
   `Used` (the `\+ memberchk(Word, Used)` filter is try-time only,
   `fill_search_inc`), so carried counts are exact regardless of placement
   history — verified in-code.
2. **Same selection.** `select_min_count` picks the minimum `c(Count,Start,Dir)`
   in standard order, identical to `select_mrv`'s `sort(0, @=<, ...)`+head.
   `Start`+`Dir` uniquely identify a slot ⇒ total order, no ties, order-
   independent — so the identical winner is chosen at every node.
3. **Same winner materialization** via the unchanged `candidates/4`.
4. **Completed slots stay in the set** — a slot fully bound by crossings
   (count 0 or 1) is not dropped, so a 0-count dead slot is still selected first
   (lowest count) and fails the branch exactly as today.

Evidence:
- **Identity oracle 11/11 OK** (`benchmarks/check_fill_identity.sh`): every
  rung's CLI stdout is byte-identical to `fill_identity.sha256`. Re-run on the
  final committed bytes — still 11/11.
- **Full test suite green** — 204 passed / 0 failed, including the `select_mrv`
  white-box pins (R6 `select_mrv_recovers_correct_direction_on_tie`, P13
  `select_mrv_leaves_no_choicepoint`), P3 `candidate_count_matches_candidates`,
  the shared-cell-variable regression, both fill goldens (3x3 + fill-15 bench).
- **Clean load** — no singleton/discontiguous/warning.

## Reproduction (3× identical)

`check_fill_baseline.pl --heavy` run three times on the final bytes; the gated
count tables (`search_inf` measured + `load_inf` measured for all 11 rungs) are
byte-identical across all three (diff clean). RESULT: PASS all three.

## Files

- `prolog/crosswordsmith/fill.pl` — `fill_search/4` reworked to the incremental
  drive; new `slot_with_count/4`, `fill_search_inc/4`, `select_min_count/3`,
  `min_count_walk/3`, `count_le/2`, `remove_slot/4`, `newly_bound_cells/3`,
  `recount_crossing/5`, `shares_cell/2`. `select_mrv/6` + `slot_candidate_count/4`
  retained unchanged (reference selector / white-box tests). +105 / −5 lines.
- `benchmarks/results/2026-07-05-f-h1-inc-counts.md` — this doc.

## Risks / anomalies

- **RSS +82% on g21_full only** (278,020 → ~507,336 KiB, stable across all
  runs). g21_full is the uniquely deep rung (max backtrack depth 153); the
  threaded counted-list is rebuilt per placement (`recount_crossing` rebuilds
  the spine) and, under depth-153 backtracking, raises the process global-stack
  high-water mark. All other full-ENABLE rungs stay ~278-308 MiB; g09_full
  (depth 31, 9,961 nodes) does NOT move — so it correlates with DEPTH, not node
  count. RSS is reporting-only and NEVER gates (whole-process peak footprint,
  not a search-memory metric); it is bounded (not a leak) and 507 MiB is modest
  natively. The term-reuse variant did not reduce it (the churn is spine
  rebuild, not `cnt/2` wrappers). Mitigation IF footprint ever matters (e.g. the
  Phase 3 WASM posture): carry counts in a `(Start,Dir)`-keyed assoc (O(log n)
  path-copy per update, no spine rebuild) — a redesign that changes inference
  counts and needs full re-validation, hence out of scope for this gated search
  win. Flagged for the orchestrator.
- No `search_inf` rise on any rung (the fail condition); all 11 are wins.

## Verdict recommendation

**ACCEPT.** Byte-identical output (11/11 identity), all tests green, `load_inf`
+0.00% (gated), `search_inf` −21.66% to −57.37% across all 11 rungs (0
regressions), 3× deterministic. The one anomaly (g21_full RSS) is a
reporting-only footprint on the deepest rung with a clear mechanism and a
scoped mitigation path. Recommend re-baseline (`--record`) at the merge commit
and history append, per campaign policy.
