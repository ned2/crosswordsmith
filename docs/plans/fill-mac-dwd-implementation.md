# Implementation plan: §8.4c `fill` search core (MAC + dom/wdeg + restarts)

**Status: IN PROGRESS (started 2026-07-16).**
Spec: design-spec **§8.4c [LOCKED — DP-8]** (plus the §8.4b seed-semantics
amendment). Evidence base: the DP-7 MAC-probe campaign in
[`../../benchmarks/fill_quality/README.md`](../../benchmarks/fill_quality/README.md)
(dom/wdeg + determinism follow-ups, 2026-07-16); reference implementation of
the recipe: [`probe_mac.pl`](../../benchmarks/fill_quality/probe_mac.pl)
(`probe_mac_dwd/7` — greedy mode is the default-path blueprint).

## What changes and what must not

The §8.4 MRV/backtracking tree is replaced by the §8.4c core: bignum
candidate-mask domains over the §8.4a score-ordered buckets, AC-3
propagation over the crossing graph, dom/wdeg slot selection with global
conflict weights, growing-cap restarts with weight aging, greedy (RNG-free)
candidate picks by default and seeded top-3 under `--seed`/`--shuffle`.

Preserved contracts (the reason this is a revision, not a new engine):

- **Outcomes** (AC-FILL-1/14): `filled` / `infeasible` (proof only:
  uncapped exhaustion or root-wipeout, naming the slot) / `not proven
  within budget`. Non-zero exit + no `--out` file on the failure arms.
- **Budget** (AC-FILL-9): inference-count semantics via
  `call_with_inference_limit/3`, default 800M, `--budget N` override.
- **Quality** (§8.4a, AC-FILL-5/12): score-desc-then-dictionary candidate
  order is structural (ascending-bit materialization over score-sorted
  buckets); `--min-score` prune unchanged (it happens at load).
- **Determinism** (AC-FILL-3/13): default path consults NO RNG; seeded path
  uses only the engine-internal splitmix64 (CLI/WASM parity).
- **Seeds/pins** (AC-FILL-2): seeded slots are excluded from search exactly
  as today; their bound letters constrain crossing slots' *initial domains*.
- **Report/emit**: `slots_to_layout`, quality report, JSON contract — all
  downstream of "slot Vars are bound", which the new core still guarantees
  (final placement binds each word onto the shared cell variables; this is
  also the built-in consistency check the probe used).

Accepted one-time churn (DP-8): every fill changes once → regenerate fill
goldens, re-baseline `benchmarks/fill_identity.sha256` and the fill perf
ratchet (`benchmarks/fill_baseline.json`), update affected plunit
expectations and fuzz outputs.

## Stages

- **C1 — spec + plan (this commit).** §8.4c, §8.4b amendment, DP-8, §8.5
  row, this plan.
- **C2 — engine core.** In `fill.pl` (house style: fill machinery stays in
  the module):
  - Setup per `fill_attempt`: crossing edges with edge ids from slot cells;
    flat letter-mask table (`lm/26` per Len-Pos) from the existing
    `build_masks/2` output; per-length bucket compounds; initial domains
    (bucket mask ∩ bound-letter masks for seed/crossing-constrained cells);
    weight compound (destructive `nb_setarg` global — document the global
    key and reset it per run so runs are independent).
  - Search: `restart_loop` (cap 500 ×1.5, weights aged ×0.99, attempts
    unbounded — the inference budget is the outer stop), dom/wdeg selection,
    greedy materialization by default, seeded top-3 twin gated on
    `current_search_seed/1` (read once per search, §8.4b house pattern).
  - Wire-in: `fill_search/5` dispatches to the new core; the old
    `fill_search_inc`(+`_seeded`) are RETAINED as reference predicates
    (select_mrv precedent) for white-box tests and documentation.
  - Duplicate-answer filter (`Used`, incl. seed words) as today.
  - Outcome plumbing: distinguish cap-abandon (→ keep restarting) from
    true exhaustion (→ infeasible) from inference-limit trip (→ not_proven).
    NEVER report infeasible off a budget/cap stop (AC-FILL-14).
- **C3 — verification layer.**
  - plunit: new white-box suites for edges/domains/dom-wdeg selection/
    restart determinism; update tests that pinned old fills.
  - `make update-golden` for fill goldens; REVIEW the diff (fills should
    change but stay valid + quality-comparable; run `lint` sanity).
  - Identity oracle + ratchet: regenerate, review magnitudes, commit as the
    new baseline (document before/after inference counts here).
  - Fuzz (`tests/determinism_fuzz.sh`): all det cases must hold (greedy
    path 3-process byte-identity; seeded cases reproduce per seed).
  - `make test-wasm`: bignum masks must behave under the LibBF (USE_GMP=OFF)
    backend; seeded parity CLI↔browser.
  - `benchmarks/fill_quality/`: re-run the harness gate (AC-FILL-12) —
    native completes `blocked_13a` @ 30 and @ 1 within the default budget,
    quality ≥ old engine on all completable grids.
- **C4 — docs sweep.** README (limitations bullet: the gap row is CLOSED;
  flags table wording), STATUS, benchmark README result columns.

## Measured facts to pin during C3 (do not hand-wave)

1. **Default-budget completion of the reference row on the greedy path.**
   The probe burned ~126s @ 30 greedy — likely > 800M inferences. If the
   default budget does not cover it: decide budget default vs documented
   `--budget`/`--seed` guidance, and record the choice + numbers HERE and
   in §8.4c (AC-FILL-12 wording allows either only if measured and
   documented; do not silently regress the AC).
2. **Easy-grid latency + mask-build overhead** (9×9/enable_25k and a
   full-dict run): the probe was FASTER than the old engine on the 9×9
   (0.45s vs 1.4s) but mask/edge setup on a 315k dict costs ~2s — measure
   end-to-end CLI wall times before/after.
3. **blocked_13b / blocked_15a**: still expected NOT to complete (they
   defeat ingrid too). Confirm they exhaust budget cleanly with the right
   outcome (not_proven), and record times.

## Rollback

Single-commit-revertible per stage; the old search predicates remain in the
tree (reference status), so a revert of C2/C3 restores the old engine
byte-identically (goldens/oracle/ratchet revert with it).
