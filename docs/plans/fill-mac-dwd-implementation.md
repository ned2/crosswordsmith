# Implementation plan: §8.4c `fill` search core (MAC + dom/wdeg + restarts)

**Status: DONE (built 2026-07-16; C1 spec/plan 625563c, C2 core 1ffc881,
C3 search-power amendments + re-baselines ce8e02d, C4 docs sweep in the
closing commit). All four stages complete; the C3 findings below amended
the locked spec (§8.4c points 4–6, AC-FILL-13) — read them before touching
the pick policy or load ordering.**
Spec: design-spec **§8.4c [LOCKED — DP-8]** (plus the §8.4b seed-semantics
amendment). Evidence base: the DP-7 MAC-probe campaign summarized in
[`../../benchmarks/fill_quality/README.md`](../../benchmarks/fill_quality/README.md)
(dom/wdeg + determinism follow-ups, 2026-07-16). The historical
`probe_mac_dwd/7` reference implementation (greedy mode is the default-path
blueprint) is recoverable at commit `4653996` via
[`../research/benchmark-probe-historical-reconstruction.md`](../research/benchmark-probe-historical-reconstruction.md).

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
- **Determinism** (AC-FILL-3/13, as amended in C3): default path is a pure
  function of the input — no OS entropy, no run-to-run variance; it drives
  the engine-internal splitmix64 from PINNED constants (load shuffle +
  attempt ≥ 2 top3), and the seeded path uses the same PRNG on the user's
  seed (CLI/WASM parity).
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

## C3 findings — the two search-power amendments (2026-07-16)

The C2 core as spec'd (greedy default, top3 only under a user seed, eager
per-node domain materialization) passed plunit/goldens but FAILED the
identity ladder's open-grid rungs (`g17_full`, `g21_full`: budget-dead
under the default 800M; old MRV engine: ~2.5M/~3.3M inferences). The C3
campaign isolated three interacting causes and landed three fixes, all in
`fill.pl` (+ `prng_seed_state/1` in `core.pl`):

1. **Lazy candidate enumeration** (`mac_candidate/4` replacing
   materialize-then-`member/2`): the eager form cost
   O(popcount × mask-width) bignum work per node — g17_full burned 800M
   inferences in ~4k nodes. Lazy `lsb`-walk enumeration made the seeded
   g17 run collapse from 7m12s to seconds.
2. **Diversified default-path restarts** (`mac_attempt_pick/3`): greedy-only
   restarts never converged on g17_full (23 attempts / ~1M nodes / 3.7M
   node-cap — weight aging alone re-treads the same prefix), while top3
   filled it on attempt 2. Default now: attempt 1 greedy, attempts ≥ 2
   top3 on splitmix64 reseeded from pinned constant 0xCC9E2D51.
3. **Pinned load shuffle for order-free dicts** (`seed_perturb_plain/2`,
   single-band arm of `seed_tie_shuffle/2`, constant 0x9E3779B9): top3
   over a lexicographic prefix is an alphabetical clump and diversifies
   nothing — g17_50k died under all-top3 AND under a greedy/top3
   alternation; with the shuffle every ladder rung fills in seconds.
   Measured BOTH directions: scored multi-band dicts keep strict
   score-desc-then-lex (in-band permutation starved the reference row —
   pinned-shuffle run @30: not-proven at 3m30s; the DP-7 probe's fast 3/3
   ran top3 over lex bands), so the shuffle applies exactly where the
   order carries no information (plain, or all scores equal — which is
   also what keeps AC-FILL-6's uniform-dict identity).

Also fixed en route: `pinned_score_group`'s O(n²) trap avoided — the
pinned shuffle is an O(n log n) draw-key keysort (`prng_shuffle/2`), NOT
`seeded_permutation/2`'s selection walk (minutes at ENABLE scale; the
seeded path keeps the O(n²) form because its draw sequence is DP-6
contract). Toy-scale quality note: on the two-square 3×3 scored fixture
the engine now lands on the mixed square (per-node band preference never
guaranteed the global optimum; the old lex order just happened to walk
into it) — the reference-row mean is the real quality gate and it IMPROVED
(45.0 vs the probe's 44.81 and ingrid's 44.4).

## Measured facts pinned during C3

1. **Reference row (AC-FILL-12), default path, default 800M budget**
   (blocked_13a × STW): `--min-score 30` **filled in 2m20s / 7 attempts,
   mean 45.0 min 30**; `--min-score 1` **filled in 19s, mean 38.7**. No
   budget-default change needed. (`--seed` variance is real: seed 7 on
   @30 budget-exhausted at 8m47s with a correct `not proven` non-zero
   exit; seeds are the variety lever, not the completion path.)
2. **Ladder / easy-grid latency** (end-to-end CLI wall, full ENABLE):
   15×15 3.2–3.6s (attempt 1), 17×17 3.5s (attempt 1), 21×21 3.9s
   (4 attempts), 17×17 × 50k 1.3s (4 attempts). Identity oracle: 11/11
   rungs recorded and verified. Ratchet re-baselined
   (`fill_baseline.json`; heavy pass): small/blocked rungs rose to a
   ~3.3–3.8M-inference floor (mask/edge setup dominates), while the
   rungs the old engine found hardest DROPPED (g09_full 14.5M → 3.7M,
   g17_50k 13.0M → 9.8M, g15_full 7.9M → 3.7M).
3. **blocked_13b / blocked_15a**: not re-run (DP-6 report-don't-chase pin
   stands; they defeat ingrid too and completion there was never promised
   by DP-8). AC-FILL-14's outcome mapping is regression-tested at unit
   level instead (`mac_exhaustion_is_infeasible_not_retried`).
4. **Full verification battery at the final build**: plunit 379/379,
   goldens (only `fill_3.json`/`fill_15_bench.json` changed vs the C2
   commit; seeded goldens byte-identical — the lazy top3 is draw-for-draw
   equivalent), `make test` ALL PASSED, fuzz 66/66, identity 11/11,
   ratchet PASS, `make test-wasm` OK (bignum masks + pinned streams
   bit-identical under LibBF).

## Rollback

Single-commit-revertible per stage; the old search predicates remain in the
tree (reference status), so a revert of C2/C3 restores the old engine
byte-identically (goldens/oracle/ratchet revert with it).
