# Fill performance program closeout (2026-07-16)

Closeout for [`../plans/fill-perf-program-2026-07.md`](../plans/fill-perf-program-2026-07.md).
The program is complete. One contract-safe engineering arm landed; every
search-behavior arm remained on a probe and failed its registered basket.

## Outcome

| item | verdict | decisive evidence |
|---|---|---|
| A1 capacity reporting | resolved before this campaign | DP-10 / AC-FILL-15 clean failures; envelope unchanged |
| A2 seeded load walk | **won and built** | exact-sequence Fenwick replay, 171.14s → 9.57s (17.88×) |
| B0 instrumentation | measured | long queues, 28% fruitful revisions, concentrated weights, >99% propagation wall |
| B1 revision ordering | lost | large STW wins, catastrophic `g17_50k` trajectory regressions |
| B2 alternate credit | lost | broad credit flattened weights; H2 overfit and failed quality/ladder |
| B3 aging/probing | lost | width-specific gains, authority/quality/ladder failures |
| C1 projection classes | lost | 1.00023× target compression, 0/5 target completions, three ladder losses |
| Track D | **not triggered** | no B/C arm cleared completion + quality + ladder gates |

## What Changed

`seeded_permutation/2` now replays the historical `V mod remaining` selection
sequence through a local Fenwick tree in O(n log n). Five seeds across lengths
0..256 matched the old permutation and following PRNG draw (1,285 cases).
The seeded CWL trivial-grid reproducer improved 17.88×; peak RSS rose 4.26%.
No seed semantics, fills, or contracts changed.

The B0 exact-replay instrumentation rig and its corrected evidence record were
landed under `benchmarks/fill_quality/` and `benchmarks/results/`. Search
variants themselves were not landed; their commits remain referenced from the
append-only experiment ledger.

## What Was Learned

The search is not short-queue or signal-starved. Authority worklists are long,
fruitful revisions are common, and baseline conflict weights are sharply
concentrated. The cost is overwhelmingly bignum propagation/support. That made
the thesis candidates plausible, and several produced very large local wins.

The failure mode was generality. Queue ordering, deleted-count credit, aging,
and probing each changed learned weights enough to redirect deterministic
restart trajectories. An arm could cut the STW reference by 40-90% and still
lose `g17_50k`, quality, or the second reference threshold. Broad credit was
worse: it smeared the sparse signal dom/wdeg needs.

The envelope hypothesis also failed at its premise. UK unchecked cells did not
produce large *exact* joint checked-letter classes in visited domains. With
classes effectively singleton, deferred all-different was pure cost and no
blocked target completed.

## Recommendation

Keep the locked §8.4c MAC + dom/wdeg + restart engine unchanged. Keep A2's
identity-clean load improvement. Do not adopt or combine F1-F6 from this pass,
do not change the budget/default dictionary/min-score posture, and do not
reopen `blocked_13b`/`blocked_15a` without a genuinely new mechanism.

The only named revisit triggers are a concrete need for segmented bignum masks
(the DP-10 capacity envelope), a new representative benchmark corpus that
changes the cross-rung verdict, or a sound all-different/class propagator paired
with measured non-singleton projection classes. None is present now.

## Explicit Non-Actions

- No Track D adoption decision was taken; there was no qualifying proposal.
- No search strategy was threaded through the shipped `mac_*` path.
- No golden, `fill_identity.sha256`, inference baseline, or history row was
  regenerated to accommodate a changed fill.
- No STW data was committed; the CC BY-NC-SA snapshot stayed local.
- No budget increase, min-score claim, default-dictionary change, approximate
  value clustering, or further chase of the DP-6 blocked rows was attempted.
