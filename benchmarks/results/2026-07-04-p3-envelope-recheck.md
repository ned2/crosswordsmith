# P3 — arrange density-envelope re-check after E-H9 (2026-07)

Product question: E-H9 (check-only counting legality + O(1) boundary grid,
commit fb23db5) cut search inferences a further **-12.9% to -50.3%** on the
12 ladder rungs. Do any of P2's *near-cliff* instances — the ones that were
budget-marginal or failed under the shipped 500M default — now FLIP to a
comfortable placement, and does the density envelope move?

Method: identical to P2.
`crosswordsmith_arrange:arrange_best_layout(Words, G, 500_000_000, _, Reward, Outcome)`
under `call_time`, each in a SEPARATE `timeout 300 swipl --stack-limit=4g`
process. Fixtures regenerated with `benchmarks/gen_mesh_fixture.py G N K 3 4 SEED`
(9x9 K=5, 15x15 K=6, 21x21 K=6) — byte-identical to P2's (verified by diff on
several instances). Budget = 5e8 = shipped default. Inference counts are
deterministic/machine-independent (only wall varies). Only P2's marginal and
failed instances are re-run; instances P2 placed solidly were not re-run (they
only got cheaper, they were never at risk of flipping *away* from placed).

## Re-check: P2 vs P3 (all @ 5e8 shipped default)

| size | N | seed | P2 outcome | P2 inferences | P3 outcome | P3 inferences | P3 wall (s) | verdict |
| --- | ---: | ---: | --- | ---: | --- | ---: | ---: | --- |
| 9x9  | 18 | 11 | placed (near ceiling) | 341,632,830 | **placed** | 275,166,940 | 15.31 | placed, -19.5% |
| 9x9  | 18 | 12 | not_proven | 500,000,213 | not_proven | 500,000,215 | 29.09 | still fails |
| 9x9  | 18 | 13 | placed (easy) | 2,022,986 | **placed** | 1,645,090 | 0.10 | placed, -18.7% |
| 9x9  | 19 | 11 | not_proven | 500,000,218 | not_proven | 500,000,217 | 29.90 | still fails |
| 9x9  | 19 | 12 | not_proven (gen planted only 18w) | 500,000,218 | not_proven (gen planted only 18w) | 500,000,212 | 29.20 | recurs: generator caps this seed at 18 words |
| 9x9  | 19 | 13 | placed MARGINAL | 500,000,513 | placed MARGINAL | 500,000,515 | 31.91 | still marginal (pinned at ceiling) |
| 9x9  | 20 | 11 | not_proven | 500,000,217 | not_proven | 500,000,215 | 30.36 | still fails |
| 15x15 | 44 | 11 | not_proven | 500,000,243 | not_proven | 500,000,242 | 30.42 | still fails |
| 15x15 | 44 | 12 | not_proven | 500,000,243 | not_proven | 500,000,244 | 28.90 | still fails |
| 15x15 | 44 | 13 | placed | 95,976,414 | **placed** | 79,042,760 | 4.60 | placed, -17.6% |
| 21x21 | 84 | 13 | not_proven (84w fails this seed) | 500,000,287 | not_proven | 500,000,282 | 28.49 | still fails |
| 21x21 | 88 | 11 | not_proven | 500,000,287 | not_proven | 500,000,287 | 28.88 | still fails |
| 21x21 | 88 | 12 | placed MARGINAL | 500,001,037 | placed MARGINAL | 500,001,036 | 29.61 | still marginal (pinned at ceiling) |
| 21x21 | 88 | 13 | placed MARGINAL | 500,001,030 | placed MARGINAL | 500,001,026 | 30.91 | still marginal (pinned at ceiling) |
| 21x21 | 92 | 11 | not_proven | 500,000,288 | not_proven | 500,000,286 | 30.34 | still fails |

**Nothing flipped.** No P2 hard failure became placed; no P2 budget-marginal
placement became robust. E-H9 splits the re-check set cleanly in two:

- **Solution found within budget** (the three that placed with headroom in P2):
  full E-H9 discount, -17.6% to -19.5% inferences — dead on the ladder's
  envelope-rung range (17w -17.5%, 40w -22.7%). These were never at risk.
- **Budget-saturated instances** (every not_proven, and the three MARGINAL
  placers): reproduce their P2 inference count to within ~10 inferences
  (19w s13 500,000,513->515; 88w s12 500,001,037->036; 88w s13 500,001,030->026;
  the failures all still land at ~500,000,2XX). E-H9's per-node speedup gives
  them **zero** movement.

## Why the near-cliff regime is immune to E-H9 (mechanism)

`construct_corners/7` shares ONE 5e8 budget across the 4 start corners,
decrementing the running budget by what each corner actually consumes. A
near-cliff instance has at least one *search-bound* corner that never completes
within any sane budget (P2 showed these stay not_proven at 2e9 = 4x). That
corner eats its whole slice of the running budget regardless of per-node cost:
E-H9 lets it visit MORE nodes per inference, but it still doesn't finish, so it
still spends the full slice. The total is therefore pinned at the 5e8 cap, and
whether the instance places at all depends only on whether a *productive* corner
gets to run before the budget is drained — a corner-ordering accident that E-H9
leaves unchanged here. Hence: robustly-placed instances (solution reached inside
one corner's slice) get the full discount; budget-saturated instances get none.

## Revised envelope per size (same terms as P2) — UNCHANGED

Because nothing flipped, P2's envelope statement stands verbatim; only the
*cost* of the instances that already placed within budget dropped ~18%.

| size | robust envelope (places across seeds, well under 5e8) | first fully-noisy N | first all-seed cliff |
| --- | --- | --- | --- |
| 9x9  | 17w (robust; ~46.7M in P2, ~38M post-E-H9) | 18w — 2/3 place (s11 now 275M, was 341M; s12 fails) | 19w (s13 marginal at ceiling, s11/s12 search-bound; 19w s12 the generator can only plant 18w) |
| 15x15 | 42w (all 3 seeds) | 44w — 1/3 places (s13 only, now 79M, was 96M; s11/s12 fail) | ~44-46w |
| 21x21 | 82w (robust, ~10.8M) | 84w — 2/3 place (s11 placed, s12 marginal in P2; s13 fails) | 88w — 2/3 place but MARGINAL (s12/s13 pinned at ceiling), s11 search-bound |

No N-level went all-3-seeds-placed under the re-check, so no opportunistic
next-N probes were warranted (per method: only probe N+2 if a level fully flips).

## (c) Does the search-bound conclusion still hold?

Yes — and P3 STRENGTHENS it. P2 established "raising the budget is not the lever"
(seed-11 first-failures stay not_proven at 4x budget). P3 adds the complement:
"**making each node cheaper is not the lever either.**" A -13% to -50% per-node
inference cut — the largest single-experiment win of the campaign — moved the
near-cliff envelope by exactly zero words. The near-cliff regime's cost is
dominated by search-bound corners that saturate whatever budget (and whatever
per-node cost) they are handed. The right future play is unchanged: opt-in
randomized restarts on not_proven for this heavy-tailed regime
(docs/research/arrange-search-algorithms.md #3) — a product feature, not a
search-speed change.

## (d) Surprises

1. **Version-invariant budget saturation.** The three MARGINAL placers came back
   with inference counts within ~10 of their P2 values (500,000,513->515;
   500,001,037->036; 500,001,030->026) despite E-H9 restructuring the entire
   per-node counting path. This is direct evidence the count on these instances
   is pinned by the budget mechanism, not by how deep the solution sits — the
   "marginal at ~5e8" label is a property of the *budget wall*, not the instance.
2. **Clean bimodal E-H9 payoff.** On this near-cliff set E-H9 is all-or-nothing:
   full ~18% discount if the solution is inside budget, precisely 0% if the
   instance is budget-saturated. This is the flip side of P1's finding that
   counting dominates the *backtrack-free* rungs (which is exactly where E-H9's
   big numbers came from) — it does not help the backtrack-*bound* tail.
3. **P2's generator caveat recurs deterministically.** `gen 9 19 5 3 4 12`
   again planted only 18 words (grid near geometric saturation), so 9x9 19w s12
   is really the 18w s12 instance and stays not_proven — recorded as-is.

The three P2 envelope-guard rungs (ladder_09x09_17w, ladder_15x15_40w,
ladder_21x21_82w) remain the right ratchet choices: all place with headroom, so
they took the full E-H9 discount and stay deterministic; none of the marginal
near-cliff instances are promotable (they would flip on any budget-adjacent
change, as their pinned-at-ceiling counts here re-confirm).
