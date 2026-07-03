# P2 — arrange density envelope probe (2026-07, post E-H1..E-H7 campaign)

Product question: within the SHIPPED default budget (`arrange_budget = 500,000,000`
inferences), how much denser can `arrange` now pack each grid size, and where do the
new difficulty cliffs sit?

Method: `crosswordsmith_arrange:arrange_best_layout(Words, G, Budget, _, Reward, Outcome)`
under `call_time`, each in a SEPARATE `timeout 300 swipl --stack-limit=4g` process.
Fixtures are planted-witness meshes from `benchmarks/gen_mesh_fixture.py G N K 3 4 SEED`
(satisfiable + reachable by construction), matching each size's dense-rung alphabet
(9x9 K=5, 15x15 K=6, 21x21 K=6). Budget 5e8 = shipped default; 2e9 re-probe on the
seed-11 first-failures to separate budget-bound from search-bound.
Inference counts are deterministic and machine-independent (only wall varies).

## Sweep @ 5e8 (shipped default)

| size | N | seed | outcome | reward | inferences | wall (s) |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| 9x9  | 17 | 11 | placed      | 182 | 46,707,515 | 2.54 |
| 9x9  | 18 | 11 | placed      | 191 | 341,632,830 | 18.39 |
| 9x9  | 19 | 11 | not_proven  |  -1 | 500,000,218 | 27.66 |
| 9x9  | 20 | 11 | not_proven  |  -1 | 500,000,217 | 29.06 |
| 15x15 | 37 | 11 | placed      | 402 | 23,031,188 | 1.21 |
| 15x15 | 38 | 11 | placed      | 424 | 5,968,765  | 0.33 |
| 15x15 | 40 | 11 | placed      | 443 | 13,467,343 | 0.74 |
| 15x15 | 42 | 11 | placed      | 464 | 4,225,447  | 0.24 |
| 15x15 | 44 | 11 | not_proven  |  -1 | 500,000,243 | 26.68 |
| 21x21 | 82 | 11 | placed      | 916 | 10,793,977 | 0.62 |
| 21x21 | 84 | 11 | placed      | 926 | 100,316,272 | 5.21 |
| 21x21 | 88 | 11 | not_proven  |  -1 | 500,000,287 | 25.93 |
| 21x21 | 92 | 11 | not_proven  |  -1 | 500,000,288 | 27.13 |

## Multi-seed cliff characterization (I4 discipline: 3 seeds at last-placed / first-failed N)

| size | N | seed | outcome | reward | inferences | note |
| --- | ---: | ---: | --- | ---: | ---: | --- |
| 9x9  | 18 | 11 | placed     | 191 | 341,632,830 | near ceiling |
| 9x9  | 18 | 12 | not_proven |  -1 | 500,000,213 | 18w FAILS this seed |
| 9x9  | 18 | 13 | placed     | 189 | 2,022,986   | easy |
| 9x9  | 19 | 11 | not_proven |  -1 | 500,000,218 | |
| 9x9  | 19 | 12 | not_proven |  -1 | 500,000,218 | gen made only 18 words (grid near saturation) |
| 9x9  | 19 | 13 | placed     | 206 | 500,000,513 | MARGINAL — placed right at the ceiling |
| 15x15 | 42 | 11 | placed     | 464 | 4,225,447   | easy |
| 15x15 | 42 | 12 | placed     | 479 | 123,455,825 | |
| 15x15 | 42 | 13 | placed     | 467 | 4,071,273   | easy — 42w robust across seeds |
| 15x15 | 44 | 11 | not_proven |  -1 | 500,000,243 | |
| 15x15 | 44 | 12 | not_proven |  -1 | 500,000,243 | |
| 15x15 | 44 | 13 | placed     | 483 | 95,976,414  | only 1/3 seeds places |
| 21x21 | 84 | 11 | placed     | 926 | 100,316,272 | |
| 21x21 | 84 | 12 | placed     | 930 | 500,001,004 | MARGINAL |
| 21x21 | 84 | 13 | not_proven |  -1 | 500,000,287 | 84w FAILS this seed |
| 21x21 | 88 | 11 | not_proven |  -1 | 500,000,287 | |
| 21x21 | 88 | 12 | placed     | 977 | 500,001,037 | MARGINAL |
| 21x21 | 88 | 13 | placed     | 974 | 500,001,030 | MARGINAL |

## 2e9 re-probe of the seed-11 first-failures (budget-bound vs search-bound?)

| size | N | seed | budget | outcome | inferences | wall (s) |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| 9x9  | 19 | 11 | 2e9 | not_proven | 2,000,000,220 | 108.9 |
| 15x15 | 44 | 11 | 2e9 | not_proven | 2,000,000,239 | 106.4 |
| 21x21 | 88 | 11 | 2e9 | not_proven | 2,000,000,285 | 104.1 |

All three stay unproven at 4x the default budget -> SEARCH-bound (completion is far past
5e8, not just past it), even though neighbouring seeds at the same/higher N place
marginally near 5e8. Classic heavy-tailed backtracking runtime: at a fixed near-cliff N
some seeds finish just under the ceiling and others don't finish at 4x it.

## Envelope (densest N completing under the 5e8 default) vs pre-campaign reference

| size | pre-campaign | post-campaign robust envelope | first fully-noisy N | first all-seed cliff |
| --- | --- | --- | --- | --- |
| 9x9  | ~16w (cliffed)            | 17w (46.7M, robust)  | 18w (1/3 fails; 341M when it places) | 19w (search-bound at 2e9) |
| 15x15 | 36w (182.6M pre; 43.9M now) | 42w (all 3 seeds place) | 44w (2/3 fail) | ~44-46w |
| 21x21 | ~80w (>=82 hard-for-all)  | 82w (10.8M, robust)  | 84w (1/3 fails) | 88w+ (2/3 marginal, s11 search-bound) |

## Recommended new ladder rungs (guard the new envelope; committed here, NOT wired into workloads.pl)

- `fixtures/ladder_09x09_17w.pl` (gen: `9 17 5 3 4 11`) — robust 46.7M / 2.5s. Strongest
  new guard: ~56x the inferences of the current 16w hard rung, pushes the 9x9 envelope
  one word past the pre-campaign ~16w cliff, completes deterministically well under 5e8.
- `fixtures/ladder_15x15_40w.pl` (gen: `15 40 6 3 4 11`) — 13.5M / 0.74s. Density sentinel
  4 words past the current 36w rung; robust, cheap, deterministic.
- `fixtures/ladder_21x21_82w.pl` (gen: `21 82 6 3 4 11`) — 10.8M / 0.62s. Density sentinel
  2 words past the current 80w rung and past the pre-campaign 82w "hard-for-all" line.

NOT recommended as rungs (unsuitable for a ratchet): every near-cliff instance that places
only MARGINALLY at ~5e8 (9x9 19w s13; 21x21 84w s12, 88w s12/s13) — a tiny solver change
would flip them to not_proven and break the ratchet. Rungs must complete with headroom.
