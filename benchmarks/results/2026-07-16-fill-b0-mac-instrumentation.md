# B0-I: shipped MAC structural instrumentation (2026-07-16)

Measurement-only result for `docs/plans/fill-perf-program-2026-07.md` B0-I.
No file under `prolog/`, no oracle, and no baseline was changed.

**Budget-accounting correction.** Adjudication found that the first rig revision
placed MAC setup before `call_with_inference_limit/3`, unlike product
`fill_attempt/8`, where only `seed_used/3` is outside and all of
`fill_search/5` is inside. The canonical values below come from the corrected
rig: `setup_core/5` (mask construction, domains, edges, weights), root AC, and
all attempts share the supplied inference budget. Setup wall is recorded from
inside that goal, and interrupted runs retain valid setup/attempt state. The
authority counters did not move; the bounded CWL row did and is replaced here.

## Provenance and method

- Base: `c7ae472a99ce31b20a99bd35b02466ff8f1500fa`, branch
  `experiment/b0-mac-instrument`, host `bunk x86_64`, SWI-Prolog 10.1.10.
- External STW: `/tmp/opencode/stw-lowercase-crossfire.txt`, 315,905 lines,
  SHA-256 `9aa4563eb82f9af45f3a9c4fbc7fe3c8c3e27da6b76cc48f5e09c499bd17e4ee`.
  It is CC BY-NC-SA and no words or generated fills are present here.
- Rig: `benchmarks/fill_quality/probe_mac_b0.pl`, invoked by
  `benchmarks/fill_quality/run_b0_instrument.sh`. Setup and unchanged helpers
  are module-qualified calls into `crosswordsmith_fill`; local twins are
  limited to restart/search/place/propagate/revise/bump. Candidate enumeration
  remains the shipped lazy `mac_candidate/4`; bit clearing remains xor; attempt
  1 is greedy, attempts 2+ use the pinned `0xCC9E2D51` stream; caps are
  `500 x 1.5`; weights persist and age `x 0.99`; root AC precedes attempts.
- Budget boundary: exactly as `fill_attempt/8` lines 1689-1707, `seed_used/3`
  is outside; the complete probe equivalent of `fill_search/5`, beginning with
  mask/setup construction, is inside both the inference and optional wall cap.
- Every accepted fill is rebound to the original shared slot cells through
  `mac_bind_fill/2` before `slots_to_layout/3`. A crossing error therefore
  fails the probe rather than producing data.
- Counter mode records exact queue length at every pop; redundant, fruitful,
  and DWO revisions; deleted popcount; placement DWO; attempts/nodes; and
  per-edge bump/final/peak vectors. Timing mode clocks only whole phase calls,
  never individual revisions. `profile/2` was also run as a self-time
  cross-check; generic arithmetic attribution made the coarse phase wrappers
  the clearer split.

Representative commands:

```sh
benchmarks/fill_quality/run_b0_instrument.sh compare \
  fixtures/fill_grid_04a.json fixtures/dict/enable_50k.txt 1 none 2000000000
benchmarks/fill_quality/run_b0_instrument.sh run blocked13a-stw-30 \
  grids/blocked_13a.json /tmp/opencode/stw-lowercase-crossfire.txt 30 none \
  800000000 0 counters
benchmarks/fill_quality/run_b0_instrument.sh run blocked13a-cwl50 \
  grids/blocked_13a.json dicts/cwl50.dict 50 none 800000000 240 counters
benchmarks/fill_quality/run_b0_instrument.sh run blocked13a-stw-30-timing \
  grids/blocked_13a.json /tmp/opencode/stw-lowercase-crossfire.txt 30 none \
  800000000 0 timing
```

## Soundness

The uninstrumented twin was compared in-process against `fill_attempt/8` on
plain 3x3, seeded split-3, scored 3x3 @50, and ladder `sq04_50k`. Outcomes,
complete `Numbered` layouts, input words, and final-attempt nodes all matched:
`9/9`, `1/1`, `7/7`, and `34/34`; each completed in attempt 1. The same replay
was then run for all five quality masks: product/twin nodes were open4 `8/8`,
open5 `20/20`, mini7 `23/23`, mini9 `28/28`, amer11 `90/90`; every complete
layout and word list matched.

The corrected multi-attempt replay also passed both authority rows against the
product under the same 800M accounting. STW @30 matched outcome, exact
`Numbered`, exact `InputWords`, and final-attempt nodes `919/919` (twin reports
7 attempts); STW @1 matched the same fields and nodes `228/228` (2 attempts).

Counter mode was run twice independently on `sq04_50k` and, after correction,
both authority rows; the corrected bounded CWL row was also repeated. After
deleting only label and wall-time fields, outcome, fill digest, every counter,
every attempt record, and every final/peak weight were exactly equal. Authority
fill digests are @30
`f1e37425d9ac5b6ba9fcbfad46a5a2416d87453224d6019a2c2ffe95d51308ee`
and @1 `82cd1956ced834415a59c3fd910670b158ccaa91bb3c32f496c81def1e887d78`.

Counter overhead on easy ladder rung `sq04_50k`, three fresh-process runs per
mode, search-window medians: off `7.655 ms` (samples 7.623/7.655/7.740),
counters `8.068 ms` (7.813/8.068/8.318), **+5.39%**. This passes the <=15%
gate. Setup is excluded from this ratio because it is identical and dominates
the sub-10ms search window.

## Row table

`<=2`, `=3`, and `>=4` are percentages of queue pops. Revisions include root AC.
Nodes sum all attempts, including the cap-crossing node. `top25` and Gini use
per-edge learned bump counts (the aging-independent learned-excess measure).

| row | outcome | att | nodes | queue pops | <=2% | =3% | >=4% | redundant | fruitful | rev DWO | fruitful% | place DWO | deleted | bumps | top25% | Gini |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| STW @30 | filled | 7 | 11,318 | 1,023,870 | 0.22 | 0.27 | 99.51 | 3,298,900 | 1,279,851 | 10,945 | 27.89 | 0 | 401,483,750 | 10,945 | 83.54 | 0.737 |
| STW @1 | filled | 2 | 729 | 51,822 | 1.04 | 0.89 | 98.08 | 159,840 | 65,311 | 599 | 28.93 | 0 | 21,254,814 | 599 | 91.49 | 0.799 |
| sq04_full | filled | 1 | 8 | 57 | 45.61 | 24.56 | 29.82 | 135 | 48 | 0 | 26.23 | 0 | 8,131 | 0 | 0.00 | 0.000 |
| g11_full | filled | 1 | 45 | 390 | 49.74 | 20.26 | 30.00 | 962 | 317 | 1 | 24.77 | 0 | 103,229 | 1 | 100.00 | 0.988 |
| g11_full_seed | filled | 1 | 41 | 314 | 53.50 | 18.15 | 28.34 | 656 | 270 | 0 | 29.16 | 0 | 44,548 | 0 | 0.00 | 0.000 |
| sq05_full | filled | 1 | 17 | 233 | 25.32 | 23.18 | 51.50 | 558 | 290 | 7 | 33.92 | 0 | 29,968 | 7 | 100.00 | 0.846 |
| g17_full | filled | 1 | 107 | 882 | 46.60 | 19.84 | 33.56 | 2,107 | 733 | 5 | 25.76 | 0 | 177,481 | 5 | 100.00 | 0.982 |
| g21_full | filled | 4 | 2,540 | 5,656 | 57.66 | 21.48 | 20.86 | 12,700 | 4,488 | 27 | 26.07 | 0 | 1,559,743 | 27 | 100.00 | 0.947 |
| g13_full | filled | 1 | 66 | 661 | 39.64 | 26.02 | 34.34 | 1,790 | 626 | 9 | 25.81 | 0 | 224,999 | 9 | 100.00 | 0.950 |
| sq04_50k | filled | 1 | 34 | 343 | 38.19 | 35.28 | 26.53 | 561 | 424 | 17 | 42.32 | 0 | 7,038 | 17 | 76.47 | 0.665 |
| g15_full | filled | 1 | 99 | 854 | 50.47 | 22.01 | 27.52 | 2,061 | 758 | 11 | 26.78 | 0 | 196,890 | 11 | 100.00 | 0.937 |
| g17_50k | filled | 4 | 2,744 | 20,084 | 29.57 | 28.23 | 42.20 | 35,373 | 24,337 | 1,393 | 39.83 | 0 | 1,122,996 | 1,393 | 87.51 | 0.774 |
| g09_full | filled | 2 | 533 | 925 | 68.65 | 16.86 | 14.49 | 1,509 | 677 | 18 | 30.72 | 0 | 138,001 | 18 | 100.00 | 0.935 |
| cwl50 @50 | not_proven | 9 | 27,342 | 1,169,547 | 0.91 | 1.11 | 97.98 | 4,914,095 | 1,552,820 | 26,381 | 23.91 | 0 | 811,392,821 | 26,381 | 93.83 | 0.808 |

The ladder used its committed 2B per-search harness budget to obtain exact
completion; the two STW rows and CWL stretch used the shipped 800M budget.
The corrected stretch charged 2.963s of MAC setup and 229.466s of root/attempt
search to the same 800M inference-limited goal. It hit that budget after
232.4s limited / 236.1s end to end, before the 240s wall cap: it is
`not_proven`, not infeasible. All placement DWOs were zero; all learned bumps
came from propagation revision DWOs.

## Attempt table

`edges+` is the count of edges bumped in that attempt. `max final/peak` shows
the weight distribution's upper endpoint; the rig's JSON output also emits the
complete per-edge bump/final/peak vectors and edge endpoint map.

| row | att | cap | stop | nodes | pops | redundant | fruitful | DWO | deleted | bumps | edges+ | max final | max peak |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| STW30 | 1 | 500 | cap | 501 | 35,412 | 118,317 | 43,338 | 461 | 26,790,582 | 461 | 41 | 46.000 | 46.000 |
| STW30 | 2 | 750 | cap | 751 | 72,670 | 242,526 | 86,418 | 726 | 23,521,301 | 726 | 40 | 72.990 | 72.990 |
| STW30 | 3 | 1,125 | cap | 1,126 | 124,160 | 407,241 | 155,599 | 1,116 | 65,988,822 | 1,116 | 47 | 132.901 | 132.901 |
| STW30 | 4 | 1,688 | cap | 1,689 | 107,040 | 345,541 | 134,624 | 1,641 | 28,438,818 | 1,641 | 37 | 405.153 | 405.153 |
| STW30 | 5 | 2,532 | cap | 2,533 | 143,485 | 447,284 | 190,287 | 2,476 | 57,871,080 | 2,476 | 48 | 401.102 | 401.102 |
| STW30 | 6 | 3,798 | cap | 3,799 | 437,012 | 1,382,658 | 546,496 | 3,688 | 146,145,940 | 3,688 | 52 | 658.077 | 658.077 |
| STW30 | 7 | 5,697 | filled | 919 | 104,023 | 355,007 | 123,067 | 837 | 52,726,780 | 837 | 39 | 651.496 | 651.496 |
| STW1 | 1 | 500 | cap | 501 | 32,248 | 89,301 | 43,619 | 449 | 8,624,217 | 449 | 42 | 46.000 | 46.000 |
| STW1 | 2 | 750 | filled | 228 | 19,513 | 70,229 | 21,681 | 150 | 12,630,264 | 150 | 36 | 46.540 | 46.540 |
| CWL | 1 | 500 | cap | 501 | 31,960 | 108,453 | 45,755 | 449 | 11,245,167 | 449 | 54 | 46.000 | 46.000 |
| CWL | 2 | 750 | cap | 751 | 36,838 | 153,604 | 48,368 | 741 | 30,846,279 | 741 | 35 | 114.950 | 114.950 |
| CWL | 3 | 1,125 | cap | 1,126 | 51,518 | 200,506 | 67,040 | 1,094 | 40,362,746 | 1,094 | 35 | 287.901 | 287.901 |
| CWL | 4 | 1,688 | cap | 1,689 | 91,491 | 387,134 | 116,627 | 1,663 | 76,110,943 | 1,663 | 39 | 461.662 | 461.662 |
| CWL | 5 | 2,532 | cap | 2,533 | 123,536 | 528,164 | 162,659 | 2,438 | 78,148,431 | 2,438 | 40 | 561.171 | 561.171 |
| CWL | 6 | 3,798 | cap | 3,799 | 150,097 | 633,155 | 203,310 | 3,705 | 105,879,800 | 3,705 | 46 | 923.560 | 923.560 |
| CWL | 7 | 5,697 | cap | 5,698 | 247,307 | 1,066,283 | 321,189 | 5,461 | 176,042,801 | 5,461 | 45 | 1,660.261 | 1,660.261 |
| CWL | 8 | 8,546 | cap | 8,547 | 311,601 | 1,305,734 | 424,478 | 8,196 | 198,311,077 | 8,196 | 48 | 1,914.181 | 1,914.181 |
| CWL | 9 | 12,819 | budget | 2,698 | 125,138 | 530,752 | 163,383 | 2,634 | 94,445,455 | 2,634 | 44 | 2,143.221 | 2,143.221 |

## Kill tests and mechanism

- **F1/B1 survives.** Queue `<=2` is 0.22% @30 and 1.04% @1, nowhere near
  90%; queue `>=4` is 99.51% and 98.08%, and many ladder rows exceed 20%.
  Revision ordering has abundant structural surface. Queue lengths routinely
  reach the high 20s/low 30s on authority attempts, so this is not merely a
  root-AC artifact.
- **F4/B2 survives.** Fruitful revisions are 27.89% and 28.93%, far above the
  2% kill threshold. Deleted-popcount totals are very large, while placement
  causes no DWO; deletion-responsible and fully-assigned credit can distinguish
  much more history than the current final-blow-only revision bump.
- **F2/F3 and B2/B3 survive.** Every authority attempt has 150-3,688 learned
  bumps, above 20. Top-quartile learned-excess share is 83.54%/91.49% and Gini
  0.737/0.799, so neither smear condition holds. Aging has concentrated,
  persistent state to act on. F5 probing-init has no direct B0 kill test; it
  also survives because greedy attempt 1 itself generates 461/449 bumps before
  restart, demonstrating a real warm-start signal.

Thus **B1 survives; B2 survives; B3 survives**. This is a structural verdict,
not evidence that a variant wins. The next warranted experiments are the
pre-registered order: F1 queue/edge ordering first; then F2+F4 credit arms with
and without the winning F1 order; then intra-attempt aging and pinned probing
initialization. No variant was built in B0.

## Wall attribution

Coarse timing runs produced the same fill digests as counter runs:

| row | total s | setup s (% total) | select/candidate s | placement s | propagation/support s (% total; % search) | other search s |
|---|---:|---:|---:|---:|---:|---:|
| STW @30 | 146.262 | 6.037 (4.13%) | 0.058 | 0.525 | 139.583 (95.43%; 99.54%) | 0.059 |
| STW @1 | 17.738 | 8.947 (50.44%) | 0.015 | 0.024 | 8.746 (49.31%; 99.49%) | 0.006 |

Search is unequivocally propagation/support dominated at both wide bands, so
subsequent search-latency claims must be described as bignum-path wins. At @1,
however, end-to-end latency is co-dominated by dictionary/grid/MAC setup; a
search-only win cannot remove that half. The current STW snapshot reproduces
the historical timing/quality closely (@30 146.3s vs 140s, @1 17.7s vs 19s),
but @1 now completes in 2 attempts rather than the historical note's 4; the
product CLI independently showed the same 2-attempt result.

## Quality replay

The five generated masks used only `/tmp`; no generated grid/fill was committed.
All are STW @50, matching the existing harness policy, and every twin assignment
matched the product assignment exactly.

| mask | outcome | nodes | n | mean | min | below50 |
|---|---:|---:|---:|---:|---:|---:|
| open4 | filled | 8 | 8 | 50.0 | 50 | 0 |
| open5 | filled | 20 | 10 | 50.0 | 50 | 0 |
| mini7 | filled | 23 | 23 | 50.0 | 50 | 0 |
| mini9 | filled | 28 | 28 | 50.0 | 50 | 0 |
| amer11 | filled | 90 | 44 | 50.0 | 50 | 0 |

## Verification

- `make unit`: PASS, 389/389 instances.
- Probe load safety (`swipl -q -l load.pl -l .../probe_mac_b0.pl -g halt`): PASS.
- `make test`: PASS, 389/389 instances, all goldens and CLI contracts.
- 11-rung CLI identity: PASS, every digest matched; no re-recording.
- 11-rung search/load inference ratchet (`--heavy`): PASS, every search and
  load count exactly **+0.00%**; no re-recording.

## Draft append-only ledger entry

> ### B0-I - shipped MAC instrumentation - MEASURED, all Track-B arms survive
> Measurement-only `probe_mac_dwd`-lineage twin at commit `c7ae472`: exact
> product replay on four easy/seeded/scored rows, both multi-attempt authority
> rows, and all five quality masks;
> duplicate authority runs matched all fills/counters/weights; easy-rung
> overhead +5.39%. STW `blocked_13a` @30/@1 filled in 7/2 attempts. Queue
> `<=2` was 0.22%/1.04%, `>=4` 99.51%/98.08%; fruitful revisions
> 27.89%/28.93%; 10,945/599 learned bumps; top-quartile excess
> 83.54%/91.49%, Gini 0.737/0.799. Therefore F1, F2/F4, and F3/F5 survive
> their B0 gates; build variants in preregistered order. Propagation/support
> consumed 99.54%/99.49% of search wall (bignum-path dominated). All 11 ladder
> rows filled at their 2B harness budgets; CWL @50 stopped `not_proven` on the
> shipped 800M budget before 240s, with MAC setup charged inside that budget.
> Quality replay stayed min/mean 50 with exact
> product assignments. Product tests, identity, and ratchets PASS at +0.00%; no
> baseline/oracle/plan/ledger change.
