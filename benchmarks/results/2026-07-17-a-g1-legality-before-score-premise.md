# A-G1 legality-before-score premise probe (2026-07-17)

Scope: measurement-only benchmark replay on
`probe/arrange-g1-premise`, based exactly on `3ce0b8d`. No product search code,
baseline, identity, ratchet, or golden was changed.

## Method

`benchmarks/greedy_subjects.pl` already contains an exact replay of the greedy
constructor. Before adding counters, its existing replay was run across every
corner x seed attempt on all seven `benchmarks/greedy_workloads.pl` rungs. Each
attempt compared product and replay setup outcome, placed layout term, ordered
dropped answers, and reward.

The replay was then extended with counters around the existing pure operations,
without changing their score-before-legality order:

- generated: each descriptor yielded by `find_intersecting_word/6`;
- scored: each descriptor for which `placement_key/8` succeeds;
- reaching legality: each subsequent call to `check_word_fits/5`;
- legality success/reject: that call's result;
- score-cell visits: exactly `2 * WLen` for each scored descriptor, accounting
  for `crossing_count/7` and the bbox-extension walk;
- avoidable score-cell visits: the same exact charge only for legality rejects;
- `Start < 1`: generated underflow descriptors, counted before scoring.

Search memos reset once at the top of each seven-rung semantic sweep, not per
construction. Product and replay attempts remain corner-major then seed-major;
all four corners and failed seed setups remain present. Grid mutation and
backtracking behavior are unchanged because only the benchmark replay carries
counters and `check_word_fits/5` remains non-binding.

## Results

The rejection percentage denominator is currently scored candidates, which is
also candidates reaching legality in the shipped score-first order. Avoided
percentage is avoidable score-cell visits divided by current score-cell visits.

| workload | generated | Start<1 | reach legality | successes | rejects | reject % | scored | score-cell visits | avoided visits | avoided % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| bundled_17_candidates | 1,574 | 914 | 660 | 344 | 316 | 47.88% | 660 | 15,028 | 7,968 | 53.02% |
| bundled_11_best_effort | 80 | 74 | 6 | 6 | 0 | 0.00% | 6 | 72 | 0 | 0.00% |
| benchmark_08_candidates | 2,119 | 1,491 | 628 | 348 | 280 | 44.59% | 628 | 16,328 | 7,280 | 44.59% |
| real_13x13_12w_best_effort | 8,462 | 2,862 | 5,600 | 928 | 4,672 | 83.43% | 5,600 | 80,280 | 67,788 | 84.44% |
| real_15x15_18w_best_effort | 20,070 | 4,474 | 15,596 | 1,442 | 14,154 | 90.75% | 15,596 | 200,460 | 183,076 | 91.33% |
| ladder_15x15_32w_best_effort | 61,701 | 6,621 | 55,080 | 3,278 | 51,802 | 94.05% | 55,080 | 412,496 | 388,884 | 94.28% |
| ladder_21x21_80w_best_effort | 305,987 | 25,601 | 280,386 | 7,074 | 273,312 | 97.48% | 280,386 | 2,050,244 | 2,000,112 | 97.55% |

The dense-target premise gate is **PASS**: both rejection rates are far above
25% (94.05% and 97.48%). This gates in a product experiment; it is not itself a
product speed or wall-time result.

## Soundness

The product comment that score-first is required to avoid an `arg/3` throw on
`Start < 1` is incorrect. `check_word_fits/5` starts with `Start >= 1`, before
`check_prev_cell/4` or any grid read, so legality-first safely fails on every
underflow. Underflows are substantial (6,621 and 25,601 on the dense targets)
but currently incur zero score-cell visits because `placement_key/8` has the
same leading guard in `crossing_count/7`.

For non-underflow descriptors, both operations are pure. A legality success
guarantees `Start >= 1` and a valid in-grid run, after which `placement_key/8`
is deterministic. Moving legality first therefore filters the same descriptor
stream before scoring, preserves the order of all surviving keyed candidates,
and preserves first-tie-wins. It neither binds the mutable grid nor changes
memo lifecycle. No flaw was found in the legality-first soundness argument once
the stale score-first safety comment is removed or corrected.

## Verification

- Pre-counter exact replay: all seven rungs passed every direct corner x seed
  comparison (20, 20, 20, 20, 16, 8, and 4 attempts by rung).
- Measured exact replay: the same all-seven sweep passed after counters were
  added.
- `benchmarks/check_greedy_identity.sh`: PASS,
  `ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29`.
- Focused `greedy_benchmark` plunit suite: 7/7 passed, including easy/heavy
  replay and counter-partition checks.
- No baseline, identity, ratchet, or golden was recorded or promoted.

## Draft ledger entry

### A-G1 premise - legality before score - MEASURED, PRODUCT EXPERIMENT GATED IN

- **Method/soundness:** benchmark-only exact replay at base `3ce0b8d`; product
  and replay matched setup outcomes, layouts, rewards, and dropped order for
  every corner x seed attempt on all seven greedy rungs before and after counter
  instrumentation. Counters preserve score-before-legality order and one memo
  reset per top-level sweep.
- **Result:** legality rejected 51,802/55,080 scored candidates (94.05%) on
  15x15/32w and 273,312/280,386 (97.48%) on 21x21/80w. Legality-first would
  avoid exactly 388,884/412,496 (94.28%) and 2,000,112/2,050,244 (97.55%)
  current score-cell visits. Both targets clear the pre-registered 25% gate.
- **Correction:** `find_intersecting_word/6` produced 6,621 and 25,601
  `Start < 1` descriptors on the dense targets, but `check_word_fits/5` guards
  `Start >= 1` before any `arg/3`; the product comment claiming score-first is
  needed for safety is stale. Legality-first preserves surviving candidate
  order and first-tie-wins because both tests are pure.
- **Verification/verdict:** full greedy identity remained
  `ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29` and
  focused tests passed 7/7. **Gate A-G1 into a separate product experiment;**
  this premise probe does not claim wall/inference acceptance. No product,
  baseline, identity, ratchet, or golden change was made.
