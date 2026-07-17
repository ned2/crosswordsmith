# A-G2 greedy transpose premise probe

Date: 2026-07-17
Base: `c584962f499f7118a1850b0254d7ea4f925874f1`
Branch: `probe/arrange-g2-premise`
Verdict: **PASS - gate transparent A-G2 product experiment in**

## Method

`benchmarks/run_arrange_g2_probe.pl` runs a benchmark-only tagged direct replay
of the unchanged four-corner greedy sweep. It resets search memos once, then
records every slot in the existing order: all `topleft_across` seeds, all
`topleft_down` seeds, all `topright` seeds, then all `bottomleft` seeds. Setup
failures and strict-ineligible completed constructions remain explicit records.

For each seed it compares `topleft_across -> topleft_down` and
`topright -> bottomleft`. For a 1-based row-major cell `s` on an `N x N` grid:

```text
r = (s - 1) // N
c = (s - 1) mod N
transpose(s) = c*N + r + 1
transpose(across) = down
transpose(down) = across
```

The comparison independently checks setup outcome, mode eligibility, placed and
dropped counts, reward, placement-list and answer order, all `pw/8` geometry
fields, literal dropped input terms/order, translation-normalized placement
assocs/direction, and the unchanged eligible raw-pool stream. It also applies
the formula twice and requires the original geometry.

## Manifest results

| workload | direct slots | setup failed | completed | eligible | ineligible | pair mismatches |
|---|---:|---:|---:|---:|---:|---:|
| bundled_17_candidates | 20 | 0 | 20 | 8 | 12 | 0 |
| bundled_11_best_effort | 20 | 12 | 8 | 8 | 0 | 0 |
| benchmark_08_candidates | 20 | 0 | 20 | 4 | 16 | 0 |
| real_13x13_12w_best_effort | 20 | 0 | 20 | 20 | 0 | 0 |
| real_15x15_18w_best_effort | 16 | 0 | 16 | 16 | 0 | 0 |
| ladder_15x15_32w_best_effort | 8 | 0 | 8 | 8 | 0 | 0 |
| ladder_21x21_80w_best_effort | 4 | 0 | 4 | 4 | 0 | 0 |
| **total** | **108** | **12** | **96** | **68** | **28** | **0** |

Both transpose pairs checked half the slots on every row (54 paired comparisons
total). Every row's tagged eligible stream matched the product raw pool in exact
corner-major/seed-major order.

## Differential sample

| deterministic set | grid | slots | setup failed | completed with drops | all-fit | mismatches |
|---|---:|---:|---:|---:|---:|---:|
| setup_failure (`ABCDE`, `AXE`, `BEE`) | 3 | 12 | 4 | 8 | 0 | 0 |
| dropped_word (`CAT`, `CAR`, isolated `DOG`) | 7 | 12 | 0 | 12 | 0 | 0 |
| all_fit (`CAT`) | 7 | 4 | 0 | 0 | 4 | 0 |

The generated entries carry metadata. Dropped lists were required to contain
the original input terms as an identity-preserving subsequence, not reconstructed
answers, and direct partners had exactly the same dropped terms/order.

## Implication and soundness

A future implementation can directly attempt only `topleft_across` and
`topright`: attempts by manifest rung become `10,10,10,10,8,4,2`, exactly 54
instead of 108 (50%). The omitted partners must still occupy their historical
raw slot blocks; the 12 failed setup slots imply 6 direct failures plus 6
symmetric omitted failures. This probe does not implement that synthesis.

Every directly compared source/partner `pw/8` clue-number field was an unbound
variable distinct within and across the paired placement lists. The independent
transpose helper creates a fresh variable per transformed record, verifies it
does not alias the source or another transformed record, and verifies two
transposes restore geometry while creating another fresh set. A product
synthesizer must preserve that behavior.

## Verification

Commands:

```sh
swipl -q benchmarks/run_arrange_g2_probe.pl
swipl -q -g "consult('load.pl'),consult('tests/greedy_benchmark.plt'),run_tests(greedy_benchmark),halt"
benchmarks/check_greedy_identity.sh
make test
```

No product source, baseline, history, manifest, identity digest, or golden file
is changed by this probe. Focused greedy tests passed 9/9; `make test` passed
414 assertions, every golden, and all CLI contracts; full greedy identity stayed
`ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29`.

## Draft ledger entry

### A-G2 premise - greedy direct transpose pairs - MEASURED, PRODUCT EXPERIMENT GATED IN

- **Method/soundness:** measurement-only benchmark replay at base `c584962`.
  Enumerated the unchanged four corner-major then seed-major blocks, retaining
  setup failures, strict eligibility, complete `pw/8` records, and original
  dropped terms. Compared direct pairs using an independent row/column transpose
  and checked involution, normalized direction, raw-pool order, and fresh clue
  variables.
- **Result:** all 108 manifest slots (54 pairs) matched across all seven greedy
  rungs: zero setup, geometry, placed-order, dropped-order, reward, eligibility,
  normalized-assoc, or raw-pool mismatches. Attempt counts were exactly
  `20,20,20,20,16,8,4`; 12 setup failures paired symmetrically. Three additional
  generated square sets checked 28 slots spanning setup failure, dropped words,
  and all-fit behavior, also with zero mismatches.
- **Implication/verdict:** direct attempts can fall exactly 108 -> 54 (50%) if a
  product candidate later searches `topleft_across`/`topright` only and inserts
  independently rebuilt transpose partners into the omitted historical slots.
  Dropped terms must remain original terms/order and every rebuilt `pw/8` must
  receive a fresh non-aliased clue variable. **Gate A-G2 product experiment in.**
