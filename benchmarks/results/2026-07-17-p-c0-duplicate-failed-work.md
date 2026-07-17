# P-C0: duplicate failed work probe

Date: 2026-07-17
Base: `1bccf47917fea074a404bf467e51865ce988b8b8`
Branch: `probe/a-c0`
Disposition: measurement-only; no product change, pruning, promotion, baseline,
identity, golden, history, or ledger edit.

## Question and controls

The shadow tests whether strict DFS repeats legal children within one parent or
reaches the same absolute failed state through different placement orders. The
two pre-registered 15x15 hard controls are `ladder_15x15_34w/topright` and
`ladder_15x15_36w/topright`: P1 identified `topright` as the thrashing corner on
both fixtures (3,072 and 3,380 unplacements in its broader event probe), whereas
the large 21x21 rung did not backtrack. They both complete normally, so no
instrumented inference cap or censored row is needed. Light sentinels cover
9x9/8w, 15x15/12w, 21x21/25w, and the real-word bundled fixture.

Two earlier hard runs under a 600-second outer health limit were discarded as
`interrupted`: a non-copying mutable-assoc trial retained branch-local keys whose
bindings were trailed away. The final implementation instead asserts copied,
ground keys into per-run benchmark-module facts, indexes fingerprints as atomic
columns, and cleans the facts after every twin. No result below comes from the
discarded rows.

## Exact observation method

Input order assigns stable integer word IDs. The canonical key is:

```text
state(GridSize, sort(RemainingWordIDs),
      sort(p(WordID, AbsoluteStartCell, Direction)))
```

This key is absolute: no translation normalization is applied. Direction is
explicit, and the full placement set plus grid size determines both letter and
boundary grids. The stale count map is deliberately not canonical state: it can
change decision order, but it cannot change whether an identical physical state
and remaining set is feasible.

Two full-depth `term_hash/4` fingerprints only select a lookup bucket. Every
candidate in that bucket is compared to the complete ground key with `==`; a
focused test injects a false fingerprint bucket and proves it is rejected. All
measured rows had zero fingerprint collisions.

A state enters the dead set only from the failure continuation reached after the
ordinary recursive goal exhausts. Semantic node/decision exceptions bypass that
continuation; a zero-node-cap test records zero dead states. Observation never
skips a branch. Parent-local failure remembers `(ParentInvocation,WordID,Start,
Dir)` only after that child's recursive subtree exhausts. Crossing descriptor
duplicates are counted before legality; legal child duplicates are counted only
after `assign_word/9` succeeds.

## Results

Each row had one warmup and two measured full replays. Outcome, reward, layout
signature, every mechanism counter, the complete decision trace (compared by
exact list equality), and instrumented inference count matched exactly between
measured replays. Runtime fields were not compared or used.

| control | nodes | proof dup / proofs | legal child dup | canonical revisits | proved-dead revisits | repeated dead nodes | dead entries | dead work |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 15x15/34w `topright` | 3,107 | 4,006 / 121,889 | 0 | 2,763 | 2,763 | 2,763 | 88.93% | 88.93% |
| 15x15/36w `topright` | 3,454 | 8,347 / 168,951 | 0 | 3,019 | 3,019 | 3,019 | 87.41% | 87.41% |
| 9x9/8w `topright` | 8 | 0 / 21 | 0 | 0 | 0 | 0 | 0.00% | 0.00% |
| 15x15/12w `topright` | 48 | 2 / 295 | 0 | 16 | 16 | 16 | 33.33% | 33.33% |
| 21x21/25w `topright` | 25 | 2 / 83 | 0 | 0 | 0 | 0 | 0.00% | 0.00% |
| bundled/6w `topleft_across` | 105 | 0 / 170 | 0 | 42 | 42 | 42 | 40.00% | 40.00% |

The hard controls contain only 344 and 435 unique canonical states respectively;
310 and 399 become dead. Every measured repeated-dead subtree contains one
recursive node: equivalent placement orders converge immediately on a state
whose normal selection/counting step wipes out. Thus global state caching has a
large node surface, but the avoided subtrees are shallow. The probe does not
claim that node percentage equals an inference or wall percentage.

## Parent-local capture

| control | failed-child hits after exhaustion | captured repeated nodes | share of global repeated nodes |
|---|---:|---:|---:|
| 15x15/34w `topright` | 0 | 0 | 0.00% |
| 15x15/36w `topright` | 0 | 0 | 0.00% |
| all light sentinels | 0 | 0 | 0.00% |

Proof duplicates exist on synthetic fixtures, but none survives legality as the
same child in the same parent. Proof multiplicity therefore does not nominate
A-C1.

## Gates

- **Global negative-state caching: premise passes.** Proved-dead revisits are
  88.93% and 87.41% of recursive entries, and account for the same percentages
  of non-overlapping repeated subtree nodes. Both exceed the 1% entry and 5%
  node-work kill thresholds. This admits a separate global-cache experiment; it
  does not accept an implementation. That experiment must still clear the 32 MiB
  RSS gate, all light inference ratchets, identity, and WASM checks.
- **A-C1 parent-local failed-child dedup: not nominated.** Capture is 0% on both
  hard controls, below the required 80%, with zero legal duplicate children.
  Implementing A-C1 from this evidence would force a candidate contrary to the
  pre-registration.

## Equivalence and validation

`verify_controls.pl` compared unchanged product authority with `none`, `lean`,
and `full` twins on easy, seeded, 9x9 hard, both strict representatives, both
15x15 hard controls, and a seeded two-corner operation. Outcome, reward, and
absolute layout signature matched. Lean and full decision traces matched exactly,
including all failed decisions; node and decision counts matched. The operation
control preserves strict's two representatives, one mutable seeded stream, one
memo reset, and shared operation semantics.

The final verification set was:

```text
focused probe plunit
benchmarks/probe_arrange/verify_controls.pl
benchmarks/probe_arrange/pc0_duplicate_work.pl
make probe-arrange-check
make test
benchmarks/check_arrange_identity.sh --heavy
make bench-check BENCH_ARGS=--heavy
make bench-greedy-identity
make bench-greedy-check BENCH_ARGS=--heavy
```

All passed. `make test` reported 422 plunit executions, all goldens, and all CLI
contract checks green. Every gated strict and greedy inference metric reproduced
at `+0.00%`; both complete identity manifests matched. `git status` showed no
product, baseline, history, identity-manifest, or golden changes.

No dedicated wall/profile/RSS adjudication was made in this parallel probe wave;
the ratchets' mandatory command wall/RSS fields were informational and ignored.
Instrumented inference counts are overhead evidence only and were not interpreted
as product cost or cutoff equivalence.

## Draft ledger entry

### P-C0 - duplicate failed work - MEASUREMENT (global premise passes; A-C1 killed)

- Observe-only exact replay on the two completing 15x15 `topright` hard controls
  found proved-dead state revisits at 2,763/3,107 (88.93%) and 3,019/3,454
  (87.41%) recursive entries, accounting for the same share of node work. Each
  repeated subtree was a one-node wipeout. Exact absolute keys use stable IDs,
  grid size, remaining IDs, start, and direction; dual-fingerprint hits are
  always full-key verified. Dead states are installed only after normal
  exhaustive failure.
- Same-parent legal duplicate children and post-failure parent-local captures
  were zero on both hard controls and every sentinel, despite 4,006/8,347
  duplicate crossing descriptors. A-C1 misses its 80% gate at 0% and is not
  nominated.
- Global negative-state caching clears the 1%/5% premise gates and may receive a
  separate bounded candidate experiment, subject to light-rung, memory, identity,
  and WASM gates. No pruning or product code was implemented by P-C0.
