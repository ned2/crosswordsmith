# Arrange performance campaign 2026-07b summary

Status: campaign complete; Track R closed without an experiment.

## Goal and method

The campaign sought output-preserving latency reductions in strict and greedy
`arrange`, measured fixed-instance restart potential without changing policy,
and bounded one topology-first prototype. Deterministic inference ratchets and
complete output identities were the acceptance authority; serialized paired
wall measurements were supporting evidence. Cutoff outcomes always came from
unchanged counter-free product authority.

The final accepted product state is A-G1, A-G2, A-D1, and A-D2. A-C1/A-C2 and
A-T0 are measured closures. No rescue/restart policy was shipped.

## Accepted interventions

| ID | intervention | intuition and result |
|---|---|---|
| A-G1 | legality before greedy scoring | 94%-97% of dense scored candidates were illegal. Rejecting first removed nearly all avoidable score-cell walks; dense sweep wall improved 59.95%-65.59%. |
| A-G2 | derive greedy transpose partners | The four visible greedy corners form two exact transpose pairs. Searching two and rebuilding fresh partners halved direct attempts while preserving raw-pool order, clue variables, dropped terms, and output; dense wall was about 1.8x faster. |
| A-D1 | stable IDs and direct trailed buckets | Fixed-arity 0/1/2 count storage removes assoc rebuild/read/sort work while preserving stale non-sharing overcounts and stable ties. All strict rows improved 0.94%-5.63%; paired dense wall improved 6.65%-8.02%. |
| A-D2 | proof-preserving newest-source delta | Previous bucket 0/1 can be classified from the newest source plus one fully validated old proof; bucket 2 always fully recounts. All strict rows improved another 5.51%-62.63%. Paired dense wall was null (ratios 1.0013/1.0024), so no wall speedup is claimed. |

## Cumulative product impact

The strict comparison uses the actual prior ratchet carried into Phase 0, with
the two Phase-0 real anchors added before any candidate. All 14 strict rows are
lower by 10.83%-63.17%, including 15x15/36w
`38,215,934 -> 14,074,891`, 9x9/17w
`38,503,164 -> 19,700,703`, and 21x21/80w
`4,246,246 -> 3,107,968` inferences.

No greedy ratchet predated the campaign, so the sound Phase-0 baseline is the
comparison point. Full greedy sweeps are lower by 35.76%-83.03%; the dense rows
are 15x15/32w `6,126,445 -> 1,241,744` and 21x21/80w
`28,942,375 -> 4,910,962`. Postprocess is unchanged.

Current-direct close-out attribution explains why strict costs do not form one
class. The 15x15/34w representatives take `143/3,107` nodes and 36w takes
`40/3,454`, with high right-corner wipeout and unplacement. Both 21x21/80w
representatives take exactly 80 nodes with no backtracking. A-D2 reduces proof
maintenance in these same trees; exact selected/decision traces are unchanged.

## Rejected and closed mechanisms

- P-D0 found 31 proof/geometry divergences. Geometry dedup would change visible
  capped counts and MRV order, so it is closed.
- P-C0 found no parent-local failed-child captures, closing A-C1. It found many
  global dead-state revisits, but each was a one-node wipeout.
- A-C2's eager and lazy exact global caches cut two hard rows 74%-80% but both
  regressed six strict rows by up to 3.70%/3.44%. Exact key/fingerprint admission
  costs lose even on zero-dead controls; revisit only with measured zero-tax
  lookup evidence.
- A-T0 passed toy legality and non-root merge gates, then exhausted 500M on
  9x9/17w after 5,645 complete relative topologies, all product-invalid. Revisit
  only with a sound nonlocal legality propagator.
- Bucket-2 residue lowering remains unsound without complete older-proof
  materialization. A-D2 deliberately leaves it on full recount.

## Envelope and corrected premises

At the unchanged unseeded operation-wide 500M policy, the robust guards still
place with substantial headroom:

| guard | authority inferences | reward |
|---|---:|---:|
| 9x9/17w | 19,699,727 | 182 |
| 15x15/40w | 4,478,919 | 443 |
| 21x21/82w | 3,763,253 | 916 |

Representative fixed cliffs 9x9/18w seed12, 15x15/44w seed11, and 21x21/88w
seed11 all exhaust 500M as right-censored `not_proven`. This does not contradict
P-R0's seeded mass: P-R0 placed 60/128 standalone seeded corner rows. It confirms
that the cliff is trajectory-sensitive/search-bound and is not moved reliably
by cheaper ordinary counting.

Key corrected premises are now durable:

- Strict product operation is two non-transpose representatives under one memo
  reset and one shared inference budget, not four independent corner budgets.
- Ordinary strict cost is often capped counting, but hard 15x15 corners can be
  shallow high-wipeout trees. Grid size alone does not determine the cost class.
- Candidate proof multiplicity is observable semantics; equal geometries are not
  interchangeable.
- The greedy score-first safety rationale was stale because legality guards
  underflow before grid reads.
- Fixed-instance search-seed variation, not generator-seed variation, is the
  evidence needed for restart policy.
- Large deterministic inference reductions do not imply a wall win when term
  and native work trade off. A-D2 is the measured example.

## Permanent verification and historical reconstruction

```sh
make test
make bench-arrange-verify BENCH_ARGS=--heavy
make bench-exact
make bench-greedy-identity
make bench-greedy-check BENCH_ARGS=--heavy
make bench-greedy-exact
make test-wasm
git diff --check
```

The closeout and mechanism commands were campaign-only and are retired. Exact
source snapshots and disposable-checkout recipes are indexed in
`docs/research/benchmark-probe-historical-reconstruction.md`.

Primary artifacts:

- `benchmarks/results/2026-07-17-arrange-campaign-closeout.md`
- accepted/rejected 2026-07-17 reports under `benchmarks/results/`
- `benchmarks/results/2026-07-17-p-r0-pilot.jsonl`
- `benchmarks/baseline.json`, `benchmarks/greedy_baseline.json`, and their
  histories
- `docs/research/benchmark-probe-historical-reconstruction.md`

## Track R closure

P-R0 clears the premise gate for a same-budget restart tournament but does not
authorize product work. Close-out declined an explicit rescue mode whose first
successful layout and reward may differ, so tuning `[16,32)` and held-out seeds
`[32,64)` remain untouched and no product machinery was built. Track R is closed,
not unfinished engineering. Reopening requires a new product decision accepting
that output trade while retaining the registered master-seed, attempt-zero
compatibility, diagnostics, reward-floor, and one-operation-wide-500M contracts.
