# Plan: the next `arrange` performance campaign

Status: IN PROGRESS (2026-07-17). The adversarial review's amendments are
folded in below: counter-free cutoff authority, a finite paired restart design
with CPU caps, pinned greedy construction subjects, exact raw-pool transpose
semantics, direct-bucket correctness locks, serial wall/profile adjudication,
and explicit topology limits.

## Execution update (2026-07-17)

Phase 0 passed. A-G1 legality-before-score and A-G2 transpose-partner synthesis
were accepted and ratcheted; the accepted product/ratchet base is `4771b97`.
Cumulatively, dense sweep inferences changed by 15x15/32w
`6,126,445 -> 1,241,744` (-79.74%) and 21x21/80w
`28,942,375 -> 4,910,962` (-83.03%). Phase 2 P-D0/P-C0/P-R0 is next, and all
three probes must share one pinned post-doc-sync HEAD.

The authoritative record is `docs/experiments.md`, with detailed evidence in
`benchmarks/results/2026-07-17-a-g1-legality-before-score-premise.md`,
`benchmarks/results/2026-07-17-a-g2-transpose-premise.md`, and
`benchmarks/results/2026-07-17-a-g2-transpose-product.md`.

This campaign investigates the leads in
`docs/research/arrange-performance-leads-2026-07.md` using the measured,
byte-identity-gated, ratcheted process from the first arrange campaign. It is a
new campaign because the targets now span three structurally different
subjects: strict deterministic search, greedy best-effort/candidates, and an
output-changing rescue policy for near-cliff instances.

The authoritative experiment ledger remains `docs/experiments.md`. Every probe,
acceptance, rejection, and correction is appended there when adjudicated.

## Goal

Find and ship output-preserving latency reductions for `arrange`, determine
whether seeded restarts improve strict completion probability under the same
total inference budget, and run one bounded alternative-model prototype. Do not
trade small-rung performance or deterministic output for an unmeasured large-
grid benefit.

The campaign has four tracks:

| Track | Product question | Primary metric | Output policy |
|---|---|---|---|
| G: greedy | Can best-effort/candidates do less repeated work? | greedy inferences, semantic work, wall | byte-identical |
| D/C: strict core | Can placement counting or duplicate failure work shrink? | strict search inferences | byte-identical |
| R: rescue | Can different seeded trees place more cliff instances within 500M? | placement probability, restricted mean inferences | deliberate policy change |
| T: topology | Does a topology-first model expose a better search tree? | nodes/inferences to known feasibility/reward | prototype only |

Track R is not allowed to masquerade as a normal speed optimization. Its goal is
completion probability and tail latency under a fixed total budget. Track T is
not a shipping rewrite unless it independently clears its bounded gates.

## Campaign-level stop condition

Pause when every ranked lead has reached one of:

- an accepted and ratcheted win;
- a measured mechanism-level rejection;
- a failed premise probe;
- a named product-evidence revisit trigger; or
- the bounded prototype's pre-declared stop condition.

Do not continue merely to produce a positive result. A clean null result closes
an avenue successfully.

## Measurement audit

The existing strict-search substrate is healthy. On 2026-07-16,
`make bench-check` reproduced all seven core rungs exactly at `+0.00%` under the
pinned SWI-Prolog 10.1.10. The current ladder covers 9/13/15/21 grids, synthetic
cost steps, real-word anchors, hard envelope guards, and a budget-saturating
latency probe. `bench_core.pl` supplies warmup, repeated samples, deterministic
inference counts, wall/cpu, and process RSS. `baseline.json`,
`check_baseline.pl`, and `history.jsonl` provide a strict 0.5% ratchet and trend
ledger.

### What is already sufficient

| Need | Existing machinery | Verdict |
|---|---|---|
| Strict deterministic metric | `arrange_best_layout/6` search sampler | ready |
| Stable cost ladder | `workloads.pl` core + heavy rungs | ready |
| Regression ratchet | `baseline.json` + `check_baseline.pl` | ready |
| End-to-end latency/RSS | command sampler under GNU time | ready, informational |
| Small output identity | four arrange CLI goldens | ready but too narrow |
| PRNG reproducibility | SplitMix64 vectors + seeded arrange tests | ready |
| Full-tree correctness spot checks | plunit solution-count comparisons | ready for escalation |
| Campaign orchestration | worktree brief/checklist + `campaign-status.sh` | ready after worktree directory creation |
| Statistical profiler | SWI `library(prolog_profile)`, proven by fill probes | reusable |

### Blocking gaps

| Gap | Why it blocks | Phase-0 action |
|---|---|---|
| No arrange ladder identity manifest | Goldens do not cover hard/real strict rungs | add strict identity digests |
| No greedy benchmark or baseline | A-G1/A-G2 cannot be ratcheted | add a separate greedy ladder and ratchet |
| No raw greedy-pool identity | Selected output can hide changed tie/pool order | digest the ordered raw pool and selected layouts |
| No reusable arrange probe rig | `probe_backtrack.pl` calls removed `probe_*` predicates and is currently non-runnable | replace it with exact-replay benchmark-only probes |
| No per-corner/censored sampler | Whole-operation sampler cannot characterize restart distributions | add a single-corner budget sampler and trace schema |
| No fixed near-cliff corpus | P2/P3 regenerated instances ad hoc | byte-freeze campaign-only cliff fixtures |
| No arrange `--promote` | acceptance reruns expensive ladders and can measure a different run | reuse fill's flow and add tested read-back verification |
| No one-command arrange adjudication | tests, identity, and ratchet are separate and easy to omit | add `bench-arrange-verify` and `bench-arrange-promote` |
| Coarse command timer for tiny greedy runs | GNU `%e` is 10ms-granularity | use heavy in-process rungs and semantic counters; command wall remains secondary |
| No search-memory metric | cache/support experiments can spend browser memory silently | fresh-process max RSS plus explicit table/state counts |
| Stale strict-corner documentation | probes could reproduce a retired four-corner model | correct strict=two/shared-budget docs in Phase 0 |

No candidate optimization starts until the Phase-0 acceptance gate passes.

## Constraints

- Default deterministic output is byte-identical for Tracks G/D/C.
- The inference ratchets retain a strict relative 0.5% tolerance and zero
  regressions. No input-size switches or absolute-floor exceptions.
- Candidate/proof enumeration order is load-bearing.
- Non-letter-sharing strict counts may remain stale overestimates. A purported
  optimization must not silently expose exact counts and change MRV ordering.
- One operation-wide inference budget covers every corner/restart attempt.
- Budget interruption means `not_proven`, never `infeasible`.
- Native/WASM seeded streams remain reproducible. Attempt stream derivation is
  versioned and domain-separated if Track R proceeds.
- New probe instrumentation stays in benchmark code or exact-replay twins. Do
  not add a permanent conditional branch to the product hot path merely to
  observe it.
- Search tables reset once per top-level request. Probe and benchmark counts
  must remain history-independent.

## Phase 0: measurement substrate (blocking)

### P0.1 Strict identity and acceptance workflow

Add `benchmarks/check_arrange_identity.sh` and a committed
`benchmarks/arrange_identity.sha256`. The manifest is generated from
`benchmarks/workloads.pl`, with core and heavy selection matching the ratchet.
It hashes the diagnostics-bearing emitted layout produced from the exact
budget-explicit benchmark subject, not an accidentally different search seam.

The identity row includes outcome and reward. A non-placed row cannot hash an
empty success. The latency-only `benchmark_16_dense` row hashes the incumbent
returned under its manifest's shipped 500M budget. It never substitutes a 2B
result for the 500M product result.

Reuse fill's one-measurement `--promote` flow in `check_baseline.pl`: run once,
reject on any regression, then record that exact result and append history.
Unlike the current fill implementation, arrange promotion must also implement
and test the read-back completeness check that fill's comments promise but its
writer does not perform. Add:

```text
make bench-arrange-verify [BENCH_ARGS=--heavy]
make bench-arrange-promote [BENCH_ARGS=--heavy]
```

`verify` runs full tests, arrange identity, and the strict ratchet. `promote`
must read the baseline back after writing and assert that every measured rung is
present with the measured value. Test core-only promotion retaining heavy rows,
new-rung insertion, and latency-only rows.

The strict identity runner is a fresh-process white-box harness. It clears seed
and check-target state, calls budget-explicit `arrange_best_layout/6`, ends the
timed region, and only then serializes through the same diagnostics layout
builder. Its canonical row contains fixture, budget, framing, outcome, reward,
and either the exact diagnostics-bearing JSON bytes or an explicit non-placement
tag. Serialization and hashing never enter the measured search goal.

### P0.2 Greedy product benchmark

Create a separate subject and files rather than adding a mode column to the
strict product bench:

```text
benchmarks/greedy_workloads.pl
benchmarks/run_arrange_greedy.pl
benchmarks/check_greedy_baseline.pl
benchmarks/greedy_baseline.json
benchmarks/greedy_history.jsonl
benchmarks/check_greedy_identity.sh
benchmarks/greedy_identity.sha256
```

Measure four layers:

| Layer | Subject |
|---|---|
| construction | one manifest-pinned seed answer and corner, with expected outcome |
| sweep | complete seed x four-corner greedy pool |
| postprocess | ranking, placement assocs, diversity, numbering |
| command | best-effort/candidates CLI wall and RSS |

The primary portable metric is sweep inferences. Construction count and the
following semantic counters are reported and identity-checked where relevant:
generated crossing descriptors, legality probes, scored candidates, score cell
visits, greedy steps, and completed constructions. Wall is a confirmation
metric because X3 showed a 17--18% wall win with only a 1--2% inference change.
Semantic counters come from a benchmark-only exact-replay twin, not permanent
branches in the product path and not the gated sweep-inference sampler.

Initial greedy ladder:

| Rung | Mode | Purpose |
|---|---|---|
| bundled fixed 17 | strict candidates K=3 | real words, 20 constructions, golden overlap |
| bundled max 11 | best-effort | dropped-word path |
| benchmark_08 fixed 13 | strict candidates K=5 | short return/diversity filtering |
| real_13x13_12w | best-effort | real-word moderate anchor |
| real_15x15_18w | best-effort | real-word size anchor |
| ladder_15x15_32w | best-effort | dense mid/heavy scorer load |
| ladder_21x21_80w | best-effort | high-fan-out heavy scorer load |

The identity manifest includes the ordered raw pool
`(corner,seed,score,placed signature,dropped signature)`, selected outputs,
candidate count, rewards, dropped order, translation-normalized candidate
assocs, and pairwise distances. CLI bytes alone are insufficient.

Characterize cold/warm inference stability and wall spread before recording the
baseline. Require two independent reproductions at `+0.00%` for deterministic
metrics.

Each construction row names `SeedAnswer`, `Corner`, and `Expected` rather than
implicitly choosing the first seed. This is load-bearing on bundled max-11:
`NARRATIVE FALLACY`, `GNOSTIC GOSPELS`, and `ETERNAL RETURN` do not fit from any
corner, while `OMEGA POINT` and `FLOW` do. Use `OMEGA POINT` for the completing
max-11 construction row and add a separate expected-failure row only if failed
seed setup is intentionally measured. `bench_core` remains success-only, so a
reified construction outcome is the sample; predicate failure is never a data
point.

### P0.3 Counter-free authority and instrumented probe utilities

Replace the stale `probe_backtrack.pl` approach with two deliberately separate
benchmark-only paths under `benchmarks/probe_arrange/`.

The **counter-free authority** calls unchanged product predicates. It alone
records cutoff outcome, success inference, censoring, reward, and layout under a
500M inference budget. No instrumentation executes inside that inference limit.

The **instrumented mechanism rig** uses exact-replay local twins for ratios,
event counts, depth/churn, duplicate work, and support transitions. It runs on
completing controls or to a matched decision/node cap outside the product
inference cutoff. It may never claim that its own 500M outcome matches product:
instrumentation consumes inferences and therefore explores less search under the
same limit.

Shared utilities provide:

- fixture loading and exact word-count assertion;
- single-corner setup and budget outcome (`ok`, `budget`, `exhausted`);
- layout/reward/signature comparison against the unchanged product path;
- per-row JSON output and stderr heartbeat;
- inference, wall, process RSS, and explicit state-size counters;
- statistical profiling separate from decision counters;
- deterministic seed installation/cleanup;
- interruption-safe reporting that never labels a cutoff infeasible.

Each instrumented twin must first prove exact replay on unbounded/completing
easy, hard, seeded, and both-corner controls: outcome, layout signature, reward,
and node/decision count where available. Instrumentation overhead is reported
separately as wall overhead and inference overhead. The 15% ceiling applies only
to wall overhead on completing controls and to whether ratios are trustworthy;
it does not certify cutoff equivalence. If overhead exceeds it, reduce the
instrumentation or use sampled counters before allowing the ratios to steer the
campaign.

Do not revive dynamic `probe_*` predicates in `core.pl`. The old
`benchmarks/probe_backtrack.pl` is historical and currently calls predicates no
longer present in the product module; either retire it with a pointer to its
recorded result or make its historical status explicit.

### P0.4 Fixed cliff corpus and restart trace schema

Byte-freeze campaign-only fixtures generated with the P2/P3 recipe, asserting
the generated word count:

| Grid | Words | Alphabet | Fixture seeds | Role |
|---:|---:|---:|---|---|
| 9 | 18 | 5 | 11,12,13 | noisy mini cliff |
| 15 | 44 | 6 | 11,12,13 | daily-size cliff |
| 21 | 84 | 6 | 11,12,13 | Sunday noisy edge |
| 21 | 88 | 6 | 11,12,13 | deeper cliff |

These fixtures are experiment subjects, not ratchet rungs. Marginal instances
must never enter `workloads.pl`.

Define one JSONL attempt schema before collecting data:

```text
rig(authority|instrumented), limit_kind(inferences|nodes|decisions|none),
operation_id, attempt_index, fixture, fixture_seed, search_seed, corner, arm,
cutoff, outcome, success_inferences, censored, max_depth, places, unplaces,
wipeouts, reward, layout_signature, swi_version, commit
```

Fields unavailable to a rig are JSON `null`, never fabricated zeroes. Authority
and instrumented rows are never pooled for outcome or cutoff analysis.

One outer wall timeout protects the batch; each row emits a heartbeat so a long
healthy run is not killed as silent. Search cutoffs are inference counts, never
wall time.

Freeze a master-seed manifest as the first 64 SplitMix64 outputs from one named
constant. The first 16 are pilot seeds, the next 16 tuning seeds, and the final
32 held-out seeds. Policies are paired on the same master seeds. This is a
finite reproducible benchmark set, not a claim of IID sampling from all 64-bit
seeds.

### Phase-0 acceptance gate

- Existing strict core and heavy baselines reproduce exactly.
- Strict identity passes every selected rung twice.
- Greedy baseline and raw-pool identity reproduce exactly twice.
- The counter-free authority matches product cutoff outcomes, while the
  instrumented twin matches completing-control layouts/rewards/decisions and
  reports calibrated wall and inference overhead separately.
- Cliff fixtures regenerate byte-identically and contain the requested word
  count.
- `make test` passes; no existing golden changes.
- Record/promote read-back checks prove every measured rung persisted.
- Strict documentation and PlDoc consistently say two non-transpose searched
  representatives under one shared operation budget; greedy remains four pool
  entries. Historical result documents stay historical but receive no new
  misleading commentary.

## Phase 1: cheap transparent greedy experiments

Run serially because A-G2 must compose on the accepted A-G1 baseline.

### A-G1: legality before score

Hypothesis: the current greedy path fully scores the large majority of crossing
candidates that `check_word_fits/5` later rejects. Move the pure legality probe
before `placement_key/8`.

Pre-registered expectation: a large reduction in scored candidates and material
wall improvement on 15x15/32w and 21x21/80w; smaller inference movement is
acceptable if semantic counters and wall agree. Null/negative is legitimate.

Acceptance:

- strict and greedy identities unchanged;
- zero strict or greedy ratchet regressions;
- scored-candidate count drops substantially on both dense targets;
- wall improvement exceeds characterized spread on at least one dense target.

### A-G2: derive transpose partners

Hypothesis: greedy transpose-pair searches are exact transposes for every seed,
including dropped order, so search two representatives and synthesize all four
old pool entries in their original order.

Pre-registered expectation: construction count falls exactly 50% and sweep wall
approaches a 2x win. Because the expected win is large, escalate verification:
raw-pool identity on the whole greedy ladder, per-seed literal transpose checks,
all CLI identities, full tests, and a fresh differential probe over additional
fixtures.

Any raw pool, tie order, dropped order, reward, or selected candidate difference
rejects the transparent version.

The implementation contract is exact. Current pool generation is corner-major
then seed-major. Derived partners occupy the omitted direct search's old slots;
failed direct constructions imply symmetric omitted entries. Every synthetic
`pw/8` is rebuilt with a fresh clue-number variable -- it must not alias the
source record. Dropped entries preserve input terms and order. A benchmark-only
tagged direct replay records corner/seed/dropped data before the production pool
currently discards them, and direct-vs-derived equality is checked before the
ordinary sort/diversity stages. Memo work will not halve exactly because the
surviving direct searches still populate pure pair tables; the hard expectation
is 50% fewer attempted searches, not exactly 50% wall.

## Phase 2: premise probes (may run in parallel)

These probes share no candidate code and may collect deterministic inference
and counter data concurrently in isolated worktrees after A-G1/A-G2 are
adjudicated. The orchestrator records that accepted commit as the Phase-2 base
SHA, and every probe worktree starts from exactly that SHA even if `main` later
moves. Wall,
profile, RSS, and instrumentation-overhead runs are serialized on a quiescent
host; worktree isolation does not isolate host load. They return data and draft
ledger entries; they do not merge product changes.

### P-D0: support delta and proof multiplicity

For every letter-sharing count refresh, measure previous/exact bucket, residue
survival, legal proofs against only the newest placed word, proof-vs-geometry
multiplicity, candidates examined to saturation, bucket transitions, and watched
cell dirty rates.

Gate A-D2 only if newest-source logic classifies at least 25% of current full
recounts on both light and dense classes without increasing candidate checks on
the light sentinels. Any proof/geometry divergence kills a geometry-only design.

### P-C0: duplicate failed work

Observe without pruning: same-parent duplicate `(WordID,Start,Dir)` children,
revisited canonical placement sets, revisited proved-dead states, and repeated
subtree nodes. Use exact equality on sampled hash hits; dual fingerprints alone
never prove equality.

Kill global caching if repeated dead states are below 1% of recursive entries on
both 15x15 hard corners or account for less than 5% of their node work. Prefer
parent-local dedup if it captures at least 80% of avoidable repetition.

### P-R0: fixed-instance seeded runtime distributions

Pilot first: the first 16 frozen master seeds per non-transpose corner on exactly
four fixed instances: 9x9/18w seed 12, 15x15/44w seed 11, 21x21/84w seed 12, and
21x21/88w seed 11, plus one easy control per grid size. Run the unchanged full
perturbation seam under the shipped 500M cap. The worst-case saturated pilot is
128 cliff rows, about 64 CPU-minutes at the recorded 30s/row, before controls.
Stop early after eight seeds only if every cliff row is censored and no corner
or seed shows depth/progress differentiation; otherwise finish all 16.

If the pilot shows useful short-run mass or complementary corners, use the next
16 master seeds for policy tuning over the seed-11/12 fixtures, and reserve the
final 32 seeds on seed-13 fixtures for held-out adjudication. Do not expand to a
12-fixture x 64-seed census by default. Track R has an eight CPU-hour campaign
cap: at most 1.5 hours pilot, 2 hours tuning, and 4.5 hours held-out tournament.
Use counter-free standalone traces to screen candidate cutoffs offline, but not
to claim exact operation inference accounting: a live rescue controller may
share pure input memos across attempts and therefore charges later attempts
differently. Live counter-free tuning under the exact proposed memo lifecycle
measures at most two nominated policy configurations. Exceeding a phase
allocation is a product-value decision, not an automatic continuation.

Analyze right-censored data correctly. The pilot must distinguish fixture seed
from search seed and report empirical success curves by cutoff. Generator-seed
variation is not evidence of a fixed-instance runtime distribution.

Gate Track R only if at least one fixed cliff instance has material completion
probability across search streams or the two corners show useful complementarity
under 500M. Perturbation-arm complementarity is not claimed by this probe.

## Phase 3: conditional strict experiments

Run one product candidate at a time against the current ratcheted baseline.

### A-D1: stable IDs and direct trailed buckets

Replace per-node count-assoc rebuild/sort work with stable integer word IDs, a
fixed-arity backtrackable bucket term, and stable 0/1/2 partitions. Preserve the
current visible stale-overcount semantics and input-order ties exactly.

Judge this storage change alone. Zero regressions is mandatory. If direct state
does not clear the light rungs, do not layer residues onto it.

Before performance adjudication, add focused correctness locks for a crafted
non-sharing stale-overcount case, equal-bucket input-order ties, sibling
backtrack restoration, proof multiplicity, and exact full-tree solution counts.
Then require complete ladder identity. Stable IDs/direct storage are a
hypothesis, not presumed semantics preservation.

### A-D2: residue plus newest-source delta

Only if P-D0 and A-D1 pass. Use an old residue plus explicit search against the
newly placed word to classify safe bucket transitions; fall back to the current
full recount whenever proof is insufficient. This is not E-H10's unsound
shrinking-domain assumption.

Expected result and exact classification surface come from P-D0, not from the
literature. One focused follow-up is allowed if the mechanism is sound and a
single measured tax blocks the ratchet; otherwise close the avenue.

### A-C1: parent-local failed-child dedup

Only if P-C0 nominates it. Record a child decision only after its subtree
exhausts normally; skip the same failed child when another crossing proof
generates it under the same parent. Never cache budget-interrupted branches.

Global state caching is a separate experiment admitted only if P-C0 shows
substantial non-parent transpositions and a bounded useful working set below the
campaign's memory ceiling.

Memory gate for any cache/support table: no more than 32 MiB additional peak RSS
on the measured native probes without an explicit browser-memory review, and no
unbounded request-to-request growth. Accepted structural changes also run the
WASM battery.

## Phase 4: restart policy tournament (conditional, output-changing)

Track R begins only after P-R0. It lives on an experiment branch until the user
accepts the product trade: first rescue placement may differ from the current
best-of-two-corners output and reward.

The product contract under test is an **opt-in rescue mode**, provisionally
spelled `--rescue`; it never replaces deterministic default search or the
existing one-shot `--seed N` behavior. Rescue requires a master `--seed N` (or
uses the recoverable seed selected by `--shuffle`). Attempt zero is the prefix of
today's coupled two-corner seeded operation: one mutable PRNG stream, with the
second corner inheriting draws consumed by the first, stopped only by its
scheduled cutoff. Domain-separated attempt/corner seeds begin with later rescue
attempts; they do not retroactively redefine attempt zero. Rescue spends only
the same operation-wide budget across all attempts.
Diagnostics record the master seed and a versioned rescue-policy
identifier. A rescued solution is `placed`; exhausting the outer budget remains
`not_proven`. Exact CLI spelling and spec text are locked at the product-policy
checkpoint before product implementation, but these compatibility semantics are
not negotiable inside the experiment.

Compare under one unchanged 500M operation-wide budget:

- empirically tuned fixed cutoffs;
- paired Luby schedules;
- paired geometric schedules;
- current full perturbation first. Root-only/equal-bucket arms are admitted only
  after a cheap diversity probe shows complementarity worth adding new seams.

For attempts after zero, seeds are derived independently from
`(master seed, policy version, corner, arm, attempt index)`. Attempt zero uses
the coupled legacy stream described above. Immutable input memos may survive;
grid, count cache, branch choices, and PRNG state reset per later attempt.
E-H6-style failure weights do not carry across attempts unless separately
justified.

The control is the current operation-wide `/6` behavior, including its one
mutable seeded stream across the two searched corners. Every candidate policy
is paired with that control on the same frozen master seed. The primary finite-
benchmark estimand is the equally weighted mean of each held-out fixture's
success rate over 32 seeds; RMST censors every failure at 500M. Report a
fixture-cluster bootstrap as descriptive only -- four held-out fixture groups
are too few for it to be the acceptance authority.

At most three candidate policies survive tuning into the held-out tournament:
best fixed cutoff, best Luby/geometric schedule, and at most one portfolio arm.
Including control, the worst-case held-out matrix is 4 fixtures x 32 seeds x 4
policies = 512 operation runs, about 4.3 CPU-hours if all saturate.

Tournament outcome, success inference, censoring, reward, and layout come only
from counter-free authority runs. Attempt-level mechanism traces are collected
in separate representative replays outside the product inference cutoff (or to
a matched decision cap) and are descriptive, never the acceptance authority.

Pre-registered acceptance on held-out fixed fixtures:

- placement probability at 500M improves by at least 15 percentage points over
  seeded no-restart;
- restricted mean inference-to-first-placement is at most 80% of control, with
  improvement on at least three of four held-out fixtures and no fixture losing
  placement probability;
- every easy/control fixture places in every stream;
- median reward drops no more than 2% and fifth-percentile reward no more than
  5% among paired successes;
- before the tournament, establish each fixture's current-control reward floor
  only from successful **operation-wide `/6` control** authority runs; single-
  corner pilot rewards are not comparable and never enter it. Every candidate
  success, including a rescue whose paired control failed, must score at least
  95% of that fixture floor. If a fixture has no operation-wide control success,
  it is not eligible to support a quality-safe shipping claim until a reference
  success or explicit product quality threshold is supplied;
- native/WASM attempt traces, seeds, outcomes, rewards, and layouts agree on
  representative rescued and exhausted cases.

If no policy clears the completion gate, do not add restart machinery. If a
policy clears performance but fails the reward policy, park it with the measured
trade and ask for a product decision.

Restart nogoods are not part of the first tournament. Admit them only after
plain diversified restarts show repeated completed work worth preserving.

## Phase 5: topology-first bounded prototype

Run A-T0 in an isolated research branch. It must not touch the shipping path.
Precompute matching-letter events, search over component merges, and maintain
orientation parity, relative coordinate potentials, component bounding boxes,
and selected-edge connectivity. The search must permit merging two components
that do not contain the anchor; otherwise it has recreated the current rooted
DFS.

This track has a four CPU-hour measurement cap and at most two implementation
iterations. Per-row limits are 10M inferences for toy cases, 500M for bundled
and 9x9/17w, and 2B for 15x15/36w, with a 512 MiB process-memory ceiling. Do not
run the 36-word gate if 9x9/17w fails. Exceeding these limits requires a new
measured premise and explicit continuation decision.

Gates, in order:

1. Match exhaustive legality/existence on several 4--8 word cases.
2. Produce a non-root component merge.
3. Reproduce the bundled and 9x9/17-word incumbent feasibility/reward within a
   bounded run.
4. Show that geometry or a reward threshold prunes before complete layouts.
5. Reproduce 15x15/36-word feasibility/reward before seeking an improvement.

Stop immediately at the first failed gate and record why. Iterated
`reward >= K` feasibility is admitted only if crossing-cardinality bounds remove
events or force crossings before leaves. A deterministic one-worker CP-SAT model
may serve as an optional external oracle; if the happy path fails, document the
gap and move on rather than introducing a product dependency.

## Experiment execution protocol

The orchestrator owns the ledger and adjudication. Every runner receives the
eight-part brief from `.claude/skills/experiment-campaign/brief-template.md`:
current HEAD SHA and self-heal command, context, measured mechanism, soundness
argument to challenge, known traps, pre-registered expectation, validation, and
a decision-ready deliverable contract.

Rules:

- Provision each worktree/branch at the explicit accepted base SHA named in its
  brief. For serial experiments this is normally current accepted `main`; a
  parallel probe wave shares one pinned campaign SHA even if `main` later moves.
- One product-code experiment at a time; serial composition preserves
  attribution. Measurement-only probes may run deterministic counters in
  parallel, but wall/profile/RSS/overhead adjudication runs serially on a
  quiescent host.
- Agents never merge, re-record a baseline, or touch the main checkout.
- Every long runner emits a per-row/trial heartbeat at least every ten minutes.
- Accept/reject using the skill's acceptance checklist.
- Tooling fixes land separately from candidate changes.
- Rejections remain on their branch and receive a main-branch ledger entry.
- An unexpectedly large win triggers stronger identity and full-tree checks
  before acceptance.

Branch naming:

```text
campaign/arrange-p0-*
experiment/a-g1-*
experiment/a-g2-*
probe/a-d0-*
probe/a-c0-*
probe/a-r0-*
experiment/a-d1-*
experiment/a-d2-*
experiment/a-c1-*
experiment/a-r1-*
experiment/a-t0-*
```

## Campaign close-out

After the last accepted transparent change:

- rerun strict and greedy core+heavy ratchets and identities;
- repeat the hard-corner attribution probe to update the cost model;
- rerun the robust density-envelope guards under the shipped budget;
- run `make test-wasm` for accepted structural core changes;
- append the final ledger entry;
- write a research summary covering wins, rejections, corrected premises,
  product impact, and named revisit triggers.

If Track R ships, its close-out separately reports success curves, tail latency,
reward distribution, deterministic seed semantics, and native/WASM parity. Do
not blend completion probability with the transparent inference-win headline.

## Recommended order

1. Red-team this plan against current code and empirically challenge Phase-0
   assumptions.
2. Build and accept Phase 0 as infrastructure-only commits.
3. Run A-G1 then A-G2.
4. Dispatch P-D0, P-C0, and P-R0 from the same accepted HEAD; parallelize only
   deterministic counter collection and serialize their wall/profile passes.
5. Run only the strict candidates nominated by those probes, serially.
6. Run the restart tournament only after P-R0 and an explicit product-policy
   checkpoint.
7. Run A-T0 as a bounded independent prototype.
8. Close out with composed product measurements and durable docs.
