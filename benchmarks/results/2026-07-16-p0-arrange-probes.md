# P0.3/P0.4 arrange probe infrastructure (2026-07-16)

Scope: benchmark/test/docs infrastructure only, based exactly on
`aa4d3c86a660b0aef4514488e04a57a684836499` on
`campaign/arrange-p0-probes`. No product search code, baseline, ratchet rung, or
golden was changed.

## Counter-free authority

`benchmarks/probe_arrange/probe_arrange.pl` provides:

- `authority_operation/6`: unchanged product `arrange_best_layout/6`, including
  strict's two representatives and one shared operation inference budget.
- `authority_corner/7`: one memo reset and one seed installation, then unchanged
  product `construct_one/7` for a standalone corner.
- exception-safe seed cleanup and explicit `ok|budget|exhausted` mapping.

No mechanism counter executes inside either product inference limit. A product
budget hit is `not_proven` and censored, never infeasible.

| comparison | budget | authority | direct product | reward/signature |
| --- | ---: | --- | --- | --- |
| bundled operation, forced cutoff | 1,000 | not_proven/budget | not_proven | match |
| bundled operation, representative shipped cap | 500,000,000 | placed/ok | placed | 60 / match |
| bundled standalone top-left | 500,000,000 | placed/ok | unchanged corner product | 60 / match |
| bundled standalone top-right | 500,000,000 | placed/ok | unchanged corner product | 60 / match |
| bundled seeded operation, seed 7 | 500,000,000 | placed/ok | exact replay operation | 60 / match |

The explicit 500M authority CLI sample completed in 113,055 measured
inferences. This is a completing control, not a strict baseline or identity
gate.

An outer 1ms health timeout on a seeded hard control emitted
`outcome=interrupted`, `cutoff=interrupted`, `censored=true`, null reward and
signature. Unit probes verified seed cleanup on normal result, failure,
exception, and interruption paths.

## Instrumented replay

The local twin preserves product seed order, equal-count bucket order, crossing
proof order/multiplicity, capped stale overcounts, grid and boundary trails, and
one-reset memo lifecycle. It invokes unchanged product legality primitives.
`none` has no events; `lean` counts depth/churn/decisions; `full` adds support
transitions, exact duplicate child/state keys, and explicit state sizes.

All controls matched authority outcome, reward, and full placement signature.
Lean and full decision/node counts were exact matches.

| control | corner | seed | reward | decisions/nodes | max depth | places | unplaces | wipeouts |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| easy bundled | top-left across | none | 60 | 105 / 105 | 6 | 105 | 99 | 84 |
| hard 9x9/16w | top-left across | none | 170 | 200 / 200 | 16 | 200 | 184 | 93 |
| seeded bundled | top-left across | 7 | 60 | 7 / 7 | 6 | 7 | 1 | 1 |
| second representative | top-right | none | 60 | 6 / 6 | 6 | 6 | 0 | 0 |

For the same rows, uninstrumented-twin / lean / full inferences were
106,461/108,512/243,526; 510,872/514,630/621,864;
12,326/12,446/21,257; and 10,248/10,348/17,900. These are replay-overhead
measurements, not product cutoff-equivalent results.

Statistical profiling is separate from counters. The uninstrumented hard
control profile completed in 0.023s with 972 profile nodes; sampled top work was
in product grid/crossing/counting predicates. The profile is descriptive only.

## Overhead calibration

Protocol: serialized on the host, one process per control, two discarded
authority/instrumented warmup pairs, then 21 alternating measured pairs. Ratios
use in-process search wall; paired p10-p90 describes observed ratio spread.

Wall overhead:

| control | mode | authority median | rig median | ratio | overhead | paired ratio p10-p90 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| bundled top-left | lean | 3.685 ms | 3.721 ms | 1.010 | +1.0% | 1.000-1.021 |
| hard 9x9/16w top-left | lean | 19.486 ms | 19.059 ms | 0.978 | -2.2% | 0.973-0.981 |
| bundled top-left | full | 3.630 ms | 8.207 ms | 2.261 | +126.1% | 2.230-2.275 |
| hard 9x9/16w top-left | full | 19.577 ms | 38.678 ms | 1.976 | +97.6% | 1.948-1.990 |

Inference overhead, reported separately:

| control | mode | authority | rig | ratio | overhead |
| --- | --- | ---: | ---: | ---: | ---: |
| bundled top-left | lean | 104,691 | 108,511 | 1.036 | +3.65% |
| hard 9x9/16w top-left | lean | 507,625 | 514,630 | 1.014 | +1.38% |
| bundled top-left | full | 104,691 | 243,525 | 2.326 | +132.61% |
| hard 9x9/16w top-left | full | 507,625 | 621,864 | 1.225 | +22.50% |

Full mode exceeds the 15% wall ceiling and must not steer ratios over broad
runs. Lean mode clears the ceiling on both controls and is the admissible mode;
full counters are for sampled or matched-decision mechanism replays only. The
ceiling says nothing about inference-cutoff equivalence.

Fresh-process informational RSS on the same hard single corner was 11,980 KiB
authority and 12,060 KiB lean; both command walls rounded to 0.06s. This
single-sample process metric is not an overhead confidence estimate.

## Fixed corpus

`check_fixtures.py` regenerated to temporary files, asserted exact generated and
committed counts, and byte-compared all files. Two complete check invocations
passed independently.

| fixture | words | SHA-256 |
| --- | ---: | --- |
| cliff_09x09_18w_seed11.pl | 18 | aad34a8655ece5b301ed607a0a199b54465018ba5c539b3cba63581e12076ae7 |
| cliff_09x09_18w_seed12.pl | 18 | 68e58e6a44ffd1456d3a028f441917404ffe1b7f25a2bee4c288060a91e05992 |
| cliff_09x09_18w_seed13.pl | 18 | 1338012d8238cec18f9808d51e911f3765f1da36d9ceb2b3b9c1b46420f2a9ca |
| cliff_15x15_44w_seed11.pl | 44 | 48db3cdf8dceb9931701dc0f23505b3906675b8c4b3879d04e1a6b00bf60ef9e |
| cliff_15x15_44w_seed12.pl | 44 | 7b4741a7184b616d14075a58d0615f5df260c9da6fd312ee85ac528eb8236f13 |
| cliff_15x15_44w_seed13.pl | 44 | f113ef1e1b8efa34a3985e8dd14475d7d31e4bf3ea441d5681eab220ea7618a0 |
| cliff_21x21_84w_seed11.pl | 84 | 02988f651691274de9e68d9664d567f1cfe02ad9aa09a1572e1527987a7e10ca |
| cliff_21x21_84w_seed12.pl | 84 | 1d8a7f6a158aeef7fb6377d9293f7d91f472d7c59dd4e4502f6dd07693ba215e |
| cliff_21x21_84w_seed13.pl | 84 | ca2c80a5805d274c3866288081fc2815f0f1007faeb3c9abfe5e146a6033204f |
| cliff_21x21_88w_seed11.pl | 88 | 02d4bb428b0f7e6b11bc44e713762c3caedeb7955bd183de01ffef1155a2d951 |
| cliff_21x21_88w_seed12.pl | 88 | 392ee6ec1d75625d042f0686e231bfbb5741237f0d43423b8c4883eb6df0bd20 |
| cliff_21x21_88w_seed13.pl | 88 | cc8e6482fda4701608d4b311ce292fd3da70a7781d25324ae2145890c00c0cc8 |

These are campaign subjects only and are absent from `workloads.pl`.

## Seed manifest and schema

The named constant is `CROSSWOR_ASCII_U64 = 0x43524f5353574f52`.
`search_seeds.json` contains its first 64 product-compatible SplitMix64 outputs
and has SHA-256
`b982c035f82d584416f63d8e2411a6b683af979c3c841a27ee582472fdc27384`.
Partitions are `[0,16)` pilot, `[16,32)` tuning, `[32,64)` held out. Two checks
passed against both independent Python and product `splitmix64/3`.

The JSONL schema includes the reviewed P-R0 fields plus `nodes`, `decisions`,
state-size/duplicate counters, measured wall/CPU/inferences, and optional process
wall/RSS. Unavailable values are JSON null. Synthetic tests and actual emitted
authority/instrumented rows passed validation. A deliberately same-group mixed
rig aggregation was rejected.

## Gate verdict and risks

P0.3/P0.4 scope: **PASS**. Counter-free authority matches direct product;
exact replay matches completing controls; lean wall overhead is below 15%; all
12 fixtures regenerate byte-identically at exact counts; seed and schema checks
pass. This is not a strict/greedy baseline or identity-gate claim.

Risks:

- Full mode is too expensive for unsampled ratios; use lean or matched sampled
  full replays.
- The replay deliberately calls private product ordering helpers. Product
  refactors must rerun all exact-replay controls before traces are trusted.
- Wall/RSS evidence is host-specific; rerun serialized on the adjudication host.
- Cliff fixtures are known-satisfiable generated subjects, not IID samples.

## Draft ledger entry

### P0.3/P0.4 - arrange authority, mechanism rig, and fixed restart corpus - INFRASTRUCTURE

- Added a counter-free unchanged-product authority and a separate benchmark-only
  exact-replay mechanism rig under `benchmarks/probe_arrange/`. Authority owns
  inference-cutoff outcomes; instrumented rows are limited to completing or
  matched node/decision controls and cannot be pooled with authority.
- Replay matched easy, hard, seeded, and both strict representatives in outcome,
  reward, layout signature, and lean/full decisions. Lean wall overhead was
  +1.0% and -2.2% on 21-pair controls; full mode was +126%/+98%, so full is
  sampled-only while lean clears the 15% steering gate.
- Froze 12 exact-count cliff fixtures and 64 SplitMix64 seeds from
  `CROSSWOR_ASCII_U64`, with two independent regeneration checks. Added a null-
  strict JSONL schema and mixed-rig aggregation rejection.
- Retired the stale P1 probe as explicitly historical; no product `probe_*`
  hook was restored. Tests and goldens pass; product Prolog and golden diffs are
  empty. Verdict: P0.3/P0.4 PASS, with full-counter overhead as a documented
  sampled-use limitation.
