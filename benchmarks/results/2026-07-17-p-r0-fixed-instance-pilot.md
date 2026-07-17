# P-R0 fixed-instance seeded runtime pilot (2026-07-17)

## A. Base and scope

- Worktree: `/tmp/opencode/crosswordsmith-ar0`; branch: `probe/a-r0`.
- Measurement base: `1bccf47917fea074a404bf467e51865ce988b8b8` exactly.
- Scope: benchmark runner, authority JSONL, analysis, report, and a benchmark
  schema CLI correction only. Product code, baselines, histories, identities,
  workloads, and goldens are unchanged.
- Question: fixed-instance search-stream variation under the unchanged full
  perturbation seam. Fixture-generator seed was held fixed and is never treated
  as runtime-distribution evidence.

## B. Protocol, provenance, and allocation

The serial launcher ran each registered cliff fixture at `topleft_across` and
`topright`, independently, over frozen master seed indices `[0,16)`. Every row
used unchanged product `construct_one/7`, a standalone request lifecycle with one
memo reset and one seed installation, and the shipped 500,000,000-inference
cutoff. The outer 120-second wall timeout was health protection only. No wall,
profile, or RSS value adjudicates the result.

The three controls used the first pilot search seed on `ladder_09x09_08w`,
`ladder_15x15_12w`, and `ladder_21x21_25w`. All rows have explicit unique
operation IDs, attempt index zero, fixture seed, search seed, seed index, arm,
SWI version 10.1.10, and base commit. Rows were fsynced individually and resume
skips existing operation IDs.

The run consumed 1,887.630 child CPU-seconds (31.461 minutes, 0.524 CPU-hours),
including process startup; authority's in-process search CPU sum was 1,873.566
seconds. This is 35.0% of the 1.5 CPU-hour pilot allocation. The hard stop was
not approached.

## C. Outcomes

Every budget termination below is right-censored `not_proven`. There were no
`infeasible`, outer-interrupted, or other rows. `success_inferences` is reported
only for authority completions; censored rows have no invented runtime.

| fixed fixture | corner | placed / 16 | censored / 16 | success inference range | reward range |
| --- | --- | ---: | ---: | ---: | ---: |
| 9x9/18w, fixture seed 12 | top-left across | 7 | 9 | 34,949,473-429,036,909 | 184-208 |
| 9x9/18w, fixture seed 12 | top-right | 6 | 10 | 386,894-60,497,736 | 184-201 |
| 15x15/44w, fixture seed 11 | top-left across | 6 | 10 | 630,592-176,530,811 | 473-492 |
| 15x15/44w, fixture seed 11 | top-right | 7 | 9 | 687,988-146,072,193 | 466-491 |
| 21x21/84w, fixture seed 12 | top-left across | 10 | 6 | 2,954,384-198,887,266 | 906-951 |
| 21x21/84w, fixture seed 12 | top-right | 11 | 5 | 2,534,511-441,866,743 | 906-945 |
| 21x21/88w, fixture seed 11 | top-left across | 10 | 6 | 2,728,480-427,670,140 | 969-1,000 |
| 21x21/88w, fixture seed 11 | top-right | 3 | 13 | 2,985,776-4,462,537 | 963-998 |

All 60 successful cliff rows had distinct layout signatures within their
fixture/corner group. The controls placed in 13,776, 45,327, and 142,816
authority inferences respectively, with zero censoring.

## D. Empirical success curves and complementarity

The empirical curves count a seed as successful by a cutoff when its counter-free
authority completion inference is at or below that threshold; all other rows
remain censored. Each cell is completions out of 16. These are fixed finite seed
counts, not IID probability estimates.

| fixture/corner | 1M | 10M | 50M | 100M | 250M | 500M |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 9x9/18w top-left | 0 | 0 | 1 | 3 | 6 | 7 |
| 9x9/18w top-right | 1 | 4 | 4 | 6 | 6 | 6 |
| 15x15/44w top-left | 3 | 4 | 4 | 5 | 6 | 6 |
| 15x15/44w top-right | 3 | 4 | 4 | 6 | 7 | 7 |
| 21x21/84w top-left | 0 | 8 | 9 | 9 | 10 | 10 |
| 21x21/84w top-right | 0 | 9 | 9 | 9 | 9 | 11 |
| 21x21/88w top-left | 0 | 5 | 5 | 5 | 8 | 10 |
| 21x21/88w top-right | 0 | 3 | 3 | 3 | 3 | 3 |

Same-seed corner pairing shows useful standalone complementarity:

| fixture | both | top-left only | top-right only | neither | either | gain over better corner |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 9x9/18w seed 12 | 2 | 5 | 4 | 5 | 11 | +4 |
| 15x15/44w seed 11 | 3 | 3 | 4 | 6 | 10 | +3 |
| 21x21/84w seed 12 | 6 | 4 | 5 | 1 | 15 | +4 |
| 21x21/88w seed 11 | 2 | 8 | 1 | 5 | 11 | +1 |

This is corner complementarity under independent standalone authority requests.
It is not perturbation-arm complementarity and does not predict exact current or
future operation-wide accounting: a live operation shares immutable memos and
today's seeded operation carries one mutable stream from the first corner into
the second.

## E. Stop and gate decisions

The eight-seed early stop did not apply: successful rows and corner
differentiation appeared in the first seed block. The pilot therefore completed
all 16 registered seeds and did not touch tuning `[16,32)` or held-out `[32,64)`.

**Track R gate: PASS.** Every fixed cliff has material completion mass at 500M;
the either-corner standalone rates range from 10/16 to 15/16, and every fixture
also has positive same-seed union gain over its better corner. This admits later
policy tuning but is not itself a rescue-policy or live-controller result. No
restart policy was implemented or run.

## F. Authority, schema, and equivalence evidence

Before the pilot:

- all 12 cliff fixtures regenerated byte-identically at exact requested counts;
- all 64 frozen seeds matched independent Python and product SplitMix64;
- the schema suite passed, including configured cutoff vs termination and
  mixed-rig rejection;
- authority matched direct product at tiny and 500M controls;
- completing standalone/replay controls matched outcomes, rewards, layouts, and
  decisions; and
- focused cleanup tests passed normal, failure, exception, and timeout exits.

The pilot contains authority rows only. Authority owns outcome, completion
inferences, censoring, reward, and layout. No instrumented row was collected or
pooled because authority outcomes already differentiated; therefore no progress
counter is called cutoff-equivalent.

Whole-file validation reports `schema: OK (131 rows; mixed-rig guard passed)`.
That invocation found and corrected a benchmark-only `argparse` typo
(`nargs="-"` to `nargs="*"`); a CLI regression test now proves file validation
starts and succeeds. The row contract itself was unchanged.

## G. Files, status, and risks

- `benchmarks/probe_arrange/run_pr0.py`: base-locked, resumable serial authority
  launcher with operation-ID deduplication, per-row heartbeats, and CPU cap.
- `benchmarks/probe_arrange/analyze_pr0.py`: exact-matrix validator and summary.
- `benchmarks/results/2026-07-17-p-r0-pilot.jsonl`: 131 authority rows.
- This report and the probe README document commands and interpretation.
- `schema.py` and `test_schema.py` contain only the standalone CLI fix/test.

Residual risks: 16 frozen streams are a finite benchmark, standalone corners do
not reproduce operation-wide memo/PRNG accounting, and success curves below
500M are derived from observed authority completion counts rather than separately
rerunning each cutoff. Wall and CPU are host-specific accounting, not search
adjudicators. No product-equivalence claim extends beyond the preflight controls.

## H. Draft ledger entry

### P-R0 - fixed-instance seeded runtime distributions - PROBE, GATE TRACK R IN

- Ran unchanged full-perturbation standalone authority at both non-transpose
  corners for four fixed cliff instances over frozen search seed indices
  `[0,16)`, plus one completing control per grid size. Fixture seed and search
  seed remained separate; 131/131 rows passed the corrected JSONL schema.
- Cliff authority produced 60/128 placements and 68/128 right-censored
  `not_proven` outcomes, with no infeasible or interrupted rows. All controls
  completed. The early stop did not apply, so all 16 pilot seeds ran.
- Every fixture showed material completion mass (either-corner 10/16 to 15/16)
  and positive paired corner union gain (+1 to +4 seeds over the better corner).
  Gate Track R in for later same-seed policy tuning. Do not claim
  perturbation-arm complementarity or a live operation-wide controller result.
- Consumed 1,887.630 child CPU-seconds (0.524 hours), below the 1.5-hour pilot
  cap. No tuning or held-out seed was used; no restart policy, product code,
  baseline, history, identity, workload, or golden changed.
