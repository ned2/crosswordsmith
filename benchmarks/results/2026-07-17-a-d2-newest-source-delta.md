# A-D2 proof-preserving newest-source delta

Date: 2026-07-17
Base: `83b7d89315e5b8107a3fdd22df2154c7de1a4501`
Branch: `experiment/a-d2`
Verdict: **KEEP candidate**

## A. Base, scope, and branch

The worktree started clean at the required A-D1 ratchet commit. A-D1's stable
IDs and direct trailed buckets were active only through strict
`corner_search/4`; fragment, enumerate, strategy, assoc, greedy, best-effort,
and candidates paths remained on their existing drivers. The baseline contained
all 14 A-D1 strict counts, and the integrated P-D0 report, rig, and tests were
present.

A-D2 changes only the strict direct driver. It does not modify protected
baselines, history, identity manifests, goldens, or the authoritative ledger.

## B. Exact classifier and residue design

The direct state now owns parallel fixed-arity `buckets/N` and `residues/N`
terms. Argument `ID` in the residue term is `none` unless the visible bucket is
1, in which case it is the sole ground proof:

```text
proof(SourceAnswer, SourceStart, SourceDir, CandidateStart, CandidateDir)
```

Both terms are local to one strict search and carry unused unshared guard
arguments. Product updates use ordinary trail-restored `setarg/3`; there are no
assoc probes, globals, dynamic facts, or non-backtrackable persistent values.
The fresh solution counter uses the existing local `nb_setarg/3` counting
pattern only while collecting up to two solutions; the ID-indexed count and
residue state is trailed.

The classifier is:

- Previous 0: scan proofs sourced through only the newest placed word. Since
  every intervening non-sharing word cannot be a crossing source, this is the
  complete exact bucket and, when 1, the exact residue.
- Previous 1: scan the newest source first. Two proofs saturate immediately. At
  0/1 newest proofs, revalidate the sole old proof against its exact source and
  current geometry; combine at Cap=2. A resulting bucket 1 retains whichever
  proof is first in the unchanged full order.
- Previous 2: always run the unchanged full recount. No residue lowering is
  attempted because unobserved older proofs may survive.
- Non-sharing words: perform no write, retaining both stale visible bucket and
  residue exactly.

The source-retaining generator duplicates the unchanged
`find_intersecting_word/6` conjunction verbatim: placed-source order, tabled
crossing order, and proof multiplicity are preserved. It never deduplicates
`(Start,Dir)`. Residue validation re-establishes source answer/start/orientation,
matching crossing positions and letter, finite-canvas fit, boundary legality,
cell letters, and adjacency before calling a proof surviving.

One attempted implementation detail was dropped: making
`find_intersecting_word/6` a thin wrapper around the source-retaining generator
added 0.06%-0.28% to untouched greedy inference metrics. Duplicating the
conjunction only for A-D2 restored every greedy metric exactly.

## C. P-D0 controls and equivalence

`benchmarks/probe_arrange/ad2_delta.pl` drives the product delta refresh and,
before any result can affect ordering, checks it against unchanged
`mrv_count/8` and a full proof recount. Its exact traces contain previous,
full, and delta buckets/residues plus ordered count snapshots, selected IDs,
and legal placement decisions. The A-D1 full driver and A-D2 delta driver match
count-for-count and decision-for-decision on every row.

| control | corner | classified | fallback | exact -> delta checks | equivalence |
|---|---|---:|---:|---:|---|
| light 09x09/08w | TLA | 26/37 (70.27%) | 11 | 251 -> 144 | exact |
| light 09x09/08w | TR | 13/18 (72.22%) | 5 | 104 -> 63 | exact |
| light 15x15/12w | TLA | 33/47 (70.21%) | 14 | 248 -> 164 | exact |
| light 15x15/12w | TR | 114/144 (79.17%) | 30 | 1,588 -> 583 | exact |
| dense 15x15/32w | TLA | 184/410 (44.88%) | 226 | 5,318 -> 3,335 | exact |
| dense 15x15/32w | TR | 444/632 (70.25%) | 188 | 17,870 -> 4,122 | exact |
| dense 21x21/80w | TLA | 1,200/2,595 (46.24%) | 1,395 | 70,562 -> 33,325 | exact |
| dense 21x21/80w | TR | 1,157/2,595 (44.59%) | 1,438 | 50,897 -> 26,790 | exact |

Aggregates reproduce P-D0 after A-D1 exactly: light `186/246` classified
(75.61%), 60 fallback, and `2,191 -> 954` candidate checks (-56.46%); dense
`2,985/6,232` classified (47.90%), 3,247 fallback, and
`144,647 -> 67,572` checks (-53.28%). All 6,478 product refreshes and all 3,171
classifications matched the unchanged full recount. The known 31
proof/geometry divergences remain proofs, not deduplicated geometries.

Focused validation also covers previous-0 transitions, previous-1 residue
survival and invalidation, source/orientation validation, proof multiplicity,
bucket-2 full fallback, stale non-sharing retention, paired bucket/residue trail
restoration, non-ground stable IDs, and full-tree counts at both strict corners.
Seed 42 full-vs-delta traces match at both corners of 15x15/12w.

## D. Strict ratchet, reproduced twice

The protected A-D1 counts are the `old` column. Both complete serialized runs
returned the same new count on every row.

| strict rung | old | new run 1 | new run 2 | absolute delta | relative delta |
|---|---:|---:|---:|---:|---:|
| 09x09/08w | 23,598 | 22,297 | 22,297 | -1,301 | -5.51% |
| 09x09/16w | 643,211 | 431,007 | 431,007 | -212,204 | -32.99% |
| 09x09/17w | 37,787,386 | 19,700,703 | 19,700,703 | -18,086,683 | -47.86% |
| 15x15/12w | 90,239 | 74,355 | 74,355 | -15,884 | -17.60% |
| 15x15/28w | 322,772 | 240,483 | 240,483 | -82,289 | -25.49% |
| 15x15/32w | 900,159 | 627,744 | 627,744 | -272,415 | -30.26% |
| 15x15/34w | 13,473,896 | 7,928,254 | 7,928,254 | -5,545,642 | -41.16% |
| 15x15/36w | 37,666,380 | 14,074,891 | 14,074,891 | -23,591,489 | -62.63% |
| 15x15/40w | 10,210,252 | 4,480,185 | 4,480,185 | -5,730,067 | -56.12% |
| 21x21/25w | 265,129 | 216,653 | 216,653 | -48,476 | -18.28% |
| 21x21/80w | 4,049,012 | 3,107,968 | 3,107,968 | -941,044 | -23.24% |
| 21x21/82w | 6,339,469 | 3,770,401 | 3,770,401 | -2,569,068 | -40.52% |
| real 13x13/12w | 3,469,476 | 2,091,663 | 2,091,663 | -1,377,813 | -39.71% |
| real 15x15/18w | 239,765 | 172,537 | 172,537 | -67,228 | -28.04% |

There are 14 wins and zero regressions. The latency-only row measured
`500,008,116` in both runs versus informational reference `500,008,136`.
All 15 strict identity digests match; the latency digest remains exactly:

```text
90289af7db529b0132bae8bb910a18e90daa6f7b1c19de7cb97508883a56c81b
```

## E. Greedy neutrality

Greedy identity remains:

```text
ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29
```

Every construction/sweep/postprocess triplet is exactly +0.00%:

| greedy rung | construction | sweep | postprocess |
|---|---:|---:|---:|
| bundled17 | 12,811 | 127,545 | 6,910 |
| bundled11 | 2,098 | 10,037 | 29 |
| benchmark08 | 43,676 | 164,166 | 5,974 |
| real13 | 36,087 | 358,826 | 90 |
| real15 | 98,467 | 627,096 | 130 |
| ladder32 | 368,867 | 1,241,744 | 233 |
| ladder80 | 239,344 | 4,910,962 | 593 |

## F. Secondary wall evidence

After deterministic gates passed, one warm pair and 21 serialized alternating
A-D1-full/A-D2-delta pairs ran per dense control. Every pair required exact
operation result equality.

| control | A-D1 full median | A-D2 delta median | paired median ratio | verdict |
|---|---:|---:|---:|---|
| 15x15/32w | 0.035402 s | 0.035561 s | 1.001270 | null (+0.13%) |
| 21x21/80w | 0.173974 s | 0.176553 s | 1.002430 | null (+0.24%) |

No wall win is claimed. Process wall readings in the ordinary ratchet generally
improved on larger rungs but are too coarse and include startup; the paired
in-process result is the authority for this secondary evidence.

## G. Tests, WASM, files, and risks

Verification completed:

- `make test`: 453 assertions, all goldens, CLI exits, and stderr contracts pass.
- A-D2 differential: all eight P-D0 controls exact for full/delta buckets,
  residues, candidate order, selections, decisions, reward, and layout.
- Escalation: both-corner small full-tree counts and seeded full/delta traces pass.
- Strict identity `--heavy`: all 15 digests exact.
- Strict core+heavy ratchet: 14 wins, zero regressions, reproduced twice.
- Greedy identity and all seven exact metric triplets pass.
- `make test-wasm`: value parity, type lock, headless SDK, worker errors, and all
  spare-worker policies pass.
- `git diff --check`: pass.
- No baseline, history, identity manifest, golden, or authoritative ledger was
  changed or regenerated.

Candidate files:

- `prolog/crosswordsmith/core.pl`
- `tests/core.plt`
- `tests/probe_arrange.plt`
- `benchmarks/probe_arrange/ad2_delta.pl`
- `benchmarks/probe_arrange/measure_ad2.pl`
- `benchmarks/probe_arrange/run_ad2_wall.pl`
- this report

Residual risks:

- Paired wall is null despite the large inference reduction; proof descriptors
  and source validation trade Prolog calls for term work not reflected uniformly
  in wall time.
- The ordinary process benchmark's 21x21/80w peak RSS rose about 5 MiB (roughly
  22%-23%) on this host, still far below the campaign's 32 MiB additional-memory
  gate. The local state is bounded by input word count and request lifetime.
- The strict path now retains one proof only for bucket 1. Any future attempt to
  lower bucket 2 remains unsound without complete older-proof materialization.
- `setarg/3` remains SWI-specific; both local terms have unshared guards and are
  covered by sibling restoration and WASM tests.

## H. Draft ledger entry

### A-D2 - proof-preserving newest-source delta - KEEP

- **Change/soundness:** strict A-D1 direct state now carries a parallel local,
  ID-indexed trailed residue term. Previous bucket 0 scans only the newest
  source; previous bucket 1 combines newest-source proofs with one fully
  revalidated old source proof; previous bucket 2 always uses the unchanged
  full recount. Proof source/order/multiplicity, stale non-sharing buckets,
  stable ties, seed stream, memo lifecycle, corners, and shared budget remain
  exact. Geometry dedup and bucket-2 lowering were not attempted.
- **Differential:** all 6,478 P-D0 refreshes matched unchanged `mrv_count/8` and
  a full proof recount. Full-vs-delta count, selection, decision, result, reward,
  and layout traces match on eight controls plus both-corner seed-42 paths.
  Classification reproduces P-D0 exactly: 75.61% light and 47.90% dense, with
  candidate checks `2,191 -> 954` and `144,647 -> 67,572`.
- **Result:** duplicate full core+heavy runs produced 14 wins and zero
  regressions, from -5.51% to -62.63%. All 15 strict identities and the latency
  digest are exact. Paired dense wall ratios are null at 1.0013 and 1.0024; no
  wall claim is made. Greedy identity and every metric remain +0.00%.
- **Verification/verdict:** full native tests/goldens/CLI contracts, focused
  residue locks, small exhaustive counts, seeded traces, duplicate ratchets,
  strict/greedy identities, greedy metrics, WASM battery, and diff checks pass.
  No protected artifact changed. **KEEP.** The unexpectedly large deterministic
  win is supported by exact old/full/delta traces rather than inferred from
  geometry or output identity alone.
