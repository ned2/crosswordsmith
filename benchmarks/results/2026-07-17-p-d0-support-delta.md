# P-D0 support delta and proof multiplicity

Date: 2026-07-17
Base: `1bccf47917fea074a404bf467e51865ce988b8b8`
Branch: `probe/a-d0`
Verdict: **NOMINATE proof-preserving A-D2; reject geometry-only; A-D1 warranted**

## Rig and method

`benchmarks/probe_arrange/d0_support.pl` is a benchmark-only twin of the strict
`mrv_inc` path. It calls unchanged product `inc_counts/8`, ordering,
`find_intersecting_word/6`, `assign_word/9`, and mutable letter/boundary grids.
The shadow observer runs only after the product count map exists and never under
an inference limit. Its residue map is an ordinary threaded term restored by
backtracking; mutable holders contain counters and decision observations only.
No shadow value enters MRV ordering, branch choice, grid state, or a cache.

The proof scanner independently reproduces product crossing order and
multiplicity but deliberately does not call or populate the product
`pair_crossings/3` table. Every observed exact bucket is checked against the
unchanged product count map. Every claimed delta classification is checked
against that full exact capped recount.

The sound classifier tested here is intentionally narrow:

- Previous visible bucket 0: intervening non-letter-sharing placements cannot
  add a crossing source, so all additions are proofs through the newest word.
- Previous visible bucket 1: its one old proof is retained as a residue.
  Intervening non-sharing placements can only remove it; the newest word is the
  only new source. Newest-source proofs plus residue survival therefore give the
  exact capped bucket.
- Previous visible bucket 2: always use the current full recount. If two
  residues fail, unobserved older proofs may still survive, so a residue-only
  lower bucket is not sound.

For bucket 0/1, hypothetical candidate checks are the newest-source scan plus
one direct residue check when a bucket-1 scan returns fewer than two proofs.
Bucket 2 incurs exactly the current full recount. This is semantic work
accounting, not a product-inference or wall claim.

Controls are deterministic completing strict fixtures. Both current
non-transpose representatives were measured for two light and two dense
fixtures. Each row ran a warmup, unchanged authority, unobserved twin, and two
observed passes. The observed passes required exact result, reward, layout,
decision sequence, counters, and measured inference count; wall and CPU were
excluded. Full-operation controls separately retained one reset and one seed
stream across both representatives.

## Transitions, classification, and multiplicity

`Transitions` is a 3x3 matrix written as
`previous-0[current 0/1/2]; previous-1[...]; previous-2[...]`. `Proofs/geom`
counts capped legal proof solutions versus unique `(Start,Dir)` geometries.
`Residue survive` counts surviving old proof residues over all old residues.

| control | corner | transitions | classified | proofs/geom | divergence events | residue survive | exact -> proposed checks |
|---|---|---|---:|---:|---:|---:|---:|
| light 09x09/08w | TLA | `7/5/2; 6/2/4; 6/3/2` | 26/37 (70.27%) | 26/25 | 1 | 8/34 | 251 -> 144 |
| light 09x09/08w | TR | `2/3/2; 1/2/3; 2/0/3` | 13/18 (72.22%) | 21/21 | 0 | 6/16 | 104 -> 63 |
| light 15x15/12w | TLA | `10/5/7; 4/5/2; 6/2/6` | 33/47 (70.21%) | 42/42 | 0 | 11/39 | 248 -> 164 |
| light 15x15/12w | TR | `46/11/1; 41/11/4; 8/12/10` | 114/144 (79.17%) | 64/63 | 1 | 30/116 | 1,588 -> 583 |
| dense 15x15/32w | TLA | `16/38/19; 23/42/46; 18/46/162` | 184/410 (44.88%) | 580/580 | 0 | 332/563 | 5,318 -> 3,335 |
| dense 15x15/32w | TR | `147/73/28; 93/58/45; 29/48/111` | 444/632 (70.25%) | 547/542 | 5 | 234/572 | 17,870 -> 4,122 |
| dense 21x21/80w | TLA | `251/233/69; 199/268/180; 63/172/1160` | 1,200/2,595 (46.24%) | 3,491/3,470 | 21 | 2,124/3,437 | 70,562 -> 33,325 |
| dense 21x21/80w | TR | `262/242/82; 197/211/163; 91/133/1214` | 1,157/2,595 (44.59%) | 3,504/3,501 | 3 | 2,183/3,447 | 50,897 -> 26,790 |

Class aggregates:

| class | classified refreshes | rate | exact -> proposed checks | check delta | divergences | proofs/geom |
|---|---:|---:|---:|---:|---:|---:|
| light | 186/246 | 75.61% | 2,191 -> 954 | -56.46% | 2 | 153/151 |
| dense | 2,985/6,232 | 47.90% | 144,647 -> 67,572 | -53.28% | 29 | 8,122/8,093 |

All 3,171 classifications matched the full exact recount. Every light sentinel
reduced candidate checks individually: -42.63%, -39.42%, -33.87%, and -63.29%.
The 31 proof/geometry divergences reject any geometry-deduplicated design.
Repeated proofs through distinct placed-word sources are part of the visible
Cap=2 bucket and must remain distinct.

## Dirty watches and saturation

For each old proof residue, `letter` watches are its run and non-edge endpoint
letter reads, `boundary` watches are its run cells in the boundary grid, and
`adjacency` watches are perpendicular letter cells read for then-empty run
cells. Dirtiness is intersection with cells newly bound by only the newest
placement: new letter bindings for letter/adjacency, and new trailed boundary
marks for boundary. The observer binds none of these cells.

| control | corner | letter dirty/watch | boundary dirty/watch | adjacency dirty/watch | bucket 0 events/checks | bucket 1 events/checks | saturated events/checks |
|---|---|---:|---:|---:|---:|---:|---:|
| light 09x09/08w | TLA | 24/168 (14.29%) | 10/122 (8.20%) | 26/166 (15.66%) | 19/162 | 10/64 | 8/25 |
| light 09x09/08w | TR | 10/83 (12.05%) | 4/55 (7.27%) | 11/72 (15.28%) | 5/26 | 5/49 | 8/29 |
| light 15x15/12w | TLA | 36/210 (17.14%) | 10/146 (6.85%) | 21/208 (10.10%) | 20/123 | 12/79 | 15/46 |
| light 15x15/12w | TR | 163/659 (24.73%) | 31/435 (7.13%) | 50/629 (7.95%) | 95/1,233 | 34/290 | 15/65 |
| dense 15x15/32w | TLA | 334/3,053 (10.94%) | 90/2,010 (4.48%) | 153/2,743 (5.58%) | 57/1,130 | 126/2,633 | 227/1,555 |
| dense 15x15/32w | TR | 468/3,152 (14.85%) | 128/2,102 (6.09%) | 227/2,923 (7.77%) | 269/11,288 | 179/5,360 | 184/1,222 |
| dense 21x21/80w | TLA | 1,956/18,248 (10.72%) | 436/11,739 (3.71%) | 711/16,160 (4.40%) | 513/24,155 | 673/25,608 | 1,409/20,799 |
| dense 21x21/80w | TR | 1,939/18,624 (10.41%) | 508/12,039 (4.22%) | 563/16,640 (3.38%) | 550/17,729 | 586/17,331 | 1,459/15,837 |

Class-level dirty-cell rates were light `20.80% / 7.26% / 10.05%` and dense
`10.90% / 4.17% / 4.30%` for letter/boundary/adjacency watches. The low rates
support selective event work conceptually, but they do not authorize pruning:
the exact source/residue classifier above is the sound seam.

## A-D1 direct-state work

These are operations performed by the unchanged count-map/order path, counted
at every interior node. Gets split into ordering-map reads and carried-count
reads; refresh proof work is separate.

| class | interior nodes/rebuilds | assoc puts | ordering gets | carry gets | sort calls/items | initial full counts | refreshes |
|---|---:|---:|---:|---:|---:|---:|---:|
| light | 81 | 356 | 356 | 63 | 81/150 | 47 | 246 |
| dense | 322 | 7,544 | 7,544 | 1,092 | 322/5,798 | 220 | 6,232 |

A-D1 remains warranted as an isolated experiment. Dense controls rebuild and
write 7,544 assoc entries, read 8,636 assoc entries, and sort 5,798 items over
only 322 interior nodes. Direct stable buckets may remove meaningful state work
without changing recount semantics. This is a nomination, not an expected win:
the strict 0.5% light-rung gate remains binding, and direct storage must preserve
stale non-sharing overcounts and input-order ties before A-D2 is layered on it.

## Equivalence and validation

- Full-operation light control: unchanged authority, unobserved twin, and
  observed twin all placed reward 79 with the same layout; the exact combined
  two-corner decision sequence had 22 entries.
- Full-operation dense control: all three placed reward 880 with the same
  layout; the exact combined two-corner decision sequence had 160 entries.
- All eight corner rows matched unchanged authority outcome, reward, and layout
  and matched the unobserved twin decision-for-decision.
- Duplicate observed passes matched every non-wall field and measured inference
  count exactly on every row.
- `make probe-arrange-check`: fixture, seed, schema, and four focused P-D0 tests
  pass.
- `make test`: 418 assertions pass; all goldens and CLI contracts pass.
- Strict core+heavy identity: all 15 selected rows exact, including latency
  digest `90289af7...6c81b`.
- Strict core+heavy ratchet: all 14 gated counts exact at `+0.00%`; latency-only
  row remains informational.
- Greedy identity: exact `ffc32b71...f1c29`.
- Greedy core+heavy ratchet: every construction/sweep/postprocess metric exact
  at `+0.00%`.
- No product, baseline, history, identity manifest, or golden file changed or
  was recorded.

## Gate verdicts and risks

**A-D2: NOMINATE, proof-preserving design only.** The pre-registered 25% gate
passes on both classes (75.61% light, 47.90% dense), all classifications are
exact, and candidate checks do not increase on any light sentinel. A production
experiment must retain proof source/multiplicity, classify only previous bucket
0/1, and use the unchanged full recount for bucket 2. It remains conditional on
A-D1 clearing its own ratchet.

**Geometry-only delta: KILLED.** Thirty-one recounts had more legal proof
solutions than unique `(Start,Dir)` geometries. Deduplication would change the
visible bucket and potentially MRV order.

**A-D1: WARRANTED.** Direct stable buckets have enough measured assoc/sort work
to justify one isolated implementation experiment. The probe does not establish
an inference win; representation overhead and backtrack restoration remain the
risks.

Residual risks:

- Candidate-check savings are semantic counts. Product inference and wall wins
  must be measured on the actual A-D2 implementation after A-D1.
- Direct residue lookup/storage can recreate E-H10b's small-rung admission tax.
  A-D1 must supply the previous bucket and residue without a new assoc probe in
  the letter-sharing loop.
- A-D2 must preserve newest-source proof order and Cap=2 early saturation. It
  cannot replace visible stale buckets with exact non-sharing counts.
- Watch dirtiness is measured against only the newest placement, not a license
  for a global event cache or a boundary sentinel in the letter grid.

## Draft ledger entry

### P-D0 - support delta and proof multiplicity - MEASUREMENT (A-D2 nominated)

- **Rig/soundness:** added a benchmark-only exact-replay observer after
  unchanged product `inc_counts/8`. Pure threaded residues and independent
  proof scans cannot influence MRV, grids, tables, or budget semantics. Previous
  buckets 0/1 are classified from newest-source proofs plus the sole old residue;
  bucket 2 always falls back. Every claim was checked against a full exact
  Cap=2 recount.
- **Result:** classified `186/246` light refreshes (75.61%) and `2,985/6,232`
  dense refreshes (47.90%). Hypothetical candidate checks changed
  `2,191 -> 954` light and `144,647 -> 67,572` dense, with reductions on all
  four light corner controls. All 3,171 classifications were exact.
- **Multiplicity:** 31 recounts had proof count greater than unique geometry
  count (`153/151` light, `8,122/8,093` dense). Geometry dedup is therefore
  rejected; distinct source proofs remain load-bearing.
- **A-D1:** measured 7,544 assoc writes, 8,636 assoc reads, and 5,798 sorted
  items over 322 dense interior nodes. A-D1 remains warranted as a standalone
  direct-state experiment, with zero light regressions mandatory.
- **Equivalence/verification:** eight two-corner-class rows and light/dense
  full-operation controls matched unchanged outcome, layout, reward, and exact
  decision order; duplicate counters were exact. Focused probe checks, full
  native tests, strict/greedy identities, and both core+heavy ratchets pass.
  No product or baseline artifact changed.
- **Verdict:** nominate proof-preserving A-D2 only after A-D1 passes. Kill the
  geometry-only design. Revisit bucket-2 deltas only with a separately proven
  complete proof representation; two watched residues are insufficient.
