# A-D1 stable IDs and direct trailed buckets

Date: 2026-07-17
Base: `362b3998f84610719b17421cfa1f18a309b15f9a`
Branch: `experiment/a-d1`
Verdict: **KEEP candidate**

## A. Scope and design

The strict two-representative search now assigns each input entry its positional
integer ID and wraps, but never copies, the original entry as `wid(ID, Entry)`.
The count state is one `buckets/N` term whose argument `ID` is the visible
0/1/2 bucket for that word. Its final argument is an unused, unshared variable,
following the local SWI-Prolog manual's `setarg/3` guidance against accidental
sharing. Every count mutation uses ordinary trail-restored `setarg/3`; there is
no `nb_setarg`, global, dynamic fact, or non-backtrackable product state.

At each interior node the remaining tagged words are scanned in input-relative
order. Bucket 0 is omitted from candidates; bucket-1 and bucket-2 difference
lists are concatenated. This is exactly the stable order produced by the old
`Count-Entry` pairs, positive filter, and stable `keysort/2`. Seeded search still
uses the same permutation for the seed and one permutation for each non-empty
equal-count bucket, in the same order.

Only `arrange.pl`'s strict `corner_search/4` calls the new driver. The existing
assoc-backed `assign_words_inc/9` remains the fragment, enumerate, strategy,
and old-reference path. Greedy/best-effort/candidates are unchanged. The one
top-level memo reset, two non-transpose representatives, shared operation
budget, proof generator, letter grid, and boundary grid are unchanged. No
A-D2 residue/newest-source logic, geometry deduplication, cache, or exact
non-sharing recount was added.

## B. Old/new correspondence

Let input position `i` identify entry `E_i`, whose unique answer is `A_i`.
For every remaining word at every selection node, the representation invariant
is:

```text
arg(i, Buckets, C)  iff  get_assoc(A_i, OldCountAssoc, C)
```

The invariant holds as follows:

1. At the first interior node after a seed, every remaining bucket argument is
   fresh and receives the same capped `mrv_count(2, ...)` as the old full assoc
   build.
2. After choosing word `W`, each remaining entry sharing a letter with `W` is
   recounted by the unchanged proof generator and its argument is overwritten
   with that exact 0/1/2 result.
3. A non-sharing entry receives no write. Its argument therefore retains the
   prior visible value, exactly matching the old assoc carry. This deliberately
   preserves stale overestimates rather than exposing the current exact count.
4. Removed IDs are irrelevant to the invariant. Remaining tagged entries retain
   input-relative order, so partitioning yields bucket-1 ties then bucket-2 ties
   in exactly old stable-sort order.
5. Child updates occur after the parent's branch choicepoint. SWI-Prolog's
   ordinary `setarg/3` trails an overwrite and restores it when backtracking
   returns before the call, so sibling branches see the parent's buckets.

IDs are derived from position before search, not by `==` lookup. A focused test
copies entries containing non-ground metadata and proves both copies receive IDs
1,2 while their metadata variables and compound entries remain non-identical.

## C. Differential and equivalence evidence

`benchmarks/probe_arrange/ad1_buckets.pl` is a benchmark-only old/new replay.
The assoc side uses unchanged `inc_counts/8`, assoc lookup, positive filtering,
`keysort/2`, and `order_candidates/3`. The direct side uses the candidate
representation. Neither side has an inference limit. Observational holders use
non-backtrackable mutation only to retain traces; no observed value enters
search, grid state, count state, or branch selection.

Every trace event compares exactly, not as a set: remaining `ID-Count` values,
ordered candidate IDs, selected IDs, legal placement decisions including failed
attempts, final placement, reward, and layout signature.

| control | corner | count nodes | selections | decisions | reward | result |
|---|---|---:|---:|---:|---:|---|
| 09x09/08w | TLA | 13 | 15 | 14 | 79 | exact |
| 09x09/08w | TR | 7 | 8 | 8 | 79 | exact |
| 15x15/12w | TLA | 14 | 15 | 15 | 122 | exact |
| 15x15/12w | TR | 47 | 40 | 48 | 122 | exact |
| 15x15/32w | TLA | 34 | 35 | 35 | 347 | exact |
| 15x15/32w | TR | 130 | 104 | 131 | 349 | exact |
| 21x21/80w | TLA | 79 | 82 | 80 | 878 | exact |
| 21x21/80w | TR | 79 | 80 | 80 | 880 | exact |

The existing P-D0 old-reference observer also passes all four focused tests
against the new product authority, including exact two-corner decision order and
proof multiplicity. A small exhaustive lock counts all solutions from the old
assoc and new direct drivers as `3 == 3`. The direct export's first solution was
probed with `deterministic/1` and correctly leaves a choicepoint, matching its
`nondet` PlDoc.

All 15 strict identity rows match. The budget-saturating latency identity remains
exactly:

```text
90289af7db529b0132bae8bb910a18e90daa6f7b1c19de7cb97508883a56c81b
```

## D. Strict inference ratchet

Two complete serialized core+heavy runs produced the same candidate count on
every row. Old values are the protected baseline; neither run recorded them.
Rows are listed in the campaign request's order.

| strict rung | old | new run 1 | new run 2 | delta |
|---|---:|---:|---:|---:|
| 09x09/08w | 25,006 | 23,598 | 23,598 | -5.63% |
| 09x09/16w | 662,138 | 643,211 | 643,211 | -2.86% |
| 09x09/17w | 38,503,164 | 37,787,386 | 37,787,386 | -1.86% |
| 15x15/12w | 95,381 | 90,239 | 90,239 | -5.39% |
| 15x15/28w | 341,666 | 322,772 | 322,772 | -5.53% |
| 15x15/32w | 929,466 | 900,159 | 900,159 | -3.15% |
| 15x15/34w | 13,601,499 | 13,473,896 | 13,473,896 | -0.94% |
| 15x15/36w | 38,215,934 | 37,666,380 | 37,666,380 | -1.44% |
| 15x15/40w | 10,338,363 | 10,210,252 | 10,210,252 | -1.24% |
| 21x21/25w | 280,528 | 265,129 | 265,129 | -5.49% |
| 21x21/80w | 4,246,246 | 4,049,012 | 4,049,012 | -4.64% |
| 21x21/82w | 6,555,246 | 6,339,469 | 6,339,469 | -3.29% |
| real 13x13/12w | 3,545,338 | 3,469,476 | 3,469,476 | -2.14% |
| real 15x15/18w | 248,068 | 239,765 | 239,765 | -3.35% |

There are 14 wins and zero regressions. The lightest win is -0.94%, beyond the
strict 0.5% threshold without an absolute-floor exception. The informational
latency row measured `500,008,136` in both runs against reference `500,000,191`;
its output identity is unchanged and its budget-pinned inference count does not
gate.

## E. Greedy neutrality

Greedy identity remains exact:

```text
ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29
```

All construction/sweep/postprocess references reproduced at +0.00%:

| greedy rung | construction | sweep | postprocess |
|---|---:|---:|---:|
| bundled17 | 12,811 | 127,545 | 6,910 |
| bundled11 | 2,098 | 10,037 | 29 |
| benchmark08 | 43,676 | 164,166 | 5,974 |
| real13 | 36,087 | 358,826 | 90 |
| real15 | 98,467 | 627,096 | 130 |
| ladder32 | 368,867 | 1,241,744 | 233 |
| ladder80 | 239,344 | 4,910,962 | 593 |

## F. Secondary wall and RSS evidence

After deterministic gates passed, `run_ad1_wall.pl` ran one warm pair and 21
serialized alternating assoc/direct pairs per dense control. Each pair required
exact old/new operation result before it entered the sample. Wall is secondary:

| control | assoc median | direct median | paired median ratio | change |
|---|---:|---:|---:|---:|
| 15x15/32w | 0.047183 s | 0.043862 s | 0.933452 | -6.65% |
| 21x21/80w | 0.233434 s | 0.213035 s | 0.919810 | -8.02% |

The normal ratchet's candidate process RSS medians were generally below the
committed values, but the committed baseline is from a different host
(`lana` versus `bunk`), so no RSS delta is adjudicated. No profile claim was
needed after deterministic counts and paired wall agreed.

## G. Verification, files, and risks

Verification completed:

- `make test`: 435 assertions, all goldens, and all CLI/stderr contracts pass.
- A-D1 differential: eight corner rows exact for counts/order/decisions/results.
- P-D0 focused observer: 4/4 pass against new product authority.
- strict identity `--heavy`: all 15 rows exact.
- strict ratchet `--heavy`: 14 wins, zero regressions, reproduced twice.
- greedy identity and core+heavy ratchet: exact.
- `git diff --check`: pass.
- No baseline, history, identity manifest, golden, or authoritative ledger was
  changed or regenerated.

Candidate files:

- `prolog/crosswordsmith/core.pl`
- `prolog/crosswordsmith/arrange.pl`
- `tests/core.plt`
- `tests/probe_arrange.plt`
- `benchmarks/probe_arrange/ad1_buckets.pl`
- `benchmarks/probe_arrange/run_ad1_wall.pl`
- this report

Residual risks:

- The strict path now has a second MRV driver beside the retained assoc driver;
  future selection changes must update both and keep the differential seam.
- `setarg/3` is extra-logical and SWI-specific. The term is deliberately local,
  guarded against accidental sharing, and covered by sibling restoration.
- The WASM battery was not part of this experiment brief and was not run; it is
  a promotion-time structural-core gate.

Variants deliberately not tried: answer-assoc ID lookup, copied-entry identity,
persistent mutable bucket lists, exact non-sharing recount, A-D2 residues,
newest-source deltas, geometry deduplication, and caches. They would confound the
registered storage-only mechanism or violate known semantics.

## H. Draft ledger entry

### A-D1 - stable IDs and direct trailed buckets - KEEP

- **Change/soundness:** strict two-representative MRV now wraps original entries
  with stable positional IDs and stores visible capped counts in one guarded
  fixed-arity term. Letter-sharing IDs receive the unchanged full recount via
  ordinary trail-restored `setarg/3`; non-sharing IDs receive no write and keep
  the old visible stale overestimate. A stable one-pass 1/2 partition exactly
  replaces positive filtering plus stable `keysort/2`. Fragment, enumerate,
  greedy, grids, proof order, memo lifecycle, and budget policy are unchanged.
- **Equivalence:** a benchmark-only assoc/direct replay matched every visible
  ID-count snapshot, candidate order, selected ID, legal placement decision,
  result, reward, and layout on eight light/dense corner controls. P-D0's old
  observer remains exact, a small full tree counts `3 == 3`, all 15 strict
  identities match, and greedy identity remains `ffc32b71...f1c29`.
- **Result:** two full strict core+heavy runs reproduced 14 wins and zero
  regressions. Improvements range from -0.94% to -5.63%; representative dense
  rows are 15x15/32w `929,466 -> 900,159` (-3.15%) and 21x21/80w
  `4,246,246 -> 4,049,012` (-4.64%). Serialized paired wall ratios are 0.9335
  and 0.9198. Every greedy construction/sweep/postprocess cell is +0.00%.
- **Verification/verdict:** full native tests/goldens/CLI contracts, P-D0,
  strict and greedy identities, duplicate strict ratchets, greedy ratchet, and
  diff checks pass. No protected benchmark or output artifact changed.
  **KEEP.** A-D2 may now be evaluated separately from this representation base;
  do not fold residues or newest-source logic into A-D1's attribution.
