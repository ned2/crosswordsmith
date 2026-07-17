# Arrange performance campaign close-out

Date: 2026-07-17

## A. Base, branch, and scope

- Required measurement base: `5dcf9754382073ddec621e8ef2e89b0590c131c6`.
- Branch: `probe/arrange-closeout`.
- Measurement code commit: `c65df0a9fbe799391c90a1931e1faf7a8fc327d9`.
- SWI-Prolog: 10.1.10.
- Product state: A-G1, A-G2, A-D1, and A-D2 are present. The last product
  changes to `core.pl`/`arrange.pl` are A-D2 `9d7055a`, A-D1 `6bb50c8`, A-G2
  `7a14082`, and A-G1 `ef2e0b0`.
- Closure state: A-C2 and A-T0 reach this base only through ledger/plan/research
  documentation. `git diff 9d7055a..5dcf975 -- prolog/crosswordsmith/core.pl
  prolog/crosswordsmith/arrange.pl` is empty, so no rejected cache or topology
  candidate is in product.
- Scope: new benchmark-only current-direct attribution, a serial authority
  launcher, this report, and the durable research summary. No product file,
  baseline, history, manifest, golden, ledger, or campaign plan was changed or
  re-recorded.

The result has two deliberately separate evidence classes. Counter-free
authority alone determines operation outcome, reward, successful inference
count, and censoring under 500M. Current-direct attribution runs outside any
inference cap and supplies only semantic tree-shape counters.

## B. Current-direct replay proof

`benchmarks/probe_arrange/closeout_direct.pl` is a local twin of the current
A-D2 direct recursion. It calls product `refresh_delta_counts/9`,
`direct_partitions/4`, `direct_order/4`, `find_intersecting_word/6`, and
`assign_word/9` verbatim. Observer mutation records events after product choices
and cannot enter counts, ordering, grid state, or branch selection. Each corner
is a standalone mechanism request with one memo reset and no inference cap.

Counter definitions are:

- `nodes`: nonterminal recursive entries, including an entry that later wipes
  out.
- `legal decisions` and `places`: successful `assign_word/9` calls; these must
  be equal.
- `unplaces`: a legal child whose recursive continuation exhausts before the
  final first solution.
- `wipeouts`: an interior selection whose current positive bucket order is
  empty.
- `max depth`: maximum number of simultaneously placed words.
- `selections`: selected word IDs before crossing/legality; stale positive
  counts can make this exceed legal placements.

Proof was layered rather than inferred from output identity:

1. Focused easy, seed-42, and dense controls compare the unchanged current
   direct product result with the close-out twin twice.
2. They also compare the full selected-ID and legal-decision lists with the
   existing A-D2 differential, which validates every delta count/residue against
   the unchanged full recount before it can affect ordering.
3. Every hard row below repeats result, counters, and complete traces exactly.
4. Every hard row matches product outcome, reward, and absolute numbered layout
   exactly. The compact layout and trace SHA-1 values identify the full ground
   terms emitted by `run_closeout_direct.pl`.

Command:

```sh
swipl -q benchmarks/probe_arrange/test_closeout_direct.pl
swipl -q benchmarks/probe_arrange/run_closeout_direct.pl
```

Both passed. The focused verifier reported `3 replay controls passed`; all six
hard rows set `product_result_exact`, `ad2_trace_exact`, and
`duplicate_counters_exact` to `true`.

## C. Hard-corner attribution and cost model

| fixture | corner | nodes | selections | legal decisions / places | unplaces | wipeouts | max depth | reward | absolute layout SHA-1 | selected+decision trace SHA-1 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| 15x15/34w | `topleft_across` | 143 | 101 | 143 / 143 | 109 | 69 | 34 | 371 | `0ff5fcf65665e21ef1e3e52aa59aaf606100a714` | `d66e4f9b141c6ba0718cafb8cb66eed670a89d28` |
| 15x15/34w | `topright` | 3,107 | 2,270 | 3,107 / 3,107 | 3,073 | 1,757 | 34 | 366 | `0b0a023ef0e14a24cf3579d2b655748d0b31a510` | `fb4fe4bd0e86fb242c7359dc747eccbdcc370609` |
| 15x15/36w | `topleft_across` | 40 | 40 | 40 / 40 | 4 | 4 | 36 | 395 | `9a003286d796e846cc8c3a738baff3e46d7e31f5` | `27027021e6989d1f9510f2bf72d5c5a4de21169f` |
| 15x15/36w | `topright` | 3,454 | 3,418 | 3,454 / 3,454 | 3,418 | 2,117 | 36 | 390 | `16ca11c0c77627ff7851e32ee6feb9cbf91b1ced` | `0488b874f2373e05e25c169d18b55c4c7d0bbfdd` |
| 21x21/80w | `topleft_across` | 80 | 82 | 80 / 80 | 0 | 0 | 80 | 878 | `262b7070ad955c123adda6edca61724ff22f1df0` | `5cc4cbb3272108a600c8bd68795e2b9b5aa49187` |
| 21x21/80w | `topright` | 80 | 80 | 80 / 80 | 0 | 0 | 80 | 880 | `9ceaee2349f0628a8bcf0da4e4b59a9e46ffb676` | `8767268d952fc5e99fd598d290c0d30a231257f9` |

The current tree confirms the registered cost model:

- 15x15 hard cost is shallow but highly corner-asymmetric. The 34-word
  `topright` tree has 21.7 times the nodes of `topleft_across`; the 36-word
  ratio is 86.4 times. Both difficult right corners reach full word depth, but
  nearly every legal placement is later undone and 56%-61% of nodes are empty-
  bucket wipeouts.
- The easier representative is fixture-specific. 34w still performs 109
  unplacements, while 36w is almost a straight descent with four unplacements.
- 21x21/80w is not a backtracking workload at either representative: exactly 80
  nodes and placements, depth 80, zero unplacements, and zero wipeouts. The two
  extra TLA selections are stale-positive candidates rejected before a legal
  placement. Its product cost is therefore predominantly candidate counting and
  proof work, not tree churn.
- A-D2 changes the cost of maintaining those choices, not the choices. Exact
  selected/decision trace equality is consistent with its proof-preserving
  design and prevents an inference reduction from being misread as a smaller
  search tree.

## D. Robust 500M envelope authority

`run_closeout_authority.py` invokes the existing counter-free
`authority-operation` runner serially. Every row is the unchanged default
unseeded two-representative operation: one memo reset and one operation-wide
500,000,000-inference budget. Committed fixtures are loaded with an exact word-
count assertion.

| guard | asserted words | outcome | reward | authority inferences | budget used | inference headroom |
|---|---:|---|---:|---:|---:|---:|
| 9x9/17w seed 11 | 17 | placed | 182 | 19,699,727 | 3.94% | 480,300,273 |
| 15x15/40w seed 11 | 40 | placed | 443 | 4,478,919 | 0.90% | 495,521,081 |
| 21x21/82w seed 11 | 82 | placed | 916 | 3,763,253 | 0.75% | 496,236,747 |

All three shipped envelope guards remain robust. Their rewards reproduce the
historical product rewards, and each places far below the cap. The 9x9 guard is
the tightest but still leaves 96.06% of the operation budget unused. No marginal
near-500M placement is promoted by this result.

## E. Representative fixed cliffs

The unseeded operation authority used the P-R0 committed fixtures requested for
9x9 and 15x15. For 21x21 it used 88w seed 11 because it is the deeper P-R0
instance and P2 already established that this exact instance stays unproven at
2B. These are observations, not pass/fail gates.

| fixed cliff | asserted words | outcome | termination | measured wrapper inferences | reward | censoring |
|---|---:|---|---|---:|---|---|
| 9x9/18w seed 12 | 18 | `not_proven` | budget | 500,000,068 | unavailable | right-censored |
| 15x15/44w seed 11 | 44 | `not_proven` | budget | 500,000,087 | unavailable | right-censored |
| 21x21/88w seed 11 | 88 | `not_proven` | budget | 500,000,130 | unavailable | right-censored |

None is infeasible. The small over-cap wrapper counts are expected bookkeeping
around the product inference limiter, not extra search allowance. Together with
P-R0's 60/128 seeded standalone placements, these unseeded observations retain
the corrected premise: the fixed near-cliff regime is trajectory-sensitive and
search-bound, while ordinary robust placements benefit from cheaper counting.

## F. Cumulative deterministic impact and wall boundaries

### Comparison bases

The strict pre-campaign ratchet was `085a44c`; its 12 synthetic counts are
carried unchanged into Phase-0 commit `3ce0b8d`. Phase 0 changed benchmark
infrastructure only and added two real-word anchors plus the latency-only row.
For a complete 14-row comparison, the table uses `3ce0b8d`, the first full
comparable ratchet immediately before A-G1. It is not an arbitrary old row: all
12 inherited values are the actual prior ratchet, and its two added real rows
are the values protected before any accepted 2026-07b product experiment.

No greedy ratchet existed before this campaign. The correct non-hypothetical
base is therefore also `3ce0b8d`, after Phase 0 fixed mode-phase accounting and
before A-G1. Using the earlier provisional `aa4d3c8` history row would include
known benchmark accounting errors, not product improvement.

### Strict search inferences

These are full operation search inferences. Greedy changes are strict-neutral,
so the cumulative strict movement is A-D1 plus A-D2.

| strict rung | protected pre-candidate base | final accepted baseline | cumulative delta |
|---|---:|---:|---:|
| 09x09/08w | 25,006 | 22,297 | -10.83% |
| 09x09/16w | 662,138 | 431,007 | -34.91% |
| 09x09/17w | 38,503,164 | 19,700,703 | -48.83% |
| 15x15/12w | 95,381 | 74,355 | -22.04% |
| 15x15/28w | 341,666 | 240,483 | -29.61% |
| 15x15/32w | 929,466 | 627,744 | -32.46% |
| 15x15/34w | 13,601,499 | 7,928,254 | -41.71% |
| 15x15/36w | 38,215,934 | 14,074,891 | -63.17% |
| 15x15/40w | 10,338,363 | 4,480,185 | -56.66% |
| 21x21/25w | 280,528 | 216,653 | -22.77% |
| 21x21/80w | 4,246,246 | 3,107,968 | -26.81% |
| 21x21/82w | 6,555,246 | 3,770,401 | -42.48% |
| real 13x13/12w | 3,545,338 | 2,091,663 | -41.00% |
| real 15x15/18w | 248,068 | 172,537 | -30.45% |

All 14 rows improve, from 10.83% to 63.17% cumulatively.

### Greedy construction and sweep inferences

Postprocess inferences are unchanged on all seven rows. A-G1 reduces both the
pinned construction and full sweep by rejecting illegal candidates before
scoring; A-G2 then reduces only the full sweep by deriving transpose partners.

| greedy rung | construction base -> final | construction delta | sweep base -> final | sweep delta |
|---|---:|---:|---:|---:|
| bundled 17 candidates | 17,267 -> 12,811 | -25.81% | 296,711 -> 127,545 | -57.01% |
| bundled 11 best-effort | 2,120 -> 2,098 | -1.04% | 15,624 -> 10,037 | -35.76% |
| benchmark 08 candidates | 52,663 -> 43,676 | -17.06% | 356,853 -> 164,166 | -54.00% |
| real 13 best-effort | 68,276 -> 36,087 | -47.14% | 1,248,410 -> 358,826 | -71.26% |
| real 15 best-effort | 218,776 -> 98,467 | -54.99% | 2,789,556 -> 627,096 | -77.52% |
| 15x15/32w best-effort | 923,529 -> 368,867 | -60.06% | 6,126,445 -> 1,241,744 | -79.73% |
| 21x21/80w best-effort | 497,243 -> 239,344 | -51.87% | 28,942,375 -> 4,910,962 | -83.03% |

### Accepted wall evidence boundaries

Wall evidence is only claimed from serialized paired runs on the same host at
the time of each experiment:

| intervention | paired wall evidence | claim boundary |
|---|---|---|
| A-G1 legality before score | dense sweep median ratios 0.401 and 0.344 | positive wall evidence; 59.95% and 65.59% faster medians |
| A-G2 transpose synthesis | dense sweep ratios 0.5436 and 0.5500 | positive wall evidence; about 1.84x and 1.82x |
| A-D1 direct buckets | dense strict ratios 0.9335 and 0.9198 | positive secondary wall evidence; 6.65% and 8.02% |
| A-D2 newest-source delta | dense strict ratios 1.0013 and 1.0024 | null; no A-D2 wall speedup claimed |

The cumulative tables are deterministic inference evidence. They are not a
claim that paired wall ratios multiply, and ordinary process timing from ratchet
runs is not used to manufacture an A-D2 wall win. No wall comparison crosses
hosts.

## G. Final verification

| check | result |
|---|---|
| focused current-direct replay | PASS, 3/3 easy/seeded/dense controls |
| duplicate hard mechanism rows | PASS, 6/6 exact result/counters/traces |
| `make test` | PASS, 453 assertions; all goldens and CLI/stderr contracts |
| strict identity `--heavy` | PASS, all 15 rows; latency digest `90289af7...6c81b` |
| strict core+heavy ratchet, run 1 | PASS, all 14 gated rows and latency row exact `+0.00%` |
| strict core+heavy ratchet, run 2 | PASS, all 14 gated rows and latency row exact `+0.00%` |
| greedy identity | PASS, `ffc32b71cae9dfa55d532f000e639cf1a4c0feb272a1926c35d0dac7440f1c29` |
| greedy core+heavy ratchet | PASS, every construction/sweep/postprocess cell `+0.00%` |
| `make test-wasm` | PASS after installing the pinned lockfile dependencies; value parity, type lock, headless SDK, worker errors, and spare policies |
| protected artifacts | unchanged; no promote/record command run |
| `git diff --check` | PASS before the probe commit and again at final close-out |

The first WASM invocation stopped before tests because `wasm/test/node_modules`
was absent. `PLAYWRIGHT_SKIP_BROWSER_DOWNLOAD=1 npm ci` installed three pinned
packages from the committed lockfile; the rerun completed with `test-wasm OK`.

## H. Remaining policy and research items

- Track R remains a product-policy decision, not an unfinished transparent
  optimization. P-R0 gates a same-budget restart tournament in, but tuning and
  held-out seeds remain untouched. Revisit only after deciding whether an
  explicit rescue mode may return a different first layout/reward and accepting
  the registered seed/diagnostic contract.
- Never promote marginal near-500M placements. Revisit the robust envelope only
  after a shipped search-policy change, or when a fixed instance places with
  material headroom across the registered seed set under the same operation-
  wide budget.
- A-C1 remains closed at zero parent-local captures. Revisit A-C2 only with a
  measured exact lookup seam that removes the per-node admission tax on zero-
  dead controls; a size/depth switch, approximate key, or relaxed ratchet is not
  such evidence.
- A-T0 remains closed at gate 3. Revisit only with a sound nonlocal legality
  propagator that changes the measured 5,645 product-invalid complete-topology
  premise, not with event reordering or a leaf-only reward bound.
- Bucket-2 delta lowering remains out of scope. Revisit only with a complete
  proof representation that can establish surviving older proofs without
  changing proof multiplicity or candidate order.

## I. Draft final ledger entry

### Arrange 2026-07b campaign close-out - MEASURED, TRANSPARENT WORK COMPLETE

- **Scope/authority:** close-out at product base `5dcf975` on
  `probe/arrange-closeout`; counter-free current product authority alone supplies
  500M outcomes, rewards, inference counts, and censoring. A benchmark-only
  current A-D2 direct twin supplies tree counters outside any inference cap.
  Product, baselines, histories, identities, goldens, workloads, and plans are
  unchanged.
- **Tree shape:** exact current-direct replay shows 15x15 hard-corner asymmetry:
  34w TLA/TR `143/3,107` nodes and 36w `40/3,454`, with right-corner
  unplacements `3,073/3,418` and wipeouts `1,757/2,117`. Both 21x21/80w
  corners are straight 80-node descents with zero unplacement/wipeout. All six
  rows match product result/reward/absolute layout and exact A-D2 traces; duplicate
  counters are exact.
- **Envelope/cliffs:** robust unseeded operation guards place at
  `19,699,727/4,478,919/3,763,253` authority inferences with rewards
  `182/443/916`. Fixed unseeded cliffs 9x9/18w seed12, 15x15/44w seed11, and
  21x21/88w seed11 exhaust 500M and remain right-censored `not_proven`, never
  infeasible.
- **Cumulative impact:** from the protected pre-candidate ratchets, all 14 strict
  rows improve 10.83%-63.17%; greedy sweeps improve 35.76%-83.03% with
  postprocess unchanged. A-G1/A-G2 and A-D1 have paired wall support. A-D2's
  paired dense wall result remains null, so its accepted claim is inference-only.
- **Verification/verdict:** focused and duplicate replay, full native tests,
  strict identity, duplicate exact strict ratchets, greedy identity/ratchet,
  final WASM battery, and diff checks pass. A-G1/A-G2/A-D1/A-D2 are the final
  transparent product state; A-C1/A-C2/A-T0 remain closed. Track R remains parked
  at its explicit output-changing product-policy checkpoint.
