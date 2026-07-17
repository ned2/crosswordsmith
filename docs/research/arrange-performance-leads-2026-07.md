# Research: next performance leads for `arrange` (2026-07)

Status: literature reconnaissance complete; Phase 0 and greedy A-G1/A-G2
accepted; Phase 2 premise probes adjudicated; A-D1 next.

## Execution update (2026-07-17)

A-G1 kept legality-first candidate filtering: dense sweep inferences fell from
`6,126,445 -> 2,335,326` (-61.88%) and `28,942,375 -> 9,362,822` (-67.65%).
A-G2 then kept exact transpose-partner synthesis: direct attempts halved and
those sweeps fell again to `1,241,744` (-46.83%) and `4,910,962` (-47.55%),
with raw-pool and output identity preserved. P-D0 then classified 75.61% of
light and 47.90% of dense recounts exactly, nominating A-D1 and proof-preserving
A-D2 while rejecting geometry dedup. P-C0 found 87-89% hard-state dead revisits
but zero parent-local captures, killing A-C1 and admitting a bounded global
cache candidate. P-R0 placed 60/128 fixed cliff rows and gates Track R; policy
tuning remains behind the explicit product checkpoint. A-D1 is the active next
serial experiment.

The authoritative record is `docs/experiments.md`, with detailed evidence in
`benchmarks/results/2026-07-17-a-g1-legality-before-score-premise.md`,
`benchmarks/results/2026-07-17-a-g2-transpose-premise.md`, and
`benchmarks/results/2026-07-17-a-g2-transpose-product.md`, plus the P-D0, P-C0,
and P-R0 reports under `benchmarks/results/2026-07-17-p-*.md`.

This note consolidates five research threads opened after the July arrange
performance campaign: incremental support maintenance, restart portfolios,
nogood/state caching, alternative formulations, and greedy/symmetry work. It
turns the literature into a ranked experiment queue while preserving the
campaign's measured negative results.

Companions:

- `arrange-perf-campaign-2026-07.md` records the shipped 71--84% inference
  reduction and the remaining policy trade around E-H10b.
- `arrange-search-algorithms.md` records the prior tree-search review.
- `../experiments.md` is the authoritative experiment ledger.
- `balafoutis-thesis-review-2026-07.md` explains why fill-style adaptive CSP
  machinery mostly does not transfer to the current arrange model.

## Load-bearing premises

Any new work must fit these measurements:

1. Ordinary strict-search cost is dominated by capped placement counting, not
   deep backtracking. P1 found 84--86% of retreats at distance 1--2 and no
   backtracking at all on the 21x21/80-word rung.
2. Near-cliff failure is a different problem. P2/P3 found that larger budgets
   and E-H9's large node-cost cut did not flip representative failures. Moving
   the envelope requires a different search trajectory.
3. The current visible count cache deliberately permits stale overcounts for
   non-letter-sharing words. They are safe because they never hide a placeable
   word, but making them exact can change golden-visible MRV ordering.
4. Candidate enumeration order is output-visible. A transparent optimization
   must preserve proof order, not merely the set of `(Start,Dir)` geometries.
5. The 0.5% relative inference ratchet remains strict. E-H10b is a measured
   large-grid win parked behind a product-evidence trigger, not a precedent for
   silently adding a small-rung tax.

## Ranked experiment queue

| Rank | ID | Lead | Scope | Evidence level |
|---:|---|---|---|---|
| 1 | A-G1 | Check greedy legality before scoring | best-effort/candidates | Strong code/profile evidence; transparent |
| 2 | A-G2 | Derive greedy transpose partners instead of searching them | best-effort/candidates | Read-only fixture probe passed; transparent |
| 3 | A-R0 | Measure fixed-instance seeded runtime distributions | strict near-cliff | Required premise probe |
| 4 | A-C0 | Measure duplicate child/state work | strict hard corners | Required premise probe |
| 5 | A-D0 | Measure newest-source deltas and placement-proof multiplicity | strict counting | Required premise probe |
| 6 | A-D1 | Stable word IDs plus direct trailed count buckets | strict counting | Literature-backed representation bet |
| 7 | A-R1 | Fixed-cutoff/Luby/geometric restart tournament | strict rescue path | Conditional on A-R0 |
| 8 | A-C1 | Parent-local failed-child deduplication | strict hard corners | Conditional on A-C0 |
| 9 | A-D2 | Residue plus newest-placed-word delta for bucket 0/1 | strict counting | Conditional on A-D0/A-D1 |
| 10 | A-T0 | Topology-first feasibility prototype | alternative strict model | Bounded research prototype |

A-G1 and A-G2 are the only current implementation candidates with enough
evidence to skip a separate instrumentation phase. Everything else begins as a
probe with an explicit stop condition.

## Thread 1: greedy construction

### What remains expensive

X3 already removed the dominant grid-copying error from the greedy path and
cut wall time by 17--18% on its heavy probes. The current residual pipeline is
still wasteful:

1. `find_intersecting_word/6` generates a crossing descriptor.
2. `placement_key/8` walks the word to count crossings, builds its cell run,
   and folds the run to compute bounding-box growth.
3. `check_word_fits/5` then rejects most candidates.

A reconnaissance profile of `ladder_21x21_80w` observed about 306,000 generated
crossing descriptors, 280,000 fully scored candidates, and only 7,074 legal
candidates. `placement_key/8` accounted for roughly 68% inclusive profile time;
only about 2.5% of fully scored candidates survived legality.

### A-G1: legality before score

Run `check_word_fits/5` before `placement_key/8`. Both are non-binding probes,
and `check_word_fits/5` already rejects negative starts, so the exception hazard
that motivated the old score-first order no longer requires paying the score
for rejected candidates.

Acceptance gate:

- Raw greedy-pool and CLI outputs remain byte-identical.
- Every inference-ratcheted greedy rung stays within 0.5%.
- Dense-rung score calls fall materially and wall time improves beyond noise.

### A-G2: derive transpose partners

The four greedy corners form the same two transpose pairs as strict search.
Best-effort/candidates must retain all four layouts because transposes are
intentionally distinct candidates, but they need not reconstruct both members.
Search one representative per pair, derive the transpose partner, and insert it
at the old position in the corner-major pool.

A read-only probe over bundled, real-word, and ladder fixtures found literal
ordered placement transposes for every seed tested, including identical dropped
word order. The production experiment must lock the complete raw pool, not only
the selected candidates.

Expected effect: halve greedy construction count while preserving candidate
semantics.

### Greedy benchmark prerequisite

The permanent arrange ratchet currently measures the strict MRV path. Add a
separate greedy ratchet with:

- Single-construction, full seed/corner sweep, pool-postprocessing, and CLI
  layers.
- Real-word anchors plus 15x15/32-word and 21x21/80-word heavy rungs.
- Inferences as the portable gate, semantic work counters for score/legality
  calls, and wall/RSS as host-specific evidence.
- SHA-256 locks for raw pool order, selected layouts, rewards, dropped order,
  candidate distances, and transpose pairs.

Direct CELF/lazy greedy is not sound here: placement scores can rise, new
placements become available, legality changes in both directions, and bbox
growth is not monotone. Later work may use admissible per-candidate bounds, but
must not treat stale scores as upper bounds.

## Thread 2: restarts and portfolios

The literature strongly supports restarts for heavy-tailed Las Vegas search:
Luby schedules are distribution-robust; geometric schedules avoid repeatedly
paying tiny cutoffs; and fixed cutoffs can win once the runtime distribution is
known. Portfolios help only when their arms are genuinely complementary.

The project has not yet established the key premise. P2 varied fixture-generator
seeds, producing different word sets. It did not hold one fixture fixed and vary
only `arrange --seed`. The density cliff is consistent with a heavy-tailed
search distribution, but does not prove one.

### A-R0: fixed-instance distribution probe

For each non-transpose corner, run at least 64 independently derived search
streams on fixed near-cliff fixtures. Record success inference, censoring,
maximum depth, unplacements, wipeouts, reward, and layout signature. Treat a
budget hit as right-censored data rather than a measured runtime.

Use easy/real controls, the envelope guards, and fixed instances just beyond
the 9x9, 15x15, and 21x21 cliffs. Keep fixture-generator seed and search seed in
separate result columns.

Proceed only if repeated search streams show useful short-run mass or material
corner/perturbation complementarity.

### A-R1: policy tournament

Under one unchanged operation-wide inference budget, compare:

- Empirically tuned fixed cutoffs.
- Paired Luby schedules over the two non-transpose corners.
- Paired geometric schedules.
- A static portfolio of root-only, equal-bucket-only, full, and unperturbed
  branch-order arms.

Derive every attempt stream independently from `(master seed, policy version,
corner, arm, attempt index)`. Do not let an earlier attempt's draw count remap
later attempts. Reset branch-local grid/count/PRNG state; retain only immutable
input memos, the incumbent, telemetry, and separately justified sound nogoods.

The acceptance objective is placement probability and restricted mean
inference-to-first-placement under the same total budget, not median node cost.
A cutoff is never evidence of infeasibility.

## Thread 3: dynamic placement support

The relevant literature is dynamic CSP and incremental view maintenance rather
than fixed-domain crossword fill. `pair_crossings/3` is a static potential
support relation; `check_word_fits/5` is its dynamic spatial filter;
`mrv_count/8` is a saturating aggregate over generated support proofs.

AC-4-style complete inverse support materialization is probably too expensive.
AC-6/AC-2001 residues, event-driven propagators, reversible sparse structures,
and delta maintenance are better conceptual matches.

### A-D0: shadow measurements

For every letter-sharing refresh, record:

- Previous visible bucket and exact recounted bucket.
- Old residue survival.
- Legal proofs generated against only the newest placed word.
- Unique `(Start,Dir)` geometries versus enumeration proofs.
- Bucket transition matrix and candidates examined before saturation.
- For event watching, which letter/boundary/adjacency cells each proof reads and
  what fraction of watches one placement dirties.

Proof multiplicity is load-bearing: the current counter may count one geometry
more than once when it is witnessed through multiple already-placed words.
Never deduplicate by geometry without measuring and intentionally changing that
contract.

### A-D1: direct trailed buckets

Assign stable integer word IDs and store visible 0/1/2 buckets in a fixed-arity
term updated with backtrackable `setarg/3`, following the successful E-H2/E-H9
pattern. Stable bucket partitions can preserve input order without rebuilding
an assoc and sorting three possible keys at every node.

Build and judge this storage change alone. It must be neutral on the light
rungs before support caching is layered onto it.

### A-D2: newest-source delta

For a previous bucket-1 word with a surviving residue, a legal proof against
the newly placed word certifies bucket 2. If no new proof exists, the bucket
remains 1. If the residue failed, fall back to full recount. Related logic may
cover bucket 0 after proof multiplicity is understood.

This differs from E-H10: it explicitly examines the only new crossing source
instead of incorrectly assuming the support set only shrinks. Stop if it cannot
classify a substantial fraction of recounts without increasing light-rung work.

## Thread 4: duplicate states and nogoods

Equivalent partial layouts are reachable because different legal placement
orders can converge and because one `(Word,Start,Dir)` child may be generated
through multiple crossed words. A small exhaustive reconnaissance case produced
64 derivations but only 22 canonical absolute layouts. That establishes the
phenomenon, not its value in first-solution production search.

### A-C0: duplicate-work probe

Thread an order-independent dual fingerprint through a shadow DFS and observe,
without pruning:

- Same-parent duplicate child decisions.
- Revisited canonical placement sets.
- Revisited states already proved dead.
- Nodes contained in the repeated failed subtrees.

Record a state as dead only after normal exhaustive failure, never after a
budget escape. On a potential production hit, verify an exact canonical key;
hash equality alone can never prune.

Kill global caching if repeated dead states are below 1% of entries or account
for less than 5% of hard-corner node work. Prefer parent-local dedup if it
captures at least 80% of the avoidable repetition.

### A-C1: parent-local failed-child dedup

Within one parent, remember a `(WordID,Start,Dir)` child only after its subtree
has exhausted. If the same child is generated through another crossing proof,
skip the repeated failed subtree. This avoids full state canonicalization and
is the best caching candidate if A-C0 confirms the expected duplication shape.

Global negative-state caching would key the sorted placed assignments plus the
remaining-word set and grid size. Translation normalization is unsound on a
finite canvas, and letter-grid-only keys omit boundary and orientation state.

## Thread 5: alternative models

Crozzle is the closest named problem family: a closed word set is geometrically
interlocked for score. Published crossword CSP work is mostly fixed-grid
dictionary fill and does not directly model this engine.

### A-T0: topology-first construction

Precompute every matching-letter event
`(WordI,PosI,WordJ,PosJ)`. Search over selected crossing events and component
merges before assigning an absolute canvas position. Maintain orientation
parity, relative coordinate potentials, component bounding boxes, and selected
crossing connectivity. Reject contradictory coordinate cycles and components
whose span exceeds the canvas.

The crucial distinction is that a branch may merge any two components. A
rooted one-word-at-a-time implementation merely recreates the current tree.
Once all words are connected, translate the relative component onto the canvas
and validate it through the existing legality core; incidental legal crossings
need not have been selected as topology edges.

Reconnaissance sizing for the 15x15/36-word rung found roughly 13,470 absolute
word placements but only 1,314 matching-letter events. This does not prove the
topology tree is smaller, but makes a bounded prototype credible.

Prototype gates:

- Match exhaustive legality on several 4--8 word cases.
- Demonstrate a merge where neither component contains the anchor.
- Reproduce known incumbent feasibility/reward before seeking improvements.
- Stop if geometric legality remains almost entirely leaf-checked.

### Iterated threshold feasibility

After topology feasibility works, ask whether a layout with reward `>= K`
exists. High `K` can make the otherwise loose placement problem tight enough
for crossing-cardinality bounds to propagate. The threshold must remove events
or force crossings before complete layouts; if it only rejects leaves, stop.

A static absolute-placement CSP and deterministic one-worker CP-SAT model are
useful comparison oracles. They are not current shipping candidates: pairwise
placement compatibility is extremely loose, support tables are large, and an
external solver conflicts with the pure SWI/WASM target.

## Closed or low-value directions

Do not reopen without evidence that changes their measured premises:

- Plain dom/wdeg, static degree, or value ordering.
- Conflict-directed backjumping.
- Letter-presence bitmask prefilters.
- The E-H10/E-H10b two-witness design as previously structured.
- A loose branch-and-bound objective over the current DFS.
- Direct CELF or stale-key priority queues for greedy selection.
- Full dihedral symmetry: only main-diagonal transpose preserves directed word
  spelling generically.
- Translation-normalized strict-state caching on a finite canvas.
- Whole-DFS tabling over live variable grids.
- Exact-cover/DLX or conflict-graph packing without a credible connectivity and
  score propagation mechanism.

## Literature starting set

- Mohr and Henderson, AC-4 (1986):
  https://doi.org/10.1016/0004-3702(86)90083-4
- Bessiere, AC-6 (1994):
  https://doi.org/10.1016/0004-3702(94)90041-8
- Bessiere et al., AC-2001/3.1 (2005):
  https://doi.org/10.1016/j.artint.2005.02.004
- Schulte and Stuckey, efficient constraint propagation engines:
  https://arxiv.org/abs/cs/0611009
- Luby, Sinclair, and Zuckerman, universal restart schedules (1993):
  https://doi.org/10.1016/0020-0190(93)90029-9
- Gomes, Selman, and Kautz, randomization and rapid restarts (1998):
  https://cdn.aaai.org/AAAI/1998/AAAI98-069.pdf
- Gomes et al., heavy-tailed SAT/CSP search (2000):
  https://doi.org/10.1023/A:1006314320276
- Lecoutre et al., nogood recording from restarts (2007):
  https://www.ijcai.org/Proceedings/07/Papers/019.pdf
- Dechter and Mateescu, AND/OR search spaces (2007):
  https://doi.org/10.1016/j.artint.2006.11.003
- McKay, isomorph-free exhaustive generation (1998):
  https://doi.org/10.1006/jagm.1997.0898
- Minoux, accelerated greedy submodular maximization (1978):
  https://doi.org/10.1007/BFb0006528
- Feo and Resende, GRASP (1995):
  https://doi.org/10.1007/BF01096763
- Forster, Harris, and Smith, the Crozzle automation problem (1992):
  https://doi.org/10.1145/143559.143598
- Gower and Wilkerson, R-by-C Crozzle NP-hardness (1996):
  https://doi.org/10.1145/331119.331148
- Wilson, crossword compilation using integer programming (1989):
  https://doi.org/10.1093/comjnl/32.3.273

## Recommended execution order

1. Build the permanent greedy benchmark and run A-G1, then A-G2.
2. Run A-R0, A-C0, and A-D0 as independent instrumentation probes.
3. Admit A-R1, A-C1, and A-D1/A-D2 only when their corresponding probe clears
   its stop condition.
4. Run A-T0 as an isolated bounded prototype, not as a product rewrite.
5. Ratchet every transparent win immediately; keep deliberate search-policy or
   output changes behind explicit product decisions.
