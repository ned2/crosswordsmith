# Spec: cryptic-style layout quality and optimisation

Status: **proposed** — design for review; not implemented. This is the third
iteration; §13 records what changed and why (two adversarial reviews). The
engine is **greedy density-construction + local improvement** (research doc
Approach B); branch-and-bound is demoted to an optional, later refinement.

Companion: [`docs/crossword-setting-and-generation-QA.md`](./crossword-setting-and-generation-QA.md) (research) and
`docs/experiments.md` (solver strategies; I2 value-ordering result; I5 maximality
fix; start-position bias).

## 1. Problem statement & use case

crosswordsmith arranges a **closed set** of a user's words (e.g. a site's links,
as a "fancy table of contents") into one interlocking grid. **No filler
dictionary** — every cell must be a letter of a user word, because every word is
a link. Desired aesthetic: the **cryptic / British lattice** look (about half the
cells checked, unchecked cells allowed, balanced/compact, no ugly adjacencies),
not the American 100%-checked block grid.

The crux (research doc Q3): an arbitrary closed set often **cannot** reach the
cryptic checking ratio (Scrabble effect / CSP phase transition). So the tool must
deliver the **best layout it can** and **report what it compromised**, while
letting a user impose hard requirements when they want them.

The current solver is a closed-set, no-filler, connected packer with the right
hard legality rules (crossing agreement, no flush adjacency, no word-inside-word
[I5], connectivity). It returns the **first** connected legal layout with **no**
notion of quality, and that first layout is structurally **tree-like** (each word
connected by a single crossing ⇒ ~2 checked cells/word regardless of length). So
the gap to "cryptic" is something to **construct toward**, not a filter to apply.

## 2. Goals

- **G1.** Construct layouts toward the cryptic look via an explicit **quality
  objective Q** (checking, compactness, unch quality), not first-solution.
- **G2.** Let a user impose any subset of cryptic rules as **hard floors**, via
  flag **groundedness** (bound ⇒ enforced; unbound ⇒ Q decides). One mechanism
  gives strict / manual / partial / auto.
- **G3.** Always return *something* useful on real (often sparse) inputs — the
  best connected layout achievable (dropping unplaceable words), with a
  **compromise report**.
- **G4.** Keep the existing fast solver + strategies intact; add quality as a
  **new driver**, not a rewrite.

## 3. Non-goals

- Filler/dictionary words (closed set; optional connectors are a later
  possibility — but see §12, the research doc ranks filler as the *first*
  relaxation, so this may need revisiting for very sparse sets).
- 180° symmetry / black-square grid construction (free-form).
- Clue generation (layout only).
- Provably-optimal Q (greedy + hill-climb is heuristic; the search is **anytime**
  under a **node/inference** budget, not wall-clock — for reproducibility).

## 4. Design decisions

- **D1. No filler.** Only user words (revisit for pathological sparsity, §12).
- **D2. Default = auto** (no hard floors; maximise Q; report). Strict/manual are
  opt-in via floors.
- **D3. Words may be dropped** when unplaceable, but **placing a word dominates
  Q** (the objective never prefers dropping a link — §5c); "place all" is also an
  optional hard floor (`all_words`).
- **D4. No symmetry.**
- **D5. Word length is fixed input.** "minLen" is **not** a hard floor (it would
  silently drop links); it is a soft Q penalty on a glut of very short words.

## 5. The model: hard invariants + optional floors + objective

Terms: a cell is **checked** iff occupied by both an across and a down word;
**unch** (for word w) iff in w only. `checkedCount(w)` = checked cells of w;
`L = len(w)`.

### 5a. Hard invariants (always; already in the solver — §10)

Crossing letters agree; unique answers; no flush end-to-end; no
parallel-adjacency; no word-inside-word (maximality, I5); connectivity (each
placed word crosses an already-placed word). Non-negotiable.

### 5b. Optional hard floors (the flags) — groundedness selects the mode

A `Floors` term (dict) carries one entry per relaxable rule, each **unbound**
(not imposed; Q handles it softly) or **bound to a value** (a hard constraint the
constructor must satisfy; placements that would make it permanently unsatisfiable
are skipped):

| floor | bound meaning | unbound |
|---|---|---|
| `min_half` | ∀ placed w: `checkedCount(w) ≥ ⌈L/2⌉` | Q rewards checking |
| `max_unch_run:k` | ∀ placed w: no run of `> k` consecutive unch cells | Q penalises long unch runs |
| `all_words` | every input word placed (no drops) | Q penalises drops (and drops are dominated anyway, §5c) |

Groundedness ⇒ mode, **no relaxation ladder/odometer** (the objective makes the
trade-offs):

- **Strict** = all floors bound. **Manual/partial** = bind what you require.
- **Auto (default)** = none bound; maximise Q; report.

Infeasibility honesty: when a *bound* floor can't be met, the tool reports it is
infeasible and gives the **cause it can cheaply establish** — specifically the
static per-word `min_half` supply check (§6) — and otherwise reports "no layout
met floor X" without over-claiming a per-constraint diagnosis the engine can't
give. (The research doc's "GAC names the failing light" needs per-light
propagators we don't have.)

Notes vs the source cryptic rules: `ends_checked` is **not** a separate rule (for
odd L it is implied by `min_half`; for even L the doc doesn't mandate it) — it is
folded into Q as a soft `uncheckedEndCount` penalty (§5c). `min_len` is a soft Q
term, not a floor (D5).

### 5c. The objective Q (always maximised; node-cost must stay cheap)

```
Q =  w_check · checkedFraction         # checked cells / filled cells -> the cryptic core
   - w_ends  · uncheckedEndCount       # words whose first/last cell is unch -> the lattice "ends-checked" look
   - w_unch  · unchRunPenalty          # Σ (runLen − 1)^2 over unch runs -> punishes 3+-in-a-row hard
   - w_bbox  · boundingBoxCost         # area + aspect penalty -> renderable on a web page
   - w_drop  · droppedWordCount        # DOMINANT: any all-words layout beats any layout that drops a link
   - w_short · shortGlut               # optional: glut of very short words
```

Three fixes baked in (from review):
- **Drops are dominant** (`w_drop` ≫ the others, effectively near-lexicographic):
  the engine prefers *any* all-words layout to *any* dropping layout; checking /
  compactness only break ties among all-words layouts, and only rank partials
  once all-words is genuinely infeasible. This kills the perverse "drop the hard
  links to raise `checkedFraction`/shrink the bbox" incentive.
- **`crossingsPerWord` removed** — redundant with `checkedFraction` (every
  crossing is a checked cell).
- **Terms normalised** so Q is comparable across layouts with different word
  counts and bounding boxes (e.g. `checkedFraction` is a ratio; `boundingBoxCost`
  is normalised area/aspect; `uncheckedEndCount`/drops are per-word counts scaled
  to be commensurable).

Weights are the tool's editorial taste; v1 ships defaults set from the v1a data
(§14). The per-node cost of evaluating Q-deltas must stay cheap (I2's lesson:
per-node cost is what kills these heuristics) — the greedy key (§6) is computed
incrementally per candidate placement, not by re-scoring whole layouts.

## 6. Engine: greedy density-construction + local improvement (anytime)

Primary engine (research doc Approach B). Chosen over branch-and-bound because it
returns a good connected partial in ~one pass, degrades gracefully *by
construction*, needs **no admissible Q upper-bound** (so the prior design's
correctness landmine disappears), and is **deterministic by node count**.

1. **Construct.** Seed the first word (ranging over `start_locs`, §8). Then
   repeatedly: among words not yet placed that **can connect now** (a legal
   crossing exists — reuse `find_intersecting_word/6`'s connectivity, so the
   connectivity rule is *respected by construction*, not fought), choose the
   placement that maximises the greedy key **(new crossings desc, bbox-growth
   asc)**, and place it. The greedy key is the density-seeking value heuristic;
   it is goal-aligned and — crucially — there is **no backtracking**, so I2's
   ~20%-per-node enumeration cost is a one-pass cost, not multiplied over a tree.
   A word with no legal placement *this round* is retried after others land; a
   word that can never connect is **dropped** (§7).
2. **Best-so-far is always defined (fixes G3).** Every connected legal state from
   the seed onward is a valid answer (remaining words implicitly dropped). The
   engine tracks the best-Q state seen, so it always has something to return —
   even if construction stalls early or the budget fires.
3. **Local improvement (hill-climb).** From the constructed layout, apply cheap
   Q-raising moves: re-place a word at a higher-Q crossing, swap which word
   occupies a contested cell, or **restart** with a different start/greedy
   tie-break and keep the best. Bounded by a node/inference budget.
4. **Optional B&B refinement (deferred; v2, only if v1a shows greedy leaves Q on
   the table).** A budgeted keep-best search over the high-value frontier;
   requires a *proven-admissible* Q upper-bound or runs without the prune. Not in
   v1b.

## 7. Dropping words (graceful degradation)

Dropping is intrinsic to greedy construction (stop placing when stuck), so G3 is
satisfied without a separate subset search:

- **Static pre-drop** (cheap): words with 0 shared letters (can never connect),
  or — under a bound `min_half` — words whose shared-letter supply `< ⌈L/2⌉`
  (statically `min_half`-infeasible), are dropped up front and listed in the
  report.
- **Constructive drop** (free): any remaining word the constructor can't place
  (it shares letters but doesn't fit the geometry that emerged — including the
  realistic "every link shares one common vowel, consumed once at a crossing"
  case) is simply left unplaced and reported.

So a word that "shares letters but can't co-exist" is handled — it just doesn't
get placed, and best-so-far (§6.2) still returns the rest.

## 8. Start position is a quality variable

Feasibility and layout shape are start-dependent (`experiments.md`: ~15× swing on
the one real-word fixture). Greedy construction is cheap, so the engine **runs
all four `start_locs/1`** (and could treat the seed as a search variable) and
keeps the best-Q result. The report records which start won. The node budget is
split across the four constructions (each is near-linear, so ×4 is affordable —
unlike the prior B&B design where ×4 split an exponential budget).

## 9. Compromise report

For the returned layout: achieved metrics vs cryptic ideals — checked fraction
(vs ~0.5), per-word `checkedCount/L`, longest unch run, bounding box / aspect,
which start won, which floors were met, and **dropped words with reasons**. Read
off the chosen layout (an output, not an input).

## 10. Integration

- New **driver** (e.g. strategy `quality`), parallel to `assign_words`/
  `assign_words_inc`; existing `baseline`/`mrv`/`mrv_capped`/`mrv_inc` and the
  golden output untouched (G4).
- Does **not** reuse the `mrv_inc` count cache (its invariant is about
  first-solution variable ordering, not checking deficits); the quality engine
  keeps its own per-word `checkedCount` state.
- Reuses legality (`find_intersecting_word`, `assign_word` incl. I5
  `no_word_merge`). The density-seeking key pushes *toward* the I5 hazard (tight
  collinear packing), so `no_word_merge` stays load-bearing; note (per I5/dense_16)
  the repeated propose-and-reject of substring placements is a real perf cost on
  dense sets, not negligible — keep the greedy key from proposing them where
  cheaply detectable.
- `assign_clue_numbers` is unaffected by dropping; protected by I5.

## 11. Edge cases (explicit)

- **1 word:** placed alone; `min_half` infeasible ⇒ if bound, report infeasible;
  else return it (low Q).
- **2 words, no shared letter:** place 1 (seed), drop the other; report.
- **All words share one common letter only** (e.g. all contain one `E`): static
  pre-drop passes them (each shares a letter), but the constructor places only as
  many as the single shared letter allows to interlock and **drops the rest**
  (§7 constructive drop) — best-so-far returns the connected core. Handled.
- **All same length:** I5-safe; floors may still conflict (`min_half` + small
  `max_unch_run`).
- **Substring word:** `no_word_merge` (I5) forbids the collinear placement; the
  greedy key may propose-and-reject it (perf cost on dense sets).
- **Large sets:** rely on the §6 node budget; greedy is near-linear per pass so
  this is far more benign than the prior B&B.

## 12. Open questions / risks

- **How close does one-pass greedy get to the cryptic ideal?** The core
  empirical question — answered by v1a measuring greedy-construction Q vs the
  current tree-like first solutions vs the ~0.5 checked target.
- **Greedy local optima / painting-into-a-corner:** ANSWERED. Multi-restart
  (grid x start x seed) suffices: the result is robustly 1- and 2-move locally
  optimal (§14 v1b), so a dense early core is not a practical problem and deeper
  search does not help.
- **Greedy-key weighting** (crossings vs bbox-growth) is a taste call; needs v1a
  data.
- **Very sparse sets** may force many drops; the research doc ranks "allow
  filler" as the least-damaging relaxation, which D1 forbids — revisit if v1a
  shows real link-sets routinely drop a lot.
- **Q-weight defaults / normalisation constants** need the v1a data.

## 13. Iteration history (what changed and why)

- **v1 (rejected):** "find a complete legal layout, accept-or-reject cryptic
  checks, relax flags via a backtracking ladder." Review: **enumerate-and-reject**,
  intractable on sparse sets; the `decide/1` ladder was lexicographic, not
  fewest-compromises; `ends_checked` mistranslated; `min_len` polarity wrong;
  compactness missing; flags conflated input and report.
- **v2 (rejected):** objective-driven **branch-and-bound**. Review: G3 not
  actually guaranteed (best-so-far only at complete layouts); the optimistic
  `Q⁺` prune-bound was **undefined and likely not cheaply admissible** (a
  silent-unsoundness landmine); the deficit heuristic was self-contradictory
  (deficit-targeted ∧ cheap ∧ connectivity-respecting can't all hold); Q had a
  **perverse drop incentive**; a time budget breaks the project's
  inference-count benchmarking discipline.
- **v3 (this revision):** **greedy density-construction + hill-climb** as the
  primary engine. This *dissolves* the v2 landmines — no admissible bound (no
  B&B in v1b), no exponential place-or-drop tree, graceful degradation by
  construction (G3 via best-so-far at every connected state), node budget for
  determinism. Surviving fixes applied: drops dominant in Q (kills the perverse
  incentive), `crossingsPerWord` removed, `uncheckedEndCount` term added,
  honest infeasibility reporting, start as a quality variable. B&B becomes an
  optional v2 refinement gated on v1a evidence.

## 14. Phased plan

### v1a — Layout-quality analyzer (measure first; cheap; no solver change)
A standalone analyzer over the emitted JSON computing the §5c/§9 metrics per
layout: checked fraction, per-word `checkedCount/L`, unch-run histogram,
crossings/word, bounding-box area + fill density + aspect, word-length stats,
connectable/dropped info. Run across the current strategies × four starts ×
fixtures **and realistic link-sets** (construct a couple). *Deliverables:*
(1) how far current first-solution layouts sit from the ~0.5 cryptic ideal — the
feasibility reconnaissance that sets Q's default weights and decides whether the
heavier engine is worth it; (2) the metric code v1b's Q reuses. *Why first:*
needed regardless, commits to nothing.

### v1b — Greedy quality engine (IMPLEMENTED — `quality.pl`, CLI `--quality`)
A `quality` engine: greedy density-construction (§6.1) with best-Q selection
(§6.2) across candidate **grid sizes** and starts and **multi-seed restarts**
(§6.3 hill-climb form), optional hard **floors** via groundedness (§5b:
`--min-half`, `--max-unch K`, `--all-words`; drop-to-satisfy + all-words reject),
**constructive drop** (§7), and the compromise report (§9, stderr). Reuses
crossword.pl legality/emit; existing strategies and golden untouched (66 tests
pass). Measured: reaches/beats the witness ceilings on the dense quality
fixtures and places all of `toc_demo`'s real words; floors give a cryptic-valid
core or honest infeasibility.

Polish done: a **compactness-weighted Q** (`quality_weights/2`: checked vs
bbox-area vs elongation) renders a near-square TOC (aspect 2.12->1.09 on
quality_22 for ~0.02 less checked); a **bounded grid sweep** (three candidate
sizes + the per-set seed cap) keeps it fast (quality_61 ~21 s -> ~5.6 s).

Polish NOT adopted: **local-move search** was tested and is a **no-op**. First a
single-word reseat (re-place a word at its now-best crossing): the multi-restart
search (grid x start x seed) already reaches *single-reseat-optimal* layouts.
Then a **2-word eject-and-reinsert** (the natural escalation): on quality_22, 448
candidate two-word moves were generated and the **best one ties the current Q**
(no improvement); identical on toc_demo. So the greedy+restart layout is robustly
**1- and 2-move locally optimal**, and under the compactness-weighted Q it is at
or above the witness ceilings. Both reverted (cf. I2).

**Conclusion on quality-guided *search*:** it is **not worth building**. The
greedy construction + multi-restart already reaches the local optimum, and 1-/2-
move search can't improve it, so a heavier optimiser (B&B; spec v2 below) has no
measured headroom on the quality fixtures. (The witness ceilings are pure-density
and ignore compactness; greedy's slightly-lower checked fraction comes *with*
better compactness, which our Q rewards - so greedy is at our Q-optimum, not
short of it.)

### v1b.1 — The construction is cut-free (the cut-based version was removed)
The greedy construction is implemented **cut-free** ("Path A" — a declarative
rewrite, *not* CLP(FD)). The original had two control constructs; both are gone:
the first-solution `!` in `word_best_placement` is replaced by "fold `assign_word`
legality INTO the `findall` generator, then take the head of the `@>=`-sorted
list"; the greedy `( Best -> place ; stop )` if-then-else is replaced by "reify
the best move (or `none`) as a term and dispatch on mutually exclusive clause
heads (`[]` vs `[_|_]`)". No `!`, `->` or `\+` remain in the construction;
`findall`/`sort`/`is` stay. `quality.pl`: `greedy_loop`, `next_move`/`best_move`/
`apply_move`, `word_best_placement`.

This was first landed as an opt-in `--quality-engine pure` variant and A-B'd
against the cut-based original. It produced **byte-identical** output on every
quality fixture (JSON + stderr report, floored and unfloored; verified by diff)
**and** ran faster:

| fixture          | cut-based            | cut-free             |
|------------------|----------------------|----------------------|
| toc_demo (16w)   | 2014 ms / 26.7M inf  | 1484 ms / 15.1M inf  |
| quality_22 (22w) | 1335 ms / 18.5M inf  |  859 ms / 10.1M inf  |
| quality_61 (61w) | 5504 ms / 69.8M inf  | 3924 ms / 37.1M inf  |

(SWI 10.0.2; CPU + inferences per `quality_layout` call, N=30/30/10.) ~26–36 %
less CPU, ~43–47 % fewer inferences.

**Why faster (the prediction was wrong):** removing the green cut was expected to
cost *more* — run `assign_word` for every candidate, not just up to the first
legal one. It costs *less* because folding `assign_word` into the generator moves
the legality filter **before** scoring: the cut-based version computed the
(heavier) `placement_key` density score for *every* `find_intersecting_word`
candidate, including the many illegal ones, then assigned only the first legal;
the cut-free version runs a cheap early-failing `assign_word` on each and scores
only the survivors. On dense meshes most crossing candidates are illegal, so
pruning scoring early wins. Since it is identical-output, faster **and** more
declarative, the cut-based construction and the `--quality-engine` flag were
**removed** — the cut-free engine is now the only one.

### v2 — Refinements
The **quality-guided B&B search** is **closed** (the local-move experiment above
shows no headroom over greedy+restart). Remaining optional refinements: richer Q
(length variety, themed/priority positioning); exposed Q-weights; optional
non-link connectors for pathologically sparse sets (§12).

## 15. v1a results (done)

Analyzer: `benchmarks/analyze_layout.py` (reads emitted JSON → metrics). Run
across current strategies × starts × fixtures plus two real-word sets
(`bundled_17_clues` 6 words; `toc_demo` 16 words). Cryptic ideal:
`checked_fraction ≈ 0.5`, `frac_words_meeting_min_half ≈ 1.0`, no 3+ unch runs.

| input (current first-solution layout) | checked_frac | meet_½ | longest_unch | words w/ 3+run | aspect |
| --- | --- | --- | --- | --- | --- |
| bundled_17 (6 real words) | **0.089** | **0.0** | 8–12 | 4–5 of 6 | 1.06 |
| toc_demo (16 real words) | **0.135** | **0.0** | 6–8 | 13–15 of 16 | 1.05–1.92 |
| benchmark_08 (8, synthetic) | 0.182 | 0.0 | 6 | all | 1.0 |
| benchmark_14 (14) | 0.259 | 0.0 | 4 | all | 1.0 |
| benchmark_16_dense (16, rigged dense) | 0.308 | 0.0 | 2 | 0 | 1.0 |
| benchmark_20 / _26 (combs) | 0.07 / 0.06 | ~0.04 | 12 / 14 | nearly all | 2.85 / 3.27 |
| benchmark_70_mesh (70, rigged dense) | 0.353 | **0.657** | 4 | 10 of 70 | 1.0 |

Findings:
- **Current layouts are nowhere near cryptic on real inputs** — ~0.09–0.14
  checked vs the ~0.5 ideal, **0% of words ≥half-checked**, long bare runs (8–12)
  and a 3+-unch run in nearly every word. The solver stops at the first
  connected tree-like layout and leaves enormous checking on the table. **The gap
  is large ⇒ v1b (greedy density engine) is warranted.**
- **Density correlates directly with checking** (rigged-dense fixtures reach
  0.31–0.35; combs ~0.06). This validates greedy *density*-construction: pushing
  toward density is pushing toward checking — and there is clear headroom over
  the current 0.09–0.14.
- **But ~0.5 is likely unreachable for arbitrary closed sets with no filler**:
  even maximally dense rigged fixtures top out ~0.31–0.35. Realistic target is
  "maximise achievable checking + compactness and **report vs the ideal**", not
  "hit 0.5". If real link-sets need higher checking, the research doc's
  least-damaging relaxation — optional non-link **connectors** (currently D1
  no-filler) — becomes the key v2 lever.
- **Compactness is an independent quality dimension** (aspect 1.0 → 3.3 across
  inputs/starts), confirming Q's `boundingBoxCost` term is doing real work.
- Start position barely moves `checked_fraction` but can swing aspect (1.05 →
  1.92 on `toc_demo`) — another reason §8 ranges over starts and keeps best Q.

## 16. Engine cost vs the solver (apples-to-apples)

The quality engine and the solver (best variant `mrv_inc`) do **different jobs**,
so there is no single "Nx slower" number: the solver finds *one* valid layout on
a *given* grid (first-solution search, can backtrack); the quality engine
*constructs* the best dense layout across a *multi-restart sweep* and *picks its
own grid*. Measured on the same word sets — the solver fed the grid the quality
engine chose, the quality cost split into a single greedy construction vs the
full sweep (SWI 10.0.2; inferences are the portable metric; every row placed all
words with no drops):

| fixture (words) | solver `mrv_inc`, 1 search | quality, 1 construction | quality, full sweep (the deliverable) |
| --- | --- | --- | --- |
| toc_demo (16) | 690k inf / 38 ms | 183k inf / 32 ms | 15.0M inf / 1734 ms — 60 builds |
| quality_22 (22) | 227k inf / 13 ms | 280k inf / 19 ms | 10.1M inf / 814 ms — 36 builds |
| quality_61 (61) | **45.9M inf / 2462 ms** | **405k inf / 42 ms** | 37.1M inf / 3693 ms — 12 builds |

Findings:
- **Like-for-like (one construction vs one search), the quality core is
  competitive-to-better and far more robust.** A single construction is the same
  order as a single solver search on the small fixtures, and ~**110x cheaper** on
  the dense one (quality_61: 405k vs 45.9M inf). Greedy construction never
  backtracks, so it scales gracefully exactly where the solver's first-solution
  search explodes (the same blow-up as `dense_16`, here on a tight 61-word grid).
- **The quality engine's real cost is the sweep, and the multiplier is the build
  count** (`grids × 4 starts × seeds` = 12–60 here). That spend is not waste — it
  buys what the solver does not: density/checking (~0.45 vs the solver's ~0.1, §15),
  automatic grid selection, and drop tolerance.
- **No single ratio — it inverts with input hardness.** On solver-*easy* inputs
  (quality_22, a trivial 227k-inf search) the full sweep is ~45x more expensive;
  on solver-*hard* inputs (quality_61) the full sweep (37.1M) is actually
  *cheaper* than one solver search, and a single construction is ~110x cheaper.
- **Takeaway:** for "a valid layout on a grid I pick", one search ≈ one
  construction and the construction will not blow up; for "the best dense
  closed-set layout, engine picks the grid", the sweep is the price of quality —
  ~1–45x a single solver search depending on whether the input is solver-easy or
  solver-hard, and the *cheaper* option precisely on the dense inputs where the
  solver struggles.
