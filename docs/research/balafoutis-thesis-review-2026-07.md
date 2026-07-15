# Research: the Balafoutis thesis vs the crosswordsmith engines (2026-07)

Distilled from a full read (2026-07-16) of Thanasis Balafoutis, *Adaptive
Strategies for Solving Constraint Satisfaction Problems*, PhD thesis,
University of the Aegean, June 2011
(https://robotics.pme.duth.gr/balafoutis/wp-content/uploads/2022/11/Thesis.pdf).
A local copy lives (untracked) at `scratch/balafoutis_thesis.md` / `.pdf`. Why we read it: the ingrit_core
maintainer cites it as the inspiration for that tool's core algorithm, and our
own DP-8 fill core (§8.4c: MAC + dom/wdeg + geometric restarts, d-way
branching) independently landed on exactly the architecture the thesis
studies — so the thesis is, in effect, a menu of measured refinements around
the solver we already ship.

Scope of this note: (a) ranked candidates for the **fill** search core;
(b) why almost none of it transfers to **arrange**, now with empirical
backing from our own campaign probes; (c) the conceptual case for — and
against — re-modelling arrange as a CSP.

Companions: `docs/plans/fill-mac-dwd-implementation.md` (the DP-8/C2/C3 build
these candidates would amend), `docs/research/fill-perf-campaign-2026-07.md`,
`docs/research/arrange-search-algorithms.md` (prior arrange literature sweep;
its E-H6/P1/P2/P3 verdicts are load-bearing below), `docs/experiments.md`.

## The load-bearing framing result

**Every contribution in the thesis is learning-from-failure machinery over a
fixed binary constraint network.** Conflict weights attached to constraints,
dom/wdeg and its variants, weight-guided propagation ordering, branching
schemes over a variable's domain, restarts that carry learned weights — all
of it presupposes (i) variables, domains, and a constraint graph known before
search starts, and (ii) frequent, informative domain wipeouts (DWOs) for the
weights to learn from. Fill satisfies both premises; arrange satisfies
neither. That single observation organizes everything below.

The thesis's headline empirical fact for us (Ch 4): **plain dom/wdeg —
exactly what fill runs today — won zero instances outright across the
thesis's comparison tables** (§4.2, Table 4.5 discussion). The wins split
among the variants (aged dom/wdeg, H1/H2/H3, alldel, fully-assigned) with no
variant dominating and real per-class variance. Standard caveats: one
2010-era academic solver, binary extensional benchmarks (RLFAP, Golomb,
quasigroup, random), d-way + lexicographic values. Its experimental setup is
strikingly close to ours (MAC-3, variable-oriented FIFO propagation,
geometric restarts at growth factor 1.5), which raises transfer confidence,
but every payoff is still a bet to measure on the fill ladder, not a
guarantee.

## Chapter map (for future readers)

| Chapter | Content | Verdict for us |
|---|---|---|
| 3.4 / 4 | Conflict-driven VOH variants (H1/H2/H3, alldel, aging, fully-assigned, random probing) + broad evaluation | Candidates F2, F3, F4, F5 |
| 5 | Weight-guided revision ordering for AC-3 queues; dom/wdeg's sensitivity to revision order | Candidate F1 (best fit) |
| 6 | 2-way vs d-way branching; adaptive branching (H_sdiff, H_cadv) | Validates our d-way; adaptive part inapplicable |
| 7 | Set branching via value clustering; cites Larrosa's crossword result | Candidate F6 (exact, non-clustered form) |
| 8.2 | Future work (global constraints, propagator-oriented solvers, tuning) | Context only |

In the markdown extraction, §3.4 sits around lines 604–720, Ch 5 at
1157–1443, Ch 6 at 1443–2048, Ch 7 at 2050–2197, §8.2 at 2216+.

## What the thesis validates about the current fill core (no action)

- **d-way branching is the right call with dom/wdeg.** Ch 6's central
  finding (1600 instances): once the VOH is adaptive, full 2-way branching is
  often *worse* than d-way on binary CSPs, contrary to folk wisdom. The
  thesis's adaptive-branching heuristics exist only to repair 2-way's
  right-branch errors — inapplicable to a d-way solver. Our fixed d-way
  scheme (`mac_search/10`) is thesis-endorsed.
- **Our restart shape matches theirs** (geometric, factor 1.5), and our
  inter-attempt ×0.99 aging (`mac_age_weights/1`) is already a mild member of
  the aging family the thesis explores.

## Ranked fill candidates

All candidates preserve the determinism/pinning contract (they are
deterministic given the pinned PRNG constants) and are contained changes to
the `mac_*` layer in `prolog/crosswordsmith/fill.pl`. Suggested arbitration:
one ratcheted experiment each on the fill benchmark ladder, in this order.

### F1. Weight-ordered revision queue (Ch 5) — effort S, best bet

Current state: `mac_propagate/7`'s worklist is a plain list with guarded
front-push — **no revision-ordering heuristic at all**. The thesis's Ch 5
proposes ordering the variable-oriented revision list by statistics we
already maintain: pop the queued slot with the smallest dom/wdeg
(`v_dom/wdeg`) or largest wdeg (`v_wdeg`), and inside a revised slot process
its edges in descending weight order. Rationale ("ASAP principle"): reach
wipeouts earlier, which both cuts propagation work *and* — because weights
feed variable selection — steers dom/wdeg toward better-informed weights
(the interaction is the thesis's own observation, Examples 2/5/7). Reported
swings up to ~5× cpu between revision orderings (Table 5.4), with node
counts cut too, not just constraint checks. Everything needed is in hand:
per-edge weights in `cw_mac_weights`, domain popcounts. Cost risk: a min-scan
per pop on a typically short worklist; cache domain sizes if it shows.

### F2. Deletion-responsible weight credit — H1/H2/H3 and alldel (§3.4.2, §3.4.4) — effort S/M

Current state: `mac_bump_edge/1` fires only at the two wipeout sites
(`mac_place/9`, `mac_revise/11`) — the edge that lands the final blow gets
all the credit. The thesis's Example 1 applies verbatim to our mask-AND
propagation: when several crossing edges shrink a slot's mask and the last
one zeroes it, the edge that removed 4/5 of the candidates gets nothing.
Fix: record the responsible edge per deletion; on wipeout bump every
recorded edge by 1 (H1), by its deletion count (H2), or normalized by domain
size (H3). Cheap in our representation — each `mac_revise` already knows the
popcount delta, so a small per-node log of `(EdgeId, Deleted)` suffices. The
`alldel` variant (credit on *every* deletion, not just DWO-terminated
episodes) is one more arm of the same experiment. Thesis evidence: no
dominance, but plain dom/wdeg never won; needs the ladder to arbitrate.

### F3. Intra-attempt weight aging (§3.4.5) — effort S

Current state: aging only *between* restart attempts (×0.99). The thesis
(following Chaff/BerkMin) ages *during* search — divide all weights by 2
every ~20 backtracks — so recent conflicts dominate. Relevant because late
attempts run thousands of nodes under a grown cap where stale early weights
accumulate influence. Their aged dom/wdeg was the best performer on several
structured classes (5× on All-Interval Series) and the *worst* on Golomb
Ruler — a pure benchmark gamble, cheap to try.

### F4. "Fully assigned" / potential-DWO weights (§3.4.6) — effort S

On a wipeout, additionally bump the edges that participated in *fruitful*
(value-deleting) revisions in the same propagation episode — they mark
"potential DWO" slots that a different revision order would have wiped out.
Distinctive property (§5.5 variance study): alldel and fully-assigned are far
**less sensitive to revision ordering** than plain dom/wdeg. That makes F4 a
hedge/companion for F1 — either order the queue (F1) or make the weights
order-robust (F2/F4) — and the combinations are worth testing together.

### F5. Pinned random-probing weight initialization (§3.4.2, after Grimes & Wallace) — effort S/M

Current weakness this targets: attempt 1 is strict greedy with *uniform*
weights, so dom/wdeg degenerates to domain-over-degree during exactly the
cheapest nodes of the run. Fix: a handful of very short, low-cap probe
searches (pinned PRNG stream, so the determinism contract holds) whose only
output is warmed conflict weights for the real attempt 1. Our restart loop
already approximates this by accident (weights persist across attempts); F5
makes it deliberate and moves it *before* the greedy attempt.

### F6. Set branching on crossing-letter projections (Ch 7 + Larrosa) — effort M/L, the envelope candidate

The thesis's generic method (cluster values by VOH score via x-means) is the
wrong tool for us, because crosswords admit an **exact** equivalence the
clustering only approximates: two candidate words that agree on all
**checked (crossing) cells** constrain neighbours identically — they are
propagation-interchangeable, differing only in unchecked cells. Branch on
the projection "letters at crossing positions", trying the highest-scored
representative per class; a failure refutes the whole class soundly. On
fully-checked open grids classes are near-singletons (little gain); on
UK-style blocked grids (~half the cells unchecked) classes are large — and
the two acknowledged unclosed rungs, `blocked_13b` and `blocked_15a` (which
defeat ingrid too), are exactly that shape. Wrinkles to design around: the
all-different filter and final word scoring *do* distinguish class members —
resolve both at solution-assembly time. This is the only candidate aimed at
the feasibility envelope rather than latency, and it echoes the thesis's own
Ch 7 citation of Larrosa, whose subproblem-merging was validated on
crossword generation.

## Skipped, with reasons

- **Adaptive branching (H_sdiff, H_cadv, §6.4)** — repairs 2-way branching's
  right-branch errors; meaningless for a d-way solver.
- **Impact-based heuristics (§3.3, §4.3)** — the thesis's own evaluation:
  they lose to conflict-driven heuristics on cpu time, worst with large
  domains. Our domains are thousands of candidate words.
- **Exhaustive multi-DWO detection ("freeze/undo", early §3.4.6)** — the
  thesis's own negative result: only ~0.2% of revisions are DWO revisions; a
  second DWO in one revision list is essentially never found.
- **2-way branching, dichotomic domain splitting, clustered set branching in
  the generic x-means form** — dominated for our problem by the exact
  interchangeability of F6.

## Why this is fill-only: the arrange verdict

Arrange (Flavour A) fails both premises of the framing result, and — unlike
when this question was purely conceptual — we now have our own measurements
saying so:

1. **No fixed constraint network.** Arrange's variables would be words, but
   which words cross which, and where, is *created* by each placement
   decision. There is no stable constraint object for a weight to attach to,
   no AC queue to reorder, no fixed domain to set-branch over, and restarts
   have no constraint-indexed state to carry.
2. **Failure-starved regime.** The thesis's machinery learns from wipeouts.
   The arrange campaign's P1 probe (`docs/experiments.md`) measured the
   dense ladder: backtracking is shallow (84–86% of retreats at distance
   1–2) and the hardest rung (21×21, 80 words) backtracks **zero** times —
   the cost is per-node counting, not tree shape. A conflict-driven
   heuristic has nothing to eat.
3. **The word-level analogy was already tried.** E-H6 implemented a
   dom/wdeg-style failure weight as an MRV tie-break and **rejected it on
   measurement**: the failure concentration is real and the weights land on
   the right words, but the low-count MRV buckets are mostly forced — there
   are no ties to break, and floating a blamed word just picks a different
   equally-doomed branch (`docs/research/arrange-search-algorithms.md`).
4. **The envelope lever is already known and isn't a heuristic.** P2/P3
   established the decoupling result: near-cliff arrange instances are
   search-bound and node-cost-insensitive; only a different tree per attempt
   (budget-triggered randomized restarts, the planned rescue path) moves the
   feasibility envelope. The thesis adds nothing to that beyond what the
   Gomes/Selman/Kautz line already contributed.

`lint` and `export` involve no search; nothing here applies to them.

## The conceptual case: arrange as a CSP

Question examined (2026-07-16 session): is there a case for redesigning
arrange on CSP-modelling lines, conceptually rather than pragmatically?
Answer: the *model* is clean and the exercise is legitimate; the thesis's
*search machinery* still mostly loses its grip; and there is one genuinely
interesting synthesis.

**The natural model.** Variables = the words (closed set). Domains = legal
placements (row, col, direction) — a few hundred values per word on a 17×17.
Binary constraints = pairwise placement compatibility (non-interacting, or
crossing with matching letter and no adjacency violation). Plus one global
constraint (single interlocked component) and the objective (the capped
checked-cells reward). This is a textbook binary COP, and formulated this
way the DP-8 machinery runs on it unchanged; crosswordsmith would become one
adaptive solver core with two models. Adjacent precedent: Botea's
hierarchical re-encoding (see `arrange-search-algorithms.md` #5, "last
resort"). Note the standing framing result from that doc still holds — our
free-canvas formulation has no direct published precedent; published
"crossword CSP" work is dictionary-fill, our fill's shape.

**Why the machinery still starves.** Feasibility in this model is nearly
trivial — almost any pair of placements is compatible because words can
avoid each other. The difficulty lives in the objective: arrange is
essentially the Crozzle problem (maximize interlock of words placed from a
closed list), optimization-hard *because* the constraint landscape is loose.
Loose constraints ⇒ rare wipeouts ⇒ flat weights ⇒ dom/wdeg degenerates to
domain-over-degree. Arrange-as-CSP sits in exactly the corner of CSP space
(loosely constrained, objective-driven, massively symmetric) the thesis
never visits. Compounding gaps, all explicitly outside the thesis: canvas
translation/reflection symmetry (E-H1 took only the transpose pair; a
placement CSP needs systematic symmetry breaking), connectivity as a global
constraint (the thesis's own future-work concession), and best-effort drops
(soft constraints / MaxCSP).

**The synthesis worth keeping: optimization as iterated satisfaction.** Pose
"does a layout with objective ≥ K exist?" and walk K (binary search or
destructive descent). Each K-feasibility instance is a *tight* CSP — a high
interlock demand makes crossings mandatory, pairwise constraints prune hard,
wipeouts become common — and the thesis's machinery is back in its natural
habitat. The terminal "no solution at K+1" run is an optimality proof: a
statement the greedy constructor can never make ("best-within-budget, not
proven-optimal" is the current contract, README). The engines' answers would
differ in epistemic status, not just quality.

**Category note.** The thesis, despite its title, is about adaptive *search*
over a model someone already chose; it says nothing about viewpoint choice,
channeling, symmetry, or soft constraints — where most of the intellectual
content of an arrange redesign would live (wider CP literature, e.g. the
Rossi/van Beek/Walsh handbook tradition).

**Verdict.** Worth a bounded design-note exercise (model + symmetry-breaking
scheme + K-iteration semantics, perhaps a toy prototype on the bundled 17
clues to locate the K phase transition) for the conceptual grounding:
arrange = Crozzle = loose COP, fill = tight CSP, one solver core + two
models, iterated satisfaction as the bridge. **Not** a shipping-engine
replacement: the anytime best-effort greedy is a product feature, the
deterministic default is a product requirement, and P2/P3 already identified
seeded restarts — not search power — as arrange's envelope lever.

## Recommended sequence

1. F1 (revision ordering) behind a strategy seam, benched on the fill ladder.
2. F2 + F4 as arms of one weight-scheme experiment (H1/H2/H3/alldel/fully
   assigned), combined with and without F1.
3. F3 (intra-attempt aging) and F5 (probing init) as cheap follow-on arms.
4. F6 (projection set branching) as a separate design spike targeted at
   `blocked_13b`/`blocked_15a` — the only envelope bet in this note.
5. The arrange-as-CSP design note only if/when the conceptual grounding is
   wanted; no benchmark pressure motivates it today.
