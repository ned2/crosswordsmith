# Research: search-algorithm candidates for the arrange engine

Distilled from a literature sweep (2026-07-04) commissioned for the arrange
performance campaign (see `docs/experiments.md`, "Crosswordsmith arrange
campaign"). Scope: techniques that could shrink the mrv_inc SEARCH TREE on
dense instances (the ladder's hard tail), where constant-factor work only
shaves linearly. Constraints honoured throughout: pure SWI-Prolog (must run
under the WASM build), deterministic default search, exact satisfaction
(place ALL words; reward rescoring stays outside the search).

## The load-bearing framing result

**Our formulation has no direct published precedent.** Nearly everything
published as "crossword construction" (Ginsberg's Dr.Fill lineage, Meehan &
Gray 1997, Botea 2007, arXiv:2308.04688) is *dictionary-fill*: variables =
fixed grid slots, values = words from a large dictionary. Crosswordsmith's
arrange is the transpose: variables = a small fixed word list, values =
geometric placements on an open canvas (word-mesh / free-canvas interlock).
The closest analogue (arXiv:2007.04663, unconstrained puzzle generation)
still grows a pattern around a vocabulary. Everything below transfers from
generic CSP literature, not crossword-specific literature — treat payoffs as
bets to measure on the ladder, not published guarantees.

Also load-bearing: the pre-pivot negative results (experiments.md I2 value
ordering, I6 static degree) show this search is intolerant of per-node
overhead — any new signal must be near-free at the node level. But their
"mrv_inc barely backtracks" premise is FALSE on the current dense ladder
rungs (15x15_36w ≈ 91.8M inferences post-E-H1), so tree-size techniques are
back in play there.

## Ranked candidates

### 1. Adaptive weighted-degree tie-breaking (dom/wdeg-style) — TRY FIRST

Boussemart et al.'s dom/wdeg: weight constraints by how often they cause
failures *during this run*; order variables by domain/weight. Unlike I6's
static degree (structural, known ahead of time — gave zero tree reduction),
this is a failure-driven signal. Key insight for a pure-Prolog fit:
`find_intersecting_word`/`assign_word` already enumerate and reject candidate
crossings inside `mrv_count` — each rejection is a forward-detected
micro-conflict generated and consumed within the same node, so the weights
can ride the existing `state(...)` thread through `select_inc` (backtracking
restores them automatically; no assert/retract impurity). Use weights to
break MRV ties within the same capped-count bucket only (preserves fail-first).
Effort S/M. Real risk of a null result given I6; must be benched on the ladder.
- Explanation-based weighted degree: https://www.researchgate.net/publication/318137036
- Evaluation of modern variable-selection strategies: https://www.researchgate.net/publication/228635686

> **P1 UPDATE (2026-07-04, measured):** the instrumentation pass ran (see
> experiments.md P1). Verdicts: **CBJ is ruled out** — retreat distances are
> 84-86% at 1-2, max ever observed 7; there is nothing to jump over. **wdeg
> has strong signal** — wipeouts concentrate on ~8-10 words per dense rung
> (top-5 = 75-91%) — but its addressable surface is only the thrashing
> corner on the 15x15 dense rungs; 21x21_80w has ZERO backtracking (its cost
> is pure per-node counting; only constant-factor work helps it). Sections
> below kept for the reasoning trail.

### 2. CBJ / dynamic backtracking — RULED OUT by P1 (was: instrument first)

Prosser's CBJ backjumps to the deepest conflict-set variable; Ginsberg's
dynamic backtracking (JAIR 1993 — demonstrated on crosswords, but the
dictionary-fill kind) keeps polynomial space via discarded eliminating
explanations. The follow-up literature is explicit that CBJ's payoff
*shrinks as variable ordering and forward checking improve* (Bacchus & van
Run; "CBJ becomes useless" under FC + good dynamic ordering) — and mrv_inc is
already both. Before building anything: instrument backtrack-target distance
on the dense rungs (log how far back each backtrack's eventually-successful
alternative sits). That single number is CBJ's plausible ceiling.
Implementation route if justified: catch/throw as the jump primitive
(`throw(backjump(Level, ConflictSet))` caught by the matching recursion
frame) bolted onto `assign_words_inc/9` — effort M; a full explicit
frame-stack rewrite (L) only if that proves too coarse.
- Ginsberg 1993: https://www.jair.org/index.php/jair/article/view/10107
- Ginsberg et al. 1990 (crossword search lessons): https://cdn.aaai.org/AAAI/1990/AAAI90-032.pdf
- Prosser CBJ: https://www.dcs.gla.ac.uk/~pat/cpM/papers/prosser_cbj.pdf
- Diminishing-returns result: https://arxiv.org/pdf/1106.0254

### 3. Budget-triggered randomized restarts (+ cross-restart nogoods) — RESCUE PATH ONLY

Gomes/Selman/Kautz heavy-tailed-runtime result: randomized restarts kill the
tail on hard combinatorial instances. Fits crosswordsmith as an escalation on
`not_proven` (budget exhausted) ONLY — the deterministic default is a product
requirement. The `--seed`/`--shuffle` machinery (`set_search_seed/1`,
`seed_word_order/2`, `order_candidates/2`) already provides reproducible
perturbation; the new piece is spending a second budget on seeded restarts,
optionally with simple cross-restart nogoods ("these two words never legally
cross") — impurity acceptable on this explicitly non-deterministic path.
Helps the currently-infeasible dense tail, not the median. Effort M.
- Heavy-tailed phenomena: https://link.springer.com/article/10.1023/A:1006314320276
- Nogood recording from restarts: https://www.ijcai.org/Proceedings/07/Papers/019.pdf

### 4. Matching-based (Hall-style) global feasibility check — EXPLORATORY, overhead-risky

A bipartite-matching feasibility bound over remaining words vs available
crossing letters would prune earlier than "some word has 0 placements", in
the spirit of alldifferent's matching consistency (Régin). Structurally the
same kind of bet as I2/I6 (pay per-node for earlier pruning) and more
expensive than the degree scan that lost in I6 — only conceivably viable
invoked periodically/at shallow depth on dense instances. Effort M/L.
- alldifferent survey: https://arxiv.org/pdf/cs/0105015

### 5. Crossing-graph decomposition (backbone-first placement) — LAST RESORT

Botea's hierarchical CSP encoding reports big wins for grid composition, but
the gain comes from a re-*encoding*, not an ordering tweak, and I6 already
showed static crossing-graph connectivity is not predictive here. Revisit
only if instrumentation shows strong core/periphery structure in dense
instances' crossing graphs. Effort L.
- Botea: https://www.cse.cuhk.edu.hk/~jlee/cp07Model/pdf/crossword.pdf
- Meehan & Gray (heuristic word-by-word beat a CLP encoding — mild evidence
  for keeping the DFS shape): https://www.gtoal.com/wordgames/meehan/cross.pdf

## Ruled out

- **Limited discrepancy search**: LDS needs a value-ordering heuristic that
  is usually right; I2 found any value-ordering signal is a net loss here,
  and there is no ranked candidate order for LDS to trust.
  https://www.ijcai.org/Proceedings/95-1/Papers/080.pdf
- **Further symmetry breaking**: the transpose symmetry is already taken
  (E-H1). Remaining symmetries would mainly speed `--enumerate`, not the
  construction path that blows up.

## Recommended sequence for the campaign

1. Instrumentation pass on the dense rungs: backtrack-distance distribution
   (CBJ ceiling) + where failures cluster (does a wdeg-style weight signal
   land on the right culprits?). Near-zero cost, de-risks both bets.
2. #1 (wdeg tie-break) behind a strategy flag, benched on the full ladder.
3. #2 only if step 1 shows long backjump distances; #3 as a product feature
   for the not_proven tail independent of raw search speed.
