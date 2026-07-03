# Crossword Setting & Generation — Research and Algorithm Notes

A structured Q&A capturing a session on how crossword setters work, what makes a
puzzle *satisfying*, the tools that fill grids from a word set, and two concrete
CSP-style generation algorithms — one for American/themed grids and one purpose-built
for cryptic grids.

> **Status — research/exploration snapshot.** Where this disagrees with
> [`design-spec.md`](./design-spec.md) or [`cryptic-layout-spec.md`](./cryptic-layout-spec.md),
> those specs win. It is retained mainly for the concrete generation **algorithms**
> (Q4–Q6) — a road-not-taken: `cryptic-layout-spec.md` adopted only the greedy
> "Approach B" constructor and deliberately rejected the GAC / signature-table
> machinery, so that design space survives nowhere else. The setter-craft and
> tooling-landscape background (much of Q1–Q3) is covered more fully, and with
> primary sources, in [`cryptic-setter-research.md`](./cryptic-setter-research.md);
> those parts are trimmed to pointers here to avoid duplicating it.

---

## Q1. What process and methodologies do crossword setters use to make *satisfying* crosswords?

Crossword construction is really two crafts. "Satisfying" means something different in
each, but both rest on the same triangle: **constraint, fairness, freshness.**

### The American themed workflow
A fixed order, and much of the skill is in *not* doing steps out of sequence:
**theme → grid → fill → clues → polish.**

- **Theme first.** The puzzle's calling card. Good themes are new, narrowly defined, and
  consistently applied. Theme entries are long phrases placed in mirror-symmetric
  positions, so they must be checked for symmetry/consistency before anything is built.
  Editors often want to see *just the theme set* first.
- **Grid.** Theme entries locked into symmetric positions; black squares ("blocks") added
  to break white space into fillable regions while preserving 180° rotational symmetry
  and connectivity. A **cheater square** adds a black square that doesn't change the word
  count but splits a hard region into two easier ones.

### What separates good fill from merely valid fill
This is where "satisfying" is won or lost.

- The non-theme words are the **fill**. A brilliant theme can still get rejected for
  clunky fill.
- The enemy is **crosswordese** — obscure/archaic words (EPEE, ERNE, OLEO) that recur only
  because of convenient letter patterns. The goal (per Will Shortz) is lively, familiar
  vocabulary: crosswords should reflect life.
- Modern method = **constraint-solving with a scored dictionary.** Every word in a
  constructor's list carries a quality score (lively = high, crosswordese = low; scales are
  personal). Software (Crossword Compiler, CrossFire, Crosshare, Sympathy) autofills,
  preferring high-scoring words; the constructor curates and overrides.
- Expert tip: **start filling the most-constrained region**, not the least.
- Resist "Scrabbliness" (forcing J/Q/X/Z) when it hurts smoothness. **Smoothness beats show.**

### Clue-writing and misdirection
Clues are written **last**, and difficulty is mostly tuned here — the same grid can be a
Monday or a Saturday depending on cluing. Conventions are tight (clue and answer agree in
tense/part-of-speech/number; a trailing "?" signals wordplay). Late-week difficulty comes
from **fair misdirection**.

### The cryptic (British) tradition
Inverts the emphasis: the grid matters less, the artistry is in individual clues
(Ximenean *definition + fair wordplay + nothing else*, the Ximenean↔Libertarian
spectrum, surface reading, the "padding" fault).

> Covered in depth — with the full setter pipeline and primary sources — in
> [`cryptic-setter-research.md`](./cryptic-setter-research.md) ("How cryptic
> setters actually work") and design principle #7 there ("Be Ximenean-configurable,
> never Ximenean-only"). Not repeated here.

### The psychology of "satisfying" (both traditions)
**Fair misdirection** (tricked, but fairly), the **aha moment** (the click of resolution),
**calibrated difficulty** (resist without frustrating), and **freshness over filler**.

### Primary sources
- Patrick Berry, *Crossword Constructor's Handbook* (American).
- Ximenes, *On the Art of the Crossword* (1966) (cryptic).
- The NYT five-part "How to Make a Crossword Puzzle" series.

---

## Q2. How are *cryptic* crosswords laid out — the grid construction process specifically?

Punchline: for standard cryptics the grid is usually the **easiest** part and often
isn't designed from scratch — the leverage is **half-checking** (only ~half of each
entry's letters are checked; the rest are "**unches**" belonging to one word, giving
the lattice / "odd-odd" look). That decoupling is what Q5's algorithm exploits, and
it's why setters spend their effort on clues, not grid-filling. The genuinely hard
layout work lives in the barred-grid variants and in hidden structure (Ninas).

> The full treatment — checking/unch fairness rules, stock-grid reuse, barred grids,
> Ninas, the Anax "few clues during fill / wordplay-rationing" workflow, and the
> observed grid-size distribution — lives in
> [`cryptic-setter-research.md`](./cryptic-setter-research.md) (checking rules, the
> "Stock-grid libraries" and "Grid sizes in the wild" addenda) and is encoded as
> normative lint/profile rules in [`design-spec.md`](./design-spec.md) §8. Not
> repeated here; Q5 below re-defines the checked/unch terms it needs independently.

---

## Q3. Are there tools that go the other way — given a target word set, fill them into a satisfying crossword *if possible*?

The answer splits by how literally you mean "given a target set."

### Easy version vs. hard version
- **Easy (mature, fast):** design/pick a grid, **pin a few words you want**, and let software
  fill the rest from a large scored dictionary. This is the standard autofill every mainstream
  tool does — see the tool table in
  [`cryptic-setter-research.md`](./cryptic-setter-research.md) (Crossword Compiler, Qxw,
  Crossfire, Crosshare, Exet, …) rather than re-listing them here.
- **Hard (the real inverse):** given a **closed set**, pack as many as possible into one interlocked
  grid, designing the grid around them with little/no outside filler. *(This closed-set inversion —
  no filler dictionary — is exactly crosswordsmith's model; see the "closed-set vs. authentic-layout
  tension" analysis in `cryptic-setter-research.md`.)*

### Tools aimed at the inverse problem
- **AmuseLabs "Magic Fill" (PuzzleMe):** built for the small-custom-list case — upload your list,
  it produces dense interlocked layouts, and supports **"Must Have Words."** It's candid about
  failure: if the grid won't fill, *relax* (bigger wordlist, lower min length, more black squares,
  larger grid).
- **Education-oriented packers** (Crossword Labs, EclipseCrossword, Armored Penguin): take *your*
  list and prioritize including all of it — but produce the sparse "criss-cross" style, not a dense
  newspaper grid. That's the trade for forcing a whole closed set.
- **Open-source CSP generators** (many on GitHub): take a grid + word list, fill via constraint
  propagation + backtracking, or report failure.

### Why "if it's possible" is the crux
- **The Scrabble effect:** a dense interlock is hard with a small set, especially American-style
  (every letter in two words). A small arbitrary set may lack enough shared letters; you get
  Scrabble-like sparseness instead.
- **Hardness / phase transition:** grid-filling is NP-complete in general and shows a classic CSP
  phase transition — a density threshold below which solutions are plentiful and above which they
  almost surely don't exist, with search hardest at the boundary. The more you force and the denser
  you demand, the closer you push to the wall.
- Tools therefore behave as **negotiation engines**: they *report* infeasibility (e.g. flagged
  unfillable slots) and *let you relax toward feasibility* (filler words, bigger grid, more black
  squares). The reliably satisfying recipe is the constructor's method in reverse: force only a few
  must-haves and let a high-quality scored dictionary do the rest.

---

## Q4. A concrete CSP-style algorithm for the inverse problem (American / general grids)

The thing that makes this harder than textbook fill: the **grid isn't given** — you co-design the
black-square pattern *and* the word placement. So it's a two-level search: an outer search over
geometry wrapped around an inner CSP over word assignment.

### 1. Frame as constrained optimization (not pure satisfaction)
- **Inputs:** target multiset `W` (each with priority `pri(w)`), scored filler dictionary `D`
  (`quality(w) ∈ [0,1]`), constraint budget `B` (size range, style, symmetry, `minLen`, max black
  ratio, must-include-all flag).
- **Variables:** the grid's slots `S` (maximal white runs).
- **Domains:** `dom(s) = { w ∈ (W ∪ D) : len(w) = len(s) }`.
- **Hard constraints:** crossing-letter agreement; all-different over words; `len ≥ minLen`;
  connectivity; symmetry; and (if must-include-all) every `w ∈ W` occupies a slot.
- **Objective `Q`** (maximize): aesthetic/quality sum (§5).

### 2. Inner CSP — filling a fixed geometry well
Precompute a position–letter index:

```
index[L][pos][letter] = bitset of word-ids of length L with `letter` at `pos`
```

Represent each `dom(s)` as a bitset; the key quantity is cheap:

```
feasibleLetters(s, pos) = { ℓ : index[len(s)][pos][ℓ] ∩ dom(s) ≠ ∅ }
```

Arc consistency (the right propagator):

```
revise(s, t):                      # s,t cross at (i,j)
    F = feasibleLetters(t, j)
    removed = false
    for w in dom(s):
        if w[i] ∉ F: dom(s).remove(w); removed = true
    return removed

AC3(queue of crossing-arcs):
    while queue:
        (s,t) = queue.pop()
        if revise(s,t):
            if dom(s) empty: return WIPEOUT
            for u in crossings(s) \ {t}: queue.push((u,s))
    return OK
```

Search = backtracking + MRV + quality-weighted value ordering + MAC + restarts:

```
search(assign):
    if complete(assign): return assign
    s = selectSlot()                       # MRV; tie-break most unfilled crossings; prefer must-target
    for w in orderValues(s):               # by valueScore desc + random jitter
        if w ∉ used and consistent(w,s):
            push domains; assign[s]=w; used.add(w)
        if AC3(arcs_into(s)) != WIPEOUT:
            r = search(assign); if r≠FAIL: return r
        pop domains; unassign s; used.remove(w)
    return FAIL
```

```
valueScore(w,s) = α·quality(w)
                + β·targetBonus(w,s)        # large if w∈W, scaled by pri(w)
                + γ·Σ_{(t,j,i) crossing s} log(1 + |{u∈dom(t): u[j]=w[i]}|)   # least-constraining
```

- **Restarts:** cap each attempt at a node budget; restart with fresh jitter to escape thrashing
  near the phase transition (what shipping fillers do when you "reshuffle").
- **Anytime branch-and-bound:** keep `bestQ`; prune when optimistic bound ≤ `bestQ`; keep the best
  fill, not the first.

### 3. Outer loop, Approach A — template + seed-matching (robust default)
1. Library of legal symmetric skeletons per size.
2. Cheap feasibility pre-filter (Hall's condition): bipartite-match target lengths to slot lengths;
   if no matching covers `W`, discard the template instantly.
3. Pin targets into compatible (and, for theme entries, symmetric) slots; fill the rest via §2.
4. Score each completed grid by `Q`; keep the best.

### 4. Outer loop, Approach B — constructive co-design (build the grid *around* the words)
Grow a layout on an unbounded lattice, then regularize.

```
expand(state):
    w = pickNextTarget(state)               # highest pri among placeable / most-anchorable
    placements = []
    for each occupied cell (r,c,ℓ):
        for each k with w[k]==ℓ:
            dir' = perpendicular(owner_dir(r,c))
            (r0,c0) = origin so w[k] lands on (r,c)
            if legal(w,r0,c0,dir',state): placements.add((r0,c0,dir'))
    rank by (newCrossings desc, bboxGrowth asc)
    return placements

legal(w, r0,c0, dir, state):
    - every cell matches an existing letter or is empty
    - terminator cells before/after w are empty/off-grid          # no accidental extension
    - each letter cell's flanking-perpendicular cells empty/off-grid
      UNLESS that cell is a real crossing                          # forbids parallel-adjacency junk
    - bounding box ≤ size budget
```

Then **regularize and densify**: arbitrary packing breaks 180° symmetry — either drop to British/
barred conventions, treat symmetry as a soft penalty, or fall back to Approach A. Uncrossed white
runs become new empty slots → hand back to §2 for dense dictionary completion.

### 5. The objective `Q` (where "satisfying" lives)

```
Q = Σ_placed quality(w)
  + λ_T · Σ_{w∈W placed} pri(w)
  + λ_C · checkedRatio          - μ_S · shortGlut
  + λ_O · openness              - μ_X · crosswordeseCount
  + σ   · symmetryBonus         - μ_D · duplicationPenalty
  + τ   · targetPlacementAesthetics  - μ_A · abbrevAndPartialCount
```

Tuning the weights *is* the tool's editorial taste.

### 6. "If it's possible" — feasibility + relaxation ladder
- **Detect:** Hall pre-filter catches length-infeasibility; if N randomized restarts all wipe out,
  declare *likely infeasible at this density*.
- **Relax (least-damaging first):** (1) allow dictionary filler; (2) demote low-pri targets;
  (3) enlarge grid; (4) raise black-square budget / cheater squares; (5) lower `minLen` / allow
  unches; (6) relax symmetry; (7) drop the hardest target, reporting which and why.

### 7. Top-level

```
solve(W, D, B):
    targets = scoreAndSort(W)
    geometries = (B.style==american) ? templateLibrary(B) : [constructiveSeeds(targets,B)]
    best = ∅
    for g in geometries:
        if not canHost(g, requiredTargets): continue
        pinTargets(g, targets)
        fill = fillCSP(g, D, B)            # §2, anytime, keeps best Q
        if fill: best = argmaxQ(best, fill)
    return best ? best : relaxAndRetry(W, D, B)
```

**Why this recipe:** bitset/feasible-letter AC makes propagation affordable on huge domains;
MRV + most-constraining attacks tight corners first; quality-weighted value ordering surfaces
high-`Q` solutions early; random restarts escape the SAT/UNSAT threshold.

---

## Q5. The same exercise, purpose-built for the *cryptic* grid (self-contained)

Built from the cryptic grid's own structure; defines its terms independently of Q4.

### What a cryptic grid is (the facts that drive the design)
Blocked, odd-sized (15×15 typical), 180°-symmetric, lights in a lattice. Defining property:
**half-checking.** Every white cell is one of two kinds:
- **Checked** — in one across light *and* one down light. The only place two lights interact.
- **Unchecked ("unch")** — in exactly one light. A **degree-1 free variable**; touches nothing else.

That decoupling is the whole leverage point. Fully-checked grids couple every cell; here the
unches let us split the problem into a small coupled core + an almost-independent completion.

### Cell model + fairness constraints on geometry

```
classify(grid):                       # minLen typically 4 (3 grudgingly)
    for white cell x:
        inAcross = horizontalRun(x).len >= minLen
        inDown   = verticalRun(x).len   >= minLen
        x.kind = CHECKED   if inAcross and inDown
               = UNCHECKED if inAcross xor inDown
               = ILLEGAL   otherwise            # isolated cell → reject grid

fair(grid):
    for each light ℓ:
        if len(ℓ) < minLen: return false
        if checkedCount(ℓ) < ceil(len(ℓ)/2): return false   # ≥ half, round up
        if two consecutive cells of ℓ are both UNCHECKED: return false
    return symmetric180(grid) and connected(white(grid))
```

"≥ half, round up" forces odd-length lights to begin and end checked.

### Key reframing: a CSP over *only the checked cells*
Don't search over words — search over **letters at checked cells**, treating each light as a
*pattern* its checked letters must complete.

For a light ℓ with checked positions `p₁<…<p_k`, a word's **checked signature** is
`sig_ℓ(w) = (w[p₁],…,w[p_k])`. Precompute realizable signatures:

```
sigTable(length, checkedMask):                 # cache by (length, mask) — masks repeat a lot
    pos = indices set in checkedMask
    return { (w[p] for p in pos) : w in Lex[length] }    # store as set / trie / MDD
```

**The skeleton CSP:**
- **Variables:** checked cells. **Domain:** 26 letters (tiny).
- **Constraint, one per light:** the letters on ℓ's checked cells (in order) must form a tuple in
  `sigTable(ℓ)` — a table (allowed-tuples) constraint.
- **Structural gift:** each checked cell is in exactly one across + one down light, so **every
  variable appears in exactly two constraints.** The network is extremely sparse.

```
GAC(csp):                              # simple tabular reduction + projection
    repeat until stable:
        for each light ℓ:
            ℓ.table = { t in ℓ.table : every letter of t still in its cell's domain }
            if ℓ.table empty: return FAIL(ℓ)          # ℓ is the failure locus
            for each cell c in scope(ℓ):
                domain(c) ∩= { letters appearing at c's position across ℓ.table }
                if domain(c) empty: return FAIL(c)
    return OK

solveSkeleton(csp):
    if GAC(csp) == FAIL: return FAIL
    if every checked cell is singleton: return assignment
    c = unfixed cell with smallest domain, tie-break most-constraining crossing
    for letter in domain(c) ordered by kindness(letter,c):
        push; domain(c) = {letter}
        r = solveSkeleton(csp); if r != FAIL: return r
        pop
    return FAIL
```

### Realization — the independent completion phase
Payoff property: **skeleton-consistent ⇒ a word fill exists** (each light's checked tuple is, by
construction, in its `sigTable`). Because unches are degree-1, lights choose words **independently**:

```
realize(grid, skeleton):
    if theme.nina: pin message letters onto chosen perimeter unch cells
    used = {}
    for each light ℓ:                     # Nina-touched lights first
        pat = checked letters from skeleton + pinned Nina unches, '?' elsewhere
        cands = { w in Lex[len(ℓ)] : matches(w,pat), w not in used }
        w = argmax_{cands} entryQuality(w)    # liveliness / clue-ability
        place w; used.add(w)
    return grid
```

Without a Nina, realization essentially never fails or backtracks (half the letters are free →
many completions). The only thing that can break it is over-constraining the unches — which is
exactly what a Nina does.

### Ninas and themes as first-class layout objects
A Nina fixes specific perimeter unch cells *before* realization. Since those cells are degree-1,
fixing them doesn't perturb the skeleton — it only narrows each affected light's candidate set.
The unch budget is spare capacity; a Nina is how you spend it.

### Seeding target words
Assign each target to a length-matching light (thematic targets in symmetric positions), then
collapse that light's table to the single signature `sig_ℓ(target)` — a unary reduction GAC
propagates outward.

### Geometry generation
- **Blocked / lattice:** primarily *select* from a library of fairness-validated symmetric grids,
  filtered so light-lengths cover the targets. For novelty, generate constructively: lay a symmetric
  black sublattice, run `classify`+`fair`, hill-climb mutations toward `fair`.
- **Barred (Azed/Listener):** no black squares; place bars so every run is a valid light and the
  stricter, higher checking ratio holds. Fewer unches → harder fill, but licenses obscure vocabulary.

### Objective (cryptic-specific)

```
Q = + w1 · entryLiveliness        - w6 · duplicationPenalty
    + w2 · unchSpread              - w7 · obscurityInChecked   # don't hide rare letters at crossings
    + w3 · lengthVariety
    + w4 · checkingKindness        # crossings yield helpful, common letters
    + w5 · ninaRealized / themeFit
```

### "If it's possible" — relaxation ladder (cryptic-native, cheapest first)
1. **Swap grid template** (a whole library exists). 2. **Abandon/shorten the Nina.**
3. **Slide checking toward the minimum.** 4. **Allow limited 3-letter lights.**
5. **Demote low-priority targets.** 6. **Switch blocked ↔ barred** (bars let you engineer where
unches fall). 7. **Go jumbo (23×23).** GAC names the exact failing light/cell each time.

### Why this is the natural cryptic algorithm
Generate a fairness-valid lattice → solve a small table-CSP over just the checked cells (26-letter
domains, every variable in two constraints) → complete each light independently → spend leftover
unch freedom on a Nina/theme. Every layer is enabled by half-checking; the existence guarantee
falls out of defining constraints over *realizable signatures* rather than whole words. None of
this layering exists when every cell is checked.

---

## Q6. Does the cryptic algorithm use any indexes?

Yes — but lighter than the American sketch, and concentrated in **setup and realization**, not the
propagation hot loop. (In Q4, a position-letter index sat *inside* propagation because every cell
coupled two words. Here the coupled core runs on 26-letter cell domains with table constraints, so
the word-level index moves out of the hot loop.)

Three indexes in play:
- **Lexicon bucketed by length** (`Lex[length]`) — used to build signature tables and drive
  realization. Trivial, but it's an index.
- **Signature-table cache keyed by `(length, checkedMask)`** — the cryptic-specific one. A
  precomputed projection (length + checked mask → realizable letter-tuples), heavily amortized
  because masks repeat across many lights.
- **Pattern-match index for realization** — a per-length trie or position-letter postings, so
  "find length-L words matching these fixed checked letters with unches free" is a lookup, not a scan.

What *doesn't* need an index: the **skeleton CSP solver itself.** With 26-letter domains and table
constraints, simple tabular reduction (STR) just rescans surviving tuples and projects to letters —
no support index required (you *can* add GAC-schema support indexes as an optimization). Also,
storing each signature table as a **trie/MDD** lets that structure double as the index for GAC's
projection step, so the table cache and the realization index can be the same object.

**Net:** the cryptic version is *less* index-dependent in its core search than the fully-checked
version; the indexing it does use is mostly amortized precompute plus an optional pattern trie for
the decoupled fill.

---

## Sources & further reading

For the **setter-craft and tooling** background (Berry, Ximenes, NYT, Anax, Crossword
Unclued, Alberich, Crossword Compiler / Qxw docs, …), see the fuller primary-source
bibliography in [`cryptic-setter-research.md`](./cryptic-setter-research.md); not
duplicated here.

Specific to this doc's unique material:
- AmuseLabs PuzzleMe "Magic Fill" docs (word-list-driven dense fill, Must-Have Words).
- **Algorithms (Q4–Q6):** CSP arc consistency, MRV / least-constraining heuristics,
  table-constraint GAC/STR, and the CSP phase-transition literature (crossword phase
  transition: Anbulagan & Botea).
