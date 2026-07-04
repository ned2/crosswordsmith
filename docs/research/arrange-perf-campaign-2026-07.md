# The arrange performance campaign (2026-07): what worked and why

A distillation of the six interventions that shipped during the July 2026
arrange performance campaign — the intuition behind each, and what it bought.
This is the "greatest hits" companion to the full experiment ledger in
`docs/experiments.md` (which also records everything that *didn't* work, with
evidence); the sibling research docs
(`arrange-search-algorithms.md`, `swi-vm-wasm-performance.md`) hold the
literature and VM cost-model groundwork.

**Headline:** search inferences on the benchmark ladder dropped **71–84% per
rung** with byte-identical output throughout. The hardest rung (15x15, 36
words) went from 182.6M inferences (~9.5s) to 38.3M (~2.2s). Inference counts
are machine-independent, so roughly the same relative speedup carries to the
WASM browser target.

## Context: the shape of the problem

`arrange --strict` places a small word list on an open canvas by DFS with MRV
ordering: at each node, every unplaced word's legal-placement count is
maintained (capped at 2 — the buckets 0/1/≥2 are all the ordering needs), the
scarcest word is placed next, and counts are incrementally refreshed for words
sharing a letter with the last placement. Determinism is a product
requirement: the same input must produce the same layout, byte for byte.

Two campaign-defining measurements framed everything:

- **P1 (instrumentation probe):** on dense rungs, backtracking is shallow
  (84–86% of retreats are distance 1–2) and one hard rung (21x21, 80 words)
  backtracks *zero* times. The cost is not a bad search tree — it is the
  per-node counting work itself. This killed the fancy tree-search ideas
  (backjumping, restarts-for-speed) before they were built, and aimed the
  campaign at per-node cost.
- **P2/P3 (envelope probes):** near-cliff instances are search-bound — they
  fail at 4x budget, and even a 50% per-node cost cut moved the feasibility
  envelope by zero words. **Node cost buys latency; only a different search
  tree per attempt (randomized restarts, a future product feature) moves the
  envelope.** Every win below is a latency/WASM-headroom win.

## The six shipped interventions

### 1. E-H1 — transpose-corner dedup (~2x, all rungs)

The strict search ran the same DFS from four start corners and kept the best
rescored result. But two of the corners produce exact transposes of the other
two: same crossings, same reward, mirrored geometry. **Intuition: never
search what you can derive by symmetry.** Run two corners instead of four.
The cheapest kind of win — half the work for a few lines of code — and a
reminder to look for symmetry *before* optimizing anything inside the search.

### 2. E-H2 — the grid as a term of unbound variables (−15..24%)

The grid was an AVL assoc (position → letter): O(log N) compares per read,
and every write rebuilt a root-to-leaf path (~15% of self-time went to the
AVL machinery). Replaced with a single compound term `grid(C1..C(N²))` whose
arguments are unbound variables: an empty cell *is* an unbound var, a letter
is a bound atom, reads are `arg/3` (O(1), allocation-free), writes are plain
unification, and undo-on-backtrack is the WAM's trail — for free.
**Intuition: in Prolog, the fastest possible "mutable store with undo" is
the one the VM already implements — unification against a big term.** Don't
build backtracking bookkeeping on top of a language whose runtime *is*
backtracking bookkeeping.

### 3. E-H3 — tabled crossing/letter memos (−8..31%)

For every candidate pairing of words, the search recomputed "at which letter
positions do these two words cross?" — millions of times across recounts,
plus a hidden tax: the old path went through `list_to_set/2`, whose
`must_be(list, ...)` type-check was the single biggest profiler line
(originally misattributed to `aggregate_all`; the profile's inclusive-time
frames lied about who was spending it). Tabled `pair_crossings/3` and
`answer_letters/2` — pure functions of the word list — computed once per
pair per search, abolished at the search boundary so counts stay
run-order-independent. **Intuition: inside an exponential search, anything
that is a pure function of the *input* rather than the *search state* should
be computed exactly once.** Corollary from the misattribution: verify which
predicate actually spends the profile time before optimizing the one at the
top of the listing.

### 4. E-H5 — a saturating counter (marginal, kept)

Placement counting used `aggregate_all(count, capped(2, ...))`. Replaced
with a hand-rolled failure-driven `nb_setarg` counter that stops at the
second solution. **Intuition: counting to a cap of 2 needs no general
aggregation machinery.** Honestly recorded as marginal on the ladder — it
shipped because it composes cleanly and simplifies the hot loop, not because
it moved the needle. Included here for completeness and calibration: not
every accepted change is a headline.

### 5. E-H7 — placed words as compact records (−1.9..−9.7%)

Placed words were SWI dicts; field reads cost ~1.5x compound-term access,
and the word's end cell was recomputed with `last/2` per candidate test.
Replaced with a `pw/8` compound (answer, letters, cells, direction, length,
start, *precomputed end*, clue number) behind accessor predicates.
**Intuition: for a fixed-shape record read millions of times, pay for field
access at construction, not per read** — and precompute derived fields (the
end cell) once at placement. Notably, this measured *better* composed on top
of the earlier wins than standalone: as the surrounding costs shrank, the
dict reads became a larger share of what remained. Constant-factor work
compounds — a reason not to dismiss "small" wins while bigger ones exist.

### 6. E-H9 — check-only counting + O(1) boundary test (−13..50%, the biggest single win)

Two composed changes to the innermost loop:

- **Don't materialize what you'll discard.** The counting path ("how many
  placements does word W have left?") ran the full placement machinery per
  candidate — binding grid cells, building the cell list, constructing the
  pw record — only for the counter to throw it all away. A non-binding twin
  (`check_word_fits`) walks the cells doing the same legality checks in the
  same order and builds nothing. Sound because nothing a single check would
  bind is re-read by that same check.
- **Replace "scan everything placed" with "read the cell you care about".**
  The no-merge rule (a new word must not touch an existing word end-to-end)
  was enforced by scanning *all* placed words per candidate — O(placed
  words), executed millions of times, and the reason the 80-word rung was
  expensive. Now a second grid-sized term of unbound vars mirrors the letter
  grid; placing a word marks its 0–2 boundary cells there, and the test is
  one `arg/3` + `nonvar` per cell. Undo is the trail again. The subtle trap:
  boundary cells must stay unbound in the *letter* grid, because adjacency
  checks must still read them as empty — hence a mirror grid, not a sentinel.

**Intuition: the grid itself is the spatial index.** Any rule phrased as
"does the new word conflict with *any* placed word?" can usually be
re-phrased as "is *this cell* marked?" — turning O(placed) into O(1) with
the marking maintained incrementally. The win was largest exactly where P1
predicted: the backtrack-free, counting-dominated rungs (21x21_80w halved).

E-H9 also carries the campaign's biggest process lesson: it landed *after*
the campaign had formally declared constant-factor work "mined out". That
verdict was scoped to the representation seam (true) but was read as
"per-node cost is minimal" (false). Distinguish "this seam is exhausted"
from "there is nothing left" — P1's data still pointed at counting cost, and
it was still right.

## Cumulative impact

Per-rung search inferences, campaign start → end (all byte-identical output):

| Rung | Start | End | Δ |
|---|---:|---:|---:|
| 9x9, 8 words | 98,058 | 28,212 | −71.2% |
| 9x9, 16 words | 3,059,551 | 674,961 | −77.9% |
| 15x15, 12 words | 448,384 | 101,588 | −77.3% |
| 15x15, 28 words | 1,655,215 | 374,521 | −77.4% |
| 15x15, 32 words | 4,804,420 | 969,070 | −79.8% |
| 15x15, 34 words | 68,508,076 | 13,650,369 | −80.1% |
| 15x15, 36 words | 182,552,102 | 38,275,505 | −79.0% |
| 21x21, 25 words | 1,264,379 | 302,257 | −76.1% |
| 21x21, 80 words | 27,815,454 | 4,487,855 | −83.9% |

Practical meaning:

- **Latency:** the hard tail moved from "seconds" to "sub-second or barely
  above" native (36 words: ~9.5s → ~2.2s); under WASM's ~3x factor, dense
  arranges land in interactive territory.
- **Feasible density under the shipped 500M budget** (measured by P2; robust
  = places across seeds): 9x9 ~16 → **17** words, 15x15 36 → **42**, 21x21
  ~80 → **82+** (the old "≥82 hard-for-all" line broke). These gains came
  from headroom — instances that already fit the budget got cheap — not from
  flipping search-bound failures (P3's decoupling result).
- **Three envelope-guard rungs** (9x9_17w, 15x15_40w, 21x21_82w) were added
  to the benchmark ladder so the widened envelope is ratcheted against
  regression.

## Cross-cutting lessons

1. **Instrument before you build.** P1 cost one probe and killed a large
   refactor (conflict-directed backjumping) that the literature recommended
   but the data showed had nothing to jump over.
2. **Determinism is a verification superpower.** Byte-identical goldens plus
   a deterministic inference-count ratchet meant every optimization could be
   proven behavior-preserving, cheaply, at every step — and the strongest
   changes (E-H9, E-H10b) were additionally verified by full-tree solution
   -count identity across strategies and corners.
3. **Serial composition keeps attribution honest.** Each accepted change
   re-recorded the baseline before the next was measured; two experiments
   that accidentally measured against stale trees produced misleading
   numbers until re-based.
4. **The strict ratchet earns its keep.** Its zero-regression rule forced a
   sloppy variant (E-H10) to be rebuilt cleanly (E-H10b) and ultimately
   surfaced a genuine product-policy trade instead of letting a hidden tax
   merge. (E-H10b — a further −17.6% on the densest rung — remains parked on
   its branch behind a recorded product-evidence trigger; see the ledger.)
5. **Failed experiments are assets when the *reason* is recorded.** The
   rejections (wdeg tie-breaking, letter bitmasks, watched witnesses) each
   closed an avenue with a measured mechanism, which is what makes the
   "what's left" statement at the end of a campaign trustworthy. Full
   entries in `docs/experiments.md`.
