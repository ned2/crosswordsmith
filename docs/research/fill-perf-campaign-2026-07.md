# The fill performance campaign (2026-07): what worked and why

A distillation of the four interventions that shipped during the July 2026
fill performance campaign — the intuition behind each, and what it bought.
This is the "greatest hits" companion to the full experiment ledger in
`docs/experiments.md` (which also records everything that *didn't* work, with
evidence); the sibling docs (`arrange-perf-campaign-2026-07.md` — the playbook
this campaign reused, `swi-vm-wasm-performance.md` — the WASM cost model)
hold the groundwork. The red-teamed plan is
`docs/plans/fill-perf-campaign.md`.

**Headline:** end-to-end CLI latency on the hardest benchmark rung (a 9x9
dense grid against the full 172k-word ENABLE dictionary) went from **5.18s to
0.87s** (0.77s with opt-in bitmasks); load-dominated rungs went from 2.5–4s
to **0.23–0.50s**. Warm dictionary load at 172k words: **3.47s → 0.099s**.
Under the machine-independent metrics: dictionary-load inferences −59.7%,
search inferences −21.7% to −57.4% per rung — with byte-identical output
throughout, proven at every step.

## Context: the shape of the problem

`fill` takes a legal blocked grid (a stock-grid mask) and a dictionary and
fills every slot by DFS with MRV ordering. The grid is a set of slots sharing
*unbound cell variables* at crossings (unification propagates letters; the
WAM trail undoes them on backtrack). The dictionary is loaded into
length-keyed buckets plus an assoc index mapping `(len, pos, char)` to the
ordset of matching bucket positions; MRV counts a slot's candidates by
intersecting the index sets for its bound cells. Determinism is a product
requirement: same input, same fill, byte for byte.

Unlike arrange (one cost: search), fill has **two costs with different
owners**: a per-invocation *load* tax (parse 172k words, build the index)
and the *search* itself. The campaign's defining discovery was about their
ratio.

Three measurements framed everything:

- **P-F1 (attribution probe):** dictionary load is **58–84% of end-to-end
  wall on 10 of 11 rungs**. The user's wait was mostly not the search at
  all — so the "startup posture" phase outranked the search-optimization
  phase, inverting the plan's emphasis. Within load, a single lambda in the
  parse loop was 73–76% of load inferences; within search, counting was
  59–90%, with MRV's per-node full recount alone at 79.5–85.6%. The probe
  also *refuted* the plan's top-ranked gap (candidate materialization
  measured ≤0.6%), closing two hypotheses (F-H4, F-H5) before anything was
  built. And it flagged a measurement blind spot that shaped the whole
  campaign: C-native work (`keysort`, `length/2`) costs wall but ~zero
  inferences — the ratchet cannot see it.
- **The WASM bignum gate probe (before F-H2, mandated by the red-team):**
  the feared "bitsets are a native win but a WASM loss" was refuted —
  the AND+popcount kernel stays ≥9x ahead of ordset intersection under
  WASM/LibBF at realistic set sizes. But the probe found the *real*
  constraint elsewhere: building the masks costs ~18% of load inferences
  and 27–36% of load wall. A counting win must never be paid for out of
  the dominant load term. This condition ended up binding (see F-H2).
- **A re-attribution before building F-H2:** after two accepted
  experiments, the original shares were stale; a fresh measurement (38–47%
  of fill-phase wall still in counting) had to pass a pre-registered 20%
  gate before construction was allowed. Shares go stale as wins land —
  re-measure before each build.

## The four shipped interventions

### 1. F-L1 — the parse-loop lambda → first-order code (−59.7% load inferences)

`normalize_word/2` filtered each line's characters with
`include([C]>>char_type(C, alpha), Cs, Letters)`. The module never imports
`library(yall)`, so that lambda was not compiled away — it was *interpreted*:
a meta-call plus a lambda copy **per character**, ~12.3 inferences each,
across the 1.57M characters of ENABLE = 19.4M of the 26.6M load inferences.
Replaced with a first-order `alpha_chars/2` recursion making the identical
`char_type(C, alpha)` decision per character. Load inferences: 26.6M → 10.7M
(−59.7%); warm load wall 3.47s → 1.96s.

**Intuition: an interpreted lambda inside a per-element loop is an
interpreter running inside your parser.** In SWI-Prolog a yall lambda in a
meta-called position is re-interpreted per call unless the module imports
yall and the compiler can expand it — check which regime you are in before
putting a lambda in a hot loop. The mechanism was nailed completely: a
control experiment (importing yall made the lambda cost collapse to
named-helper price) plus a per-character cost model that predicted the
product's measured totals *to one inference*. A pleasing detail: the minimal
fix (named helper passed to `include/3`, the variant the probe had measured
at −53%) was implemented, measured, and *rejected* — the first-order loop
also deletes `include/3`'s per-element `call/2` and won by another 14%.

### 2. F-H1 — incremental MRV counts (−21.7..−57.4% search inferences)

Every search node recounted candidates for *every* unfilled slot to pick the
most-constrained one — 79.5–85.6% of all search inferences. But a slot's
count is a **pure function of its cells' binding state**, and placing a word
binds only the placed slot's previously-free cells. So: pair each slot with
its count once at the root, thread the paired list down the search, and
after each placement recount *exactly* the slots that cross a newly-bound
cell — measured at most 5 slots per placement, against slot lists up to 154.
Backtracking restores the old counts for free, because the threaded list is
pure data: the caller's copy was never touched.

**Intuition: don't recompute a pure function whose inputs didn't change —
and in Prolog, "undo on backtrack" is free if state is threaded rather than
destructive.** (The same VM-does-the-bookkeeping insight as arrange's E-H2
grid-of-variables, applied to derived state instead of the store itself.)
Equivalence was engineered, not hoped for: the selection key
`c(Count, Start, Dir)` is a total order (Start+Dir is unique), completed
slots stay in the counted set so a 0-count dead slot still fails its branch
exactly as before, and crossing detection compares ground cell *numbers* so
no shared cell variable is ever inspected or copied (the findall/copy trap
that fill's own header warns about). Hardest rung: 34.9M → 14.9M search
inferences. A bonus the ratchet couldn't see: fewer `candidate_count` calls
also means fewer C-native `length/2` bucket walks, so one rung's wall
dropped 40.7% against a 29.9% inference drop. One honest cost: the deepest
rung (backtrack depth 153) pays +82% RSS for the per-placement list spine —
recorded, bounded, with an assoc-based mitigation specified if the WASM
footprint ever needs it.

### 3. F-L2 — the precomputed index artifact (load 1.96s → 0.099s, ~20x)

After F-L1 the load residue was `build_index` — dominated by a C `keysort`
the inference ratchet is blind to. All of it recomputes a **pure function of
a frozen file** on every invocation. So compute it once: `fill --dict D
--save-index F` serializes the *exact runtime terms* (the length buckets and
the assoc index — serializing the built term preserves even the AVL shape,
making "loaded == built" a checkable `==`), and `fill --index F` loads them
back in 0.099s instead of 1.96s. End-to-end, load-dominated rungs went
2.5s → 0.23s (11x); even the search-dominated worst case gained 3.7x.

**Intuition: a cache is just a pure function you refuse to re-run — but
formats must be measured, not assumed, and staleness must refuse, not
guess.** The documented fast path (`.qlf`) *lost* the bake-off to
`fast_read` on every axis (2.7x load, 46x build, smaller file); it stays
documented as the WASM-aligned alternative since fastrw's binary format is
SWI-version-bound. That version-boundedness is why the artifact embeds the
dictionary's SHA-256 and the SWI version and **refuses** on any mismatch
with a clear error — a build-time cache that silently rebuilt or silently
served stale data would be worse than no cache. The raw text path was left
byte-for-byte untouched (the builder *calls* it), so the 11 inference-gated
rungs still pin it, and a second identity oracle proves artifact-mode output
against the same pinned digests.

### 4. F-H2 — bitset counting, opt-in (fill-phase wall −24..−33% where it matters)

With counting still 38–47% of fill-phase wall post-F-H1 (the pre-registered
gate passed), the counting kernel was swapped: per `(len, pos, char)` a
bignum whose bit *i* is bucket position *i* (derived from the ordsets, so
agreement is by construction), and a count is `popcount(M1 /\ M2 /\ ...)`
instead of chained `ord_intersection` + `length`. Candidate *enumeration*
stays on ordsets — the gate probe had measured bit-iteration as ~7x slower
than an ordset walk, so only the count moved. Masks ride in the artifact
(schema v2) and are built at `--save-index` time, per the gate probe's
binding condition. Fill-phase wall on the hard rungs: −24% to −33%.

**Intuition: a set intersection whose only consumer is its cardinality is a
machine-word AND away from free — but measure the deployment target before
building, and check what the *carrier* costs.** The WASM probe justified
building at all (the kernel survives LibBF at ≥9x). The honest finding came
after: masks inflate the artifact +32%, and unlike construction, **size is
not amortized** — every load pays ~+0.09s at 172k, which out-costs the
counting win on load-dominated grids. So masks shipped *opt-in*
(`--save-index --masks`, named for the search-heavy/WASM use case); the
default artifact is byte-identical to the maskless format. This corrected
the gate probe's "pure win everywhere" deployment clause on the record. Note
the metric inversion: this win is almost invisible to the inference ratchet
(bignum ops are ~3 inferences flat) — it is a wall-only win and was gated
and evidenced as such.

## Cumulative impact

Search inferences per rung, campaign start → end (all byte-identical output;
the search-layer change is F-H1 — F-H2's kernel is inference-blind by
design and opt-in):

| Rung | Start | End | Δ |
|---|---:|---:|---:|
| 4x4 sq, ENABLE 172k | 529,275 | 371,343 | −29.8% |
| 11x11, ENABLE 172k | 671,484 | 514,725 | −23.4% |
| 11x11 seeded, ENABLE 172k | 699,377 | 419,582 | −40.0% |
| 5x5 sq, ENABLE 172k | 1,461,408 | 876,182 | −40.1% |
| 17x17, ENABLE 172k | 2,534,409 | 1,822,010 | −28.1% |
| 21x21, ENABLE 172k | 3,307,580 | 1,853,583 | −44.0% |
| 13x13, ENABLE 172k | 3,754,832 | 2,934,984 | −21.8% |
| 4x4 sq, 50k | 7,738,070 | 6,062,379 | −21.7% |
| 15x15, ENABLE 172k | 10,734,114 | 7,967,168 | −25.8% |
| 17x17, 50k | 19,637,890 | 13,773,804 | −29.9% |
| 9x9 dense, ENABLE 172k | 34,880,750 | 14,870,446 | −57.4% |

Dictionary-load inferences (F-L1): 26,602,930 → 10,724,707 (−59.7%) at 172k;
7,757,017 → 3,164,117 (−59.2%) at 50k. Both layers are now **gated** by the
ratchet at 0.5% (load was promoted from reported-only when it started
carrying wins worth defending).

What the user actually feels (end-to-end CLI, warm):

| | campaign start | raw path now | `--index` | `--index` + masks |
|---|---:|---:|---:|---:|
| load-dominated rungs | 2.5–4s | ~2.5–3s | **0.23–0.50s** | (masks cost here) |
| 9x9 dense (hardest) | 5.18s | 3.23s | 0.87s | **~0.77s** |

- **Interactive use is now real**: with a prebuilt index, every 15x15/21x21
  rung fills in under half a second natively; the hardest rung in the
  ladder is under a second. Inference counts are machine-independent and
  were verified bit-identical native↔WASM by the gate probe, so the search
  wins carry to the browser; the load win carries as the artifact
  (pending the fastrw-under-WASM check, with `.qlf` as the fallback).
- **The product grew a small, explicit surface**: `--save-index [--masks]`
  / `--index`, integrity-checked, off by default, raw path untouched.

## Cross-cutting lessons

1. **Attribute end-to-end, not just the algorithm.** The plan ranked search
   gaps; P-F1 showed the user's wait was 58–84% dictionary load. The
   biggest perceived-latency wins (F-L1, F-L2) were in code the search
   campaign almost skipped. Perceived latency is the product metric;
   the clever algorithm is not where the seconds were.
2. **The inference ratchet has a C-shaped blind spot.** `keysort`,
   `length/2`, GMP bignum ops, and `fast_read` all cost wall and ~zero
   inferences. Two of the four wins were partly or wholly invisible to the
   gated metric and were evidenced by wall medians riding alongside it.
   Trust the deterministic gate for *equivalence and regressions*; trust
   attribution and wall for *what to do next*.
3. **Gate expensive builds on cheap measurements.** Two do-not-skip probes
   (P-F1, the WASM kernel probe) plus a pre-registered re-attribution gate
   before F-H2 killed two hypotheses pre-build, refuted a plausible fear
   with numbers, and surfaced the real constraint (the size tax) that
   shaped F-H2's shipped form. Shares go stale after every accepted win —
   re-measure before committing to the next build.
4. **Byte-identity plus a two-layer ratchet made aggressive changes cheap
   to trust.** Every intervention — including a full search-drive rewrite
   (F-H1) and a second dictionary-load path (F-L2) — was proven
   output-identical on all 11 rungs (F-L2/F-H2 in *three* modes against the
   same pinned digests), and equivalence arguments were engineered into the
   designs (total-order selection keys, ground-only comparisons,
   loaded-`==`-built) rather than asserted.
5. **Honest negatives and self-caught errors are what make the record
   usable.** The ledger holds: a variant rejected for being merely good
   (F-L1's −53% named helper), a hypothesis rescinded by its own probe
   (materialization), a deployment clause corrected by the experiment it
   gated (F-H2's size tax), and a runner's self-caught measurement
   contamination (a git-stash A/B cycle timing error paths; a 100ms timer
   quantization inflating a regression 5x). Full entries, including
   everything closed with a measured mechanism, in `docs/experiments.md`.
