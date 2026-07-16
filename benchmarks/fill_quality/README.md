# Fill-quality benchmark — crosswordsmith `fill` vs ingrid_core

A head-to-head measuring **fill quality**, not just legality: crosswordsmith's
`fill` (plain-dict scoreless mode AND the §8.4a native scored mode) against
[`ingrid_core`](https://github.com/rf-/ingrid_core) (scored backtracking CSP),
all filling the **same grids** from the **same word list** — Spread the
Wordlist (STW), CC BY-NC-SA 4.0. Every placed entry is scored *post hoc*
against STW, independently of anything the engine reports; the engine's own
`--report-json` is cross-checked against that post-hoc scoring. This harness
is the standing regression for `fill --min-score` (design-spec §8.4a).

This is the measured answer to the "load-bearing open question" in
[`docs/research/setter-tool-landscape-2026.md`](../../docs/research/setter-tool-landscape-2026.md)
§D: *does scoreless fill produce materially worse fill, or just occasionally
worse corners?* **Answer: materially worse — down to nonsense strings — plus a
separate completion-rate gap on hard grids** (the completion gap has since
been closed by the §8.4c engine, DP-8 — see "Adopted: the engine result").

> **Not part of `make bench`.** It needs two external, non-bundled deps:
> `ingrid_core` (`cargo install ingrid_core`, needs a Rust toolchain) and the STW
> `word;score` list (download from <https://www.spreadthewordlist.com/>). Run:
> `STW=/path/to/spread-the-wordlist.txt ./run.sh`. Any `word;score` list works
> as `$STW` — the CWL measurement below reused the same driver.

## Files

- `gen_grids.py` — emits the benchmark grids in both crosswordsmith JSON and
  ingrid_core ASCII (`#`/`.`) form, so both tools fill identical masks. The
  easy quality grids plus `amer11`, an authored standard American-style 11×11
  (fully checked, 180° symmetric, spanning 11-letter slot), asserted valid at
  generation time by `check_american`.
- `score_fill.py` — scores a fill's entries against STW (n, mean, min, #below-50,
  #junk). Reconstructs ingrid's entries from its output grid + the mask. Also
  cross-checks the native fill's `--report-json` sidecar against its own
  post-hoc stats (n/min/belowThreshold exact; mean within 0.1 — SWI rounds
  half away from zero, Python half to even, so a .x5 tie differs by exactly
  0.1 in the last decimal) and exits non-zero on disagreement. Scores are
  keyed by the engine's own normalization (A–Z fold, digits/punct squeezed,
  duplicate keys keep the max score) so post-hoc lookups match what the
  engine actually played.
- `run.sh` — driver: builds the plain + `score>=50` dicts, fills each grid four
  ways (crosswordsmith full / `>=50`-prefiltered dict / **native
  `--min-score 50`** / ingrid_core), prints the comparison. This one GATES
  (non-zero exit on a report-json disagreement).
- `matrix.sh` — the FS-3(b) **completion × min-score frontier** (the
  search-power axis): every mask — the easy four, the American 11×11, and the
  three bundled blocked stock grids — × `--min-score` 1/30/50, for both
  engines (ingrid capped at 90s). Non-gating by design: an incomplete fill
  here is a data point. Evidence for the FS-4 search-power decision pass.

## Results (STW snapshots 2026-07-01 / 2026-07-15; score 50 = the "clean" benchmark)

### Fill quality, on grids both tools complete

The **native `--min-score 50`** column (§8.4a scored fill, added 2026-07-15)
is the FS-1 acceptance gate: it must — and does — match the `>=50`-dict
column (mean/min 50, zero below-clean) on every grid, and its
`--report-json` sidecar agreed with `score_fill.py`'s post-hoc stats on all
four fills. Full-dict numbers were identical across both snapshots.

| grid | crosswordsmith (full dict) | crosswordsmith (`score>=50` dict) | crosswordsmith (`--min-score 50`, native) | ingrid_core (`--min-score 50`) |
|---|---|---|---|---|
| open4 (4×4) | mean **41.2**, min 20, **5**/8 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below ✅ | mean 50, min 50, 0 below |
| open5 (5×5) | mean **27.0**, min 20, **10**/10 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below ✅ | mean 50, min 50, 0 below |
| mini7 (7×7) | mean **42.2**, min **0**, **9**/23 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below ✅ | mean 50, min 50, 0 below |
| mini9 (9×9) | mean **41.8**, min 20, **14**/28 below-clean | mean 50, min 50, 0 below | mean 50, min 50, 0 below ✅ | mean 50, min 50, 0 below |

Scoreless MRV grabs the alphabetically-first fitting strings. The `open5` fill is
the vivid case:

```
crosswordsmith (full)   : AAAAA AAAAQ ALEDA ABAAB AMEAL ASALE ALBMS AEAEA ADAAL QABLE
crosswordsmith (>=50)   : AAHED AIART ANNIE HEDDA EAGER DROSS INEAR ANDGO RIDES TEARS
ingrid_core (min 50)    : ASSTS SHERA SAMUS ERIES SISSY
```

**Two findings:**

1. **Scoreless fill is materially worse** — `AAAAA`, `AAAAQ`, `ABAAB`, `AEAEA`
   are non-words. Mean scores sit at 27–42 with entries at score 0.
2. **A crude `score>=50` *dict prefilter* fully recovers ingrid-parity quality**
   (mean/min 50, zero below-clean) on every completable grid — the cheapest form
   of scored fill (filter `--dict` by score, or an in-engine `--min-score`)
   already closes the quality gap here.

### Completion / search power, on a hard grid (`blocked_13a`, full 13-length slots)

> **Historical (pre-§8.4c MRV engine).** This gap is CLOSED as of the DP-8
> build — see "Adopted: the engine result" below. The table is kept as the
> baseline the DP-6/7/8 arc measured against.

| tool | result |
|---|---|
| crosswordsmith MRV engine (any dict) | **did not complete** — budget-exhausted at ~20s |
| crosswordsmith MRV `--min-score 50` (native, 2026-07-15) | **did not complete** — budget-exhausted at ~21s (the prune removes 195,725 of 315,903 words and shrinks every domain; scoring never lifts the search ceiling, per §8.4a) |
| ingrid_core `--min-score 20` | completes in **10.7s** (mean 45.0, 17 below-50) |
| ingrid_core `--min-score 30` | completes in **9.7s** (mean 44.4, 22 below-50) |
| ingrid_core `--min-score 50` | **times out** (>90s) — even the scored CSP can't hit "clean" on this grid |

This was a **second, distinct** competitiveness axis: the fixed-budget MRV
tree could not complete a standard 13×13 with full-length slots that ingrid
solves in ~10s — orthogonal to scoring, and the harder engineering gap. It
drove the DP-6 → DP-7 → DP-8 arc below, which replaced the MRV tree with
the §8.4c MAC + dom/wdeg core. (Note even ingrid can't fill `blocked_13a`
at the clean-50 bar, so that grid is genuinely hard, not just hard for
crosswordsmith.)

### The completion × min-score frontier (FS-3(b), 2026-07-15)

`matrix.sh` sweeps every mask × `--min-score` ∈ {1, 30, 50} for both engines
(STW snapshot 2026-07-15; ingrid capped at 90s; crosswordsmith rows are
deterministic, ingrid rows are an informational reference):

| grid | threshold | crosswordsmith (native) | ingrid_core |
|---|---|---|---|
| open4 / open5 / mini7 / mini9 | 1 | ok — means 43.0–47.8 | ok |
| open4 / open5 / mini7 / mini9 | 30 | ok — means 46.0–49.1 | ok |
| open4 / open5 / mini7 / mini9 | 50 | ok — mean/min 50, 0 below | ok |
| **amer11** (fully-checked American 11×11, spanning 11-slot) | 1 | **ok** — mean 44.3, min 20 | ok |
| **amer11** | 30 | **ok** — mean 46.4, min 30 | ok |
| **amer11** | 50 | **ok** — mean/min 50, 0 below | ok |
| blocked_13a (UK 13×13, full-length slots) | 1 | NOT completed | ok |
| blocked_13a | 30 | NOT completed | ok |
| blocked_13a | 50 | NOT completed | **NOT completed** (>90s) |
| blocked_13b (UK 13×13) | 1 / 30 / 50 | NOT completed | **NOT completed** (>90s) |
| blocked_15a (UK 15×15) | 1 / 30 / 50 | NOT completed | **NOT completed** (>90s) |

**What the frontier shows (FS-4 sizing evidence):**

1. **The ceiling is NOT "standard masks fail."** `amer11` — fully checked,
   with a spanning 11-letter slot — completes at every threshold, clean at
   50. The budget ceiling bites specifically on the blocked UK-style grids.
2. **Two of the three bundled hard grids are hard for ingrid too.**
   `blocked_13b` and `blocked_15a` complete for *neither* engine at any
   threshold. The only measured crosswordsmith-vs-ingrid search-power gap on
   this whole set is **`blocked_13a` at `--min-score` ≤ 30** — that row is
   FS-4's concrete reference target (ingrid: ~10s).
3. **Pruning never flips completion, in either direction, anywhere on this
   set** — consistent with §8.4a's "never present `--min-score` as helping
   completion" (and it didn't hurt completion here either).
4. **Score-descending ordering alone buys quality.** At `--min-score 1` (no
   meaningful prune) the native means beat the scoreless full-dict column of
   the quality table by 5–16 points (open4 41.2 → 46.3, open5 27.0 → 43.0,
   mini7 42.2 → 47.8, mini9 41.8 → 47.5) — the FS-1 *ordering* benefit,
   before any threshold is applied.

### Budget & ordering probes on the reference row (FS-4 / DP-6 evidence, 2026-07-15)

Two white-box probes against the frontier's sole measured gap row
(`blocked_13a`, STW, `--min-score` ≤ 30 — ingrid_core ~10s), both run in a
throwaway git worktree against the engine at `23e3772`:

- **Budget ladder** — `fill_budget` scaled ×2 / ×4 / ×10 / ×20 over the
  hardcoded 800M inferences, at `--min-score` 30 and 1: **not completed at
  any rung**. ×20 = 16×10⁹ inferences ≈ 354s of search, wall time scaling
  linearly with budget and zero completion signal. Budget buys latency,
  never completion (the perf campaign's finding, now confirmed from the
  product side): the tail is exponential.
- **Ordering-diversity probe** (seeded-restart proxy) — the plain lex
  tiebreak replaced by 8 different salted deterministic reorders *within
  equal-score bands* (score-descending preserved — exactly the perturbation
  a `fill --seed` would apply), each at the default budget: **0/8
  completed** (all ~22s budget-exhausted). No lucky-ordering escape exists
  on this instance.

**Conclusion (DP-6, design-spec §10):** only crossing-aware
forward-checking — ingrid's pruning class — can close the `blocked_13a`
row. `--budget` and `fill --seed` are honest control/variety levers
(spec'd in §8.4b) but measurably not completion fixes.

### MAC probes on the reference row (DP-7 evidence, 2026-07-15)

DP-6 named "crossing-aware forward-checking" as the path to the row; DP-7
tested that claim with [`probe_mac.pl`](probe_mac.pl) — a throwaway
maintaining-arc-consistency filler (AC-3 to fixpoint over the crossing
graph, bignum bitset domains in the F-H2 `build_masks/2` bit-space) run
against the engine at `fd8768b`, plus a source read of ingrid_core 1.3.1
(the reference implementation: MAC + **dom/wdeg conflict-weight learning**
+ a **fillability-dominated value order** (weight ×900 vs score ×5) +
weighted-random top-3 picks + randomized restarts with weight carryover —
a full Balafoutis adaptive-CSP stack, *not* plain forward checking).

Probe validity: fills the 3×3 and 9×9 fixtures (engine parity; the engine's
own 9×9 solution replays cleanly through the probe's propagation), and MAC
collapses node counts exactly as advertised on completable grids — the 9×9
that costs the engine 14.4M inferences takes the probe 4,114 nodes (1,114
with the fillability order).

On the reference row (`blocked_13a`, STW, `--min-score 30`; ingrid_core
**8.8s** on this machine, re-timed with the same wordlist):

| variant | wall cap | result |
|---|---|---|
| deterministic MAC, §8.4a score-desc value order | 900s (naive), 1800s (flat masks) | **timeout**, 191,297 nodes / 22.3M domain-updates |
| MAC + restarts (top-3 [4,2,1] random picks, cap 500 ×1.5, fresh seed per attempt) | 1800s | **timeout**, 14 attempts, 211,115 nodes |
| MAC + restarts + **fillability value order** (Σ log10(letter-support popcount+1) desc — ingrid's dominant signal) | 1800s | **timeout**, 13 attempts, 186,954 nodes |

**Conclusion (DP-7, design-spec §10) — since AMENDED by the dom/wdeg
follow-up below (2026-07-16), which overturns the "new search core" sizing:**
MAC-class propagation is necessary
but measurably NOT sufficient in this engine: at the probe's ~10ms/node
(full-AC sweeps over ~30k-bit masks; already 3× faster than the naive
form), 30 minutes buys ~2×10⁵ nodes and zero completions under three
different search strategies. Closing the row needs the *rest* of ingrid's
stack — dom/wdeg learning and, decisively, incremental-AC node rates 2–3
orders of magnitude beyond the probe — i.e. a new search core, not a
pruning patch. Two honesty pins: (1) the fillability value order that
powers ingrid's completion *inverts* the §8.4a quality-first candidate
order (score-desc is measured to lift fill means 5–16 pts — the two
contracts pull in opposite directions); (2) on completable grids MAC is
already *slower* end-to-end than the engine's counting search unless the
per-node cost is engineered down (the 9×9: engine ~1.4s vs naive-probe
6.1s — the node-count win does not survive the per-node cost).

Caveat for re-runners: SWI's `\` arithmetic misevaluates at the int64
boundary (`\(1<<63)` yields `2^63-1` — a 64-bit complement — instead of the
documented unbounded one's complement `-(2^63)-1`, silently wiping bits
≥ 63 of a mask AND-chain). Clear bits with `xor`, never `/\ \(1<<B)`;
found the hard way. Candidate upstream report.

### dom/wdeg follow-up: the missing ingredient (DP-7 amendment, 2026-07-16)

The one component of ingrid's stack the table above never isolated was
**dom/wdeg conflict-weight learning** — probed the next day as
`probe_mac_dwd/6` (same MAC machinery; every crossing carries a weight,
init 1.0, +1 whenever its revision wipes a domain, held in a global
`nb_setarg` compound so learning survives backtracking, persisting across
restarts with ×0.99 aging; variable order picks the slot minimizing
`popcount(dom) / Σ weights of crossings to unfilled neighbors` instead of
bare MRV). Each completed fill is **hard-verified** by binding the placed
words onto the grid's shared cell variables (crossing inconsistency fails
unification) and is byte-reproducible per seed.

Same reference row (`blocked_13a`, STW, ingrid_core 8.8s):

| variant | threshold | result | quality |
|---|---|---|---|
| MAC + restarts + **dom/wdeg**, **§8.4a score-desc value order** | 30 | **completes: 7.1s / 726 nodes** (seed 12345; seeds 777 / 424242: 11.4s / 908 and 9.2s / 451 — **3/3 seeds**) | **mean 44.81**, min 30 (ingrid @30: mean 44.4) |
| MAC + restarts + dom/wdeg, score-desc order | 1 | **completes: 12.6s / 876 nodes** | — |
| MAC + restarts + dom/wdeg, fillability order | 30 | completes: 24–25s / 2,010 nodes | mean 43.89, min 30 |
| MAC + restarts + dom/wdeg, score-desc, **greedy pick — RNG-free** | 30 | **completes: 125.7s / 10,517 nodes / 7 attempts** | mean 44.63, min 30 |
| MAC + restarts + dom/wdeg, score-desc, greedy pick — RNG-free | 1 | **completes: 9.1s / 604 nodes / 2 attempts** | — |

(For scale on the learning effect alone: the 9×9 drops from 14,512 nodes
(restarts, no learning) to 2,175 with dom/wdeg on the same score-desc
order, and from 1,114 to 369 on the fillability order.)

**What this overturns in the DP-7 conclusion:** two of its three "you need
the whole ingrid stack" claims are falsified — (1) **native node rates are
NOT required**: ~726 nodes at the probe's ~10ms/node is seconds, not
hours; (2) **the fillability quality-inversion is NOT required**: the
§8.4a score-first value order not only completes, it *beats* fillability
here on both wall time and fill quality, at quality parity with
ingrid_core on the same row. What stands: (3) plain MAC propagation alone
is insufficient (the 2026-07-15 table), and MAC still regresses easy
completable grids end-to-end at spike per-node costs. The corrected
sizing: the row is closable by a **bounded upgrade** — MAC propagation +
dom/wdeg + seeded restarts, keeping the quality-first candidate order —
not by adopting ingrid's engine wholesale. Adoption remains its own
decision (the §8.4/AC-FILL-3 contract churn priced in DP-6/DP-7 is
unchanged); see the amended DP-7 entry in the design spec.

**Determinism follow-up (2026-07-16, same day):** the `greedy` rows above
answer "does the recipe need randomness at all?" — **no.** With strict
best-first candidate picks the PRNG is never consulted (verified: `Pick =
greedy` short-circuits before any draw), and restart diversity comes only
from the learned weights reordering slot selection between attempts. The
RNG-free variant still closes both gap thresholds; the seeded top-3
randomness is an *accelerator*, not a requirement (~18× on the hardest
row: 7.1s randomized vs 125.7s greedy @ 30). So an adopted engine could
keep the strongest form of today's contract — deterministic with no RNG on
the default path — and treat `--seed`/`--shuffle` as the opt-in fast/variety
lever, exactly §8.4b's existing shape.

### Adopted: the engine result (DP-8 → §8.4c, built 2026-07-16)

The recipe above IS the shipped `fill` engine as of the DP-8 build
(design-spec §8.4c; plan `docs/plans/fill-mac-dwd-implementation.md`).
Product-CLI numbers on this machine, **default budget (800M), default
path** (no seed):

| run | result | quality |
|---|---|---|
| `blocked_13a` × STW `--min-score 30` (the reference row) | **filled: 2m20s / 7 attempts** | **mean 45.0**, min 30 (ingrid 8.8s / mean 44.4; probe 7.1s / 44.81) |
| `blocked_13a` × STW `--min-score 1` | **filled: 19s / 4 attempts** | mean 38.7 |
| 15×15 / 17×17 / 21×21 × full ENABLE (identity-ladder rungs) | filled: 3.2s / 3.5s / 3.9s | — |
| 17×17 × ENABLE-50k | filled: 1.3s / 4 attempts | — |

Two C3 build findings the probes above could not see (they only ran
blocked grids; the pure-greedy engine failed the OPEN-grid ladder rungs
until both landed): per-node candidate enumeration must be **lazy** (the
eager materialize-the-domain form burns the budget at O(popcount ×
mask-width) bignum work per node), and the default path needs **pinned
diversification** — top3 restarts on an engine-constant PRNG stream from
attempt 2, plus a pinned load shuffle for dicts with no ordering
information (plain / single-band; lexicographic clumps defeat top3
diversification, while scored multi-band lists must KEEP their lex tail —
in-band permutation starves this reference row). The engine default is
therefore "deterministic as a pure function of the input" rather than
"draw-free": full detail in §8.4c and the plan's C3-findings section.
Wall-clock note: the engine's 2m20s vs the probe's 7.1s on the reference
row is the priced cost of the deterministic pinned stream vs the probe's
lucky random seeds (per-seed variance is real — engine `--seed 7` on the
same row budget-exhausts at 8m47s with a clean `not proven` exit), plus
CLI load/mask-build overhead; the AC-FILL-12 gate is completion under the
default budget at ≥ old-engine quality, and it passes with quality above
ingrid's.

## CWL measurement (2026-07-15) — evidence for the CWL-bundling decision pass

The FS-6(b) research confirmed the [Collaborative Word
List](https://github.com/Crossword-Nexus/collaborative-word-list) (CWL) is
**MIT + fully scored** (567,657 `word;score` entries, integer 0–100), making
it the first bundleable scored candidate (design-spec §8.5). The harness was
run with CWL as `$STW` to convert its two recorded caveats — staleness
(frozen 2023-02) and undocumented band semantics — into data:

| grid | crosswordsmith (`--min-score 50`, native) | ingrid_core (`--min-score 50`) |
|---|---|---|
| open4 | mean **81.2**, min 50, 0 below, 0 junk, report-json agrees | mean 78.8, min 50, 0 below |
| open5 | mean **78.5**, min 50, 0 below, 0 junk, report-json agrees | mean 77.0, min 50, 0 below |
| mini7 | mean **82.4**, min 50, 0 below, 0 junk, report-json agrees | mean 77.6, min 50, 0 below |
| mini9 | mean **78.0**, min 50, 0 below, 0 junk, report-json agrees | mean 77.1, min 50, 0 below |

**Findings, in decreasing order of importance for DP-6:**

1. **Capacity ceiling (engine finding, list-independent):** the full 567k-word
   CWL dict **crashes `build_index`** at SWI's default 1GB stack limit — a raw
   stack-overflow ERROR at dict load, not a clean §5.1-style failure report.
   STW's 315k loads fine, so the ceiling sits somewhere in between. Any list
   this size needs `--min-score` (CWL at 50 → 252,564 words, loads cleanly);
   the unpruned "full dict" column is therefore **unmeasurable for CWL**.
   Worth an engine robustness item: fail cleanly (or stream the index build)
   instead of overflowing.
2. **Quality: drop-in clean.** Native `--min-score 50` over raw CWL fills all
   four grids at min 50 / 0 below-clean, with `--report-json` agreeing with
   post-hoc scoring, and means (78–82) at or above ingrid_core on the same
   list. On these small grids the 2023 freeze does not hurt: CWL's clean band
   behaves like STW's.
3. **Digit-entry hygiene:** CWL keeps digit-bearing entries (`MP3J;100`,
   `0AD;25`, 992 such lines). The engine's documented A–Z fold squeezes
   digits, so `MP3J;100` becomes a 3-letter `MPJ` at score **100** — which
   score-descending ordering then plays *eagerly* (the open4 fill opens with
   `MPJS`). A CWL-derived dict should pre-filter digit-bearing entries (or
   the bundling pass should ship a cleaned derivative). This is also why
   `score_fill.py` keys scores by the engine's fold — raw-keyed lookups
   mislabel such entries as junk.
4. **Search power unchanged (expected):** `blocked_13a` with CWL at
   `--min-score 50` is still not proven within budget (~26s) — the ceiling is
   the search, not the list.

### CWL on the §8.4c engine (2026-07-16) — DP-9 grounding

Re-run for the DP-9 bundling decision after DP-8 replaced the search core.
Same harness (`run.sh` with CWL as `$STW`, plus the matrix rows for `amer11`
and `blocked_13a` only — `blocked_13b`/`15a` stay untouched per the DP-6
report-don't-chase pin), against a **pinned snapshot**: upstream commit
`2efe76e11ef3` (`xwordlist.dict` last changed 2023-02-12; repo re-verified
MIT at source, not archived), raw sha256
`a945a839a5f1e6f48caf9c8de446e5cd85f3567d7f62afcf54c6b738e8906ff4`. Also
measured: the cleaned derivative DP-9 proposes to ship — digit-bearing
entries dropped (exactly the 992 lines finding 3 predicted; 567,657 →
566,665), `cwl50` = its `score ≥ 50` subset (252,200 entries, 3.7MB).

| grid | crosswordsmith (`--min-score 50`, native, §8.4c) | ingrid_core (`--min-score 50`) |
|---|---|---|
| open4 | mean **83.1**, min 50, 0 below, report-json agrees | mean 78.8, min 50, 0 below |
| open5 | mean **78.5**, min 50, 0 below, report-json agrees | mean 77.0, min 50, 0 below |
| mini7 | mean **81.7**, min 50, 0 below, report-json agrees | mean 77.6, min 50, 0 below |
| mini9 | mean **80.9**, min 50, 0 below, report-json agrees | mean 77.1, min 50, 0 below |
| amer11 | mean **78.1**, min 50, 0 below (44 slots, **6.4s end-to-end** with the `cwl50` artifact) | ok |

**Findings (supersede the 2026-07-15 quality table for the current engine):**

1. **Quality: §8.4c + CWL beats ingrid_core on all five completable masks**
   (the §8.4b engine had tied/edged it on four). `run.sh`'s report-json gate
   passes on all four gated grids.
2. **Capacity: the engine has TWO stack envelopes, and both fail as raw
   crashes, not clean §5.1 reports.** (a) *Load*: the full 566k cleaned list
   still blows the default 1GB stack at `build_index` (finding 1 stands on
   §8.4c). (b) *Search* — NEW: at `--min-score 30` the dict loads (437,400
   words survive) but `blocked_13a` dies ~17s in with a global-stack
   overflow inside `mac_support`/`mac_revise` — the §8.4c bignum candidate
   masks over full-13 domains at this band width. `amer11` @30 is fine
   (mean 67.7, min 30). The ≥50 band (252k) is inside both envelopes.
   *(Since resolved — DP-10 (§10), same day: both modes now report as an
   ordinary one-line §5.1 failure with remedies (AC-FILL-15) instead of the
   raw dump. The envelope NUMBERS here are unchanged and still current.)*
3. **The DP-8 reference-row completion is list-specific.** `blocked_13a` @30
   with CWL is envelope crash (b) above — not comparable to STW's 2m20s
   completion. At the clean floor the row stays hard for *both* engines with
   CWL exactly as with STW: native budget-exhausts cleanly (`not proven
   within budget`, 3m57s under the default 800M budget with the `cwl50`
   artifact — the ordinary §5.1 outcome, no crash), ingrid_core
   NOT-COMPLETED at the 90s cap. Do not present the CWL bundle as a
   search-power demo on the UK stock grids: **no loadable CWL floor fills
   `blocked_13a`**.
4. **Digit filter validated:** dropping `[0-9]`-bearing lines removes exactly
   the 992 recorded entries; the eager `MPJS`-artifact from finding 3 is gone
   from the cleaned derivative.
5. **The `amer11` + `cwl50` pairing is the demo:** 6.4s load+fill, mean 78.1,
   0 below-clean, and the filled layout PASSes `lint --profile blocked-uk`
   (222 checks); the `amer11` mask itself PASSes `stockgrid_report` — the
   basis for promoting it into `grids/` (DP-9).

### Seeded `cwl50` load after A2-KS (2026-07-16)

The §8.5 seeded-load follow-up retained the historical permutation exactly
while replacing its O(n²) survivor-list selection with Fenwick-indexed rank
selection. Same-machine end-to-end command on the trivial 3×3 grid, seed 7:

| implementation | wall | peak RSS | outcome |
|---|---:|---:|---|
| historical selection walk | 171.14s | 785,912 KiB | filled |
| A2-KS exact-sequence replay | **9.57s** | 819,404 KiB | filled |

That is **17.88×** faster with a 4.26% peak-RSS increase. The sequence gate
compared the old and new permutation plus the following PRNG draw for five
seeds at every length 0..256 (1,285 cases); all matched. A seeded full-ENABLE
(172,823 words) smoke run completed in 4.91s. Goldens, the 11-rung identity
oracle, both inference ratchets, fuzz, and CLI/WASM parity all remained clean;
nothing was re-baselined.

### B0 MAC instrumentation (2026-07-16)

The Balafoutis follow-up campaign first measured the shipped §8.4c policy
through an exact-replay benchmark twin. Full tables, commands, correction
record, and the runnable rig live in
[`../results/2026-07-16-fill-b0-mac-instrumentation.md`](../results/2026-07-16-fill-b0-mac-instrumentation.md).

| signal | STW `blocked_13a @30` | STW `blocked_13a @1` |
|---|---:|---:|
| attempts / nodes | 7 / 11,318 | 2 / 729 |
| queue pops at length >=4 | 99.51% | 98.08% |
| fruitful revisions | 27.89% | 28.93% |
| learned bumps | 10,945 | 599 |
| top-quartile bump share / Gini | 83.54% / 0.737 | 91.49% / 0.799 |
| propagation/support share of search wall | 99.54% | 99.49% |

None of B0's kill-tests fired: queue ordering, alternate credit, aging, and
probing-init all have measured structural surface. This is permission to
probe them, not evidence that any wins. Product/twin fills matched on both
authority rows and all five quality masks; all 11 ladder rows filled; the
corrected `cwl50 @50` stretch remained `not_proven` at 800M inferences.
Product goldens, identity, and inference ratchets stayed unchanged.

## Caveats / how to extend

- **`amer11` is also shipped** as `grids/amer11.json` (promoted at DP-9 as
  the bundled-CWL demo pairing) — `gen_grids.py` stays the benchmark's own
  source; if either mask is ever edited the two must be re-synced by hand.
- **Mask sample.** Eight masks now (easy four + American 11×11 + three blocked
  stock grids); a fuller benchmark still wants an American 15×15 (the
  next authored-mask candidate — `check_american` in `gen_grids.py` makes
  authoring safe) and several STW snapshots.
- **`junk`=0 by construction** — crosswordsmith's dict here is *derived from* STW,
  so every placeable word has an STW score; the quality signal is mean/min/below-50.
- **Post-hoc scoring is the independent check** — `score_fill.py` scores fills
  without trusting the engine, which is exactly how the native
  `cs_minscore_<g>.json` variant and its `--report-json` are gated (the two
  must agree; `score_fill.py` exits non-zero when they don't).
- **CI subset:** this head-to-head needs the two external deps and stays
  on-demand; the no-deps regression for the same §8.4a behavior (AC-FILL-5/-7)
  lives in `make test` — plunit + goldens over the bundled original scored
  fixture (`fixtures/dict_scored_sample.txt`, `tests/golden/fill_scored*`).
- **FS-3(b) is done** (2026-07-15): `matrix.sh` is the completion × min-score
  frontier over eight masks. It is a search-power *measurement* feeding the
  FS-4 decision pass, not part of the scored-fill (FS-1) acceptance gate, and
  deliberately non-gating.
