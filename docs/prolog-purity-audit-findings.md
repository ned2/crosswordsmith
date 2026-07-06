# SWI-Prolog Purity & Cut Audit — Findings & Remediation Tracker

A **purity / extra-logical-control-flow** audit of the whole current core, driven by one
empirical question: *are our cuts appropriate, and does removing them help performance?*
Scope: `prolog/crosswordsmith/{core,arrange,fill,lint,export,stockgrid,metrics}.pl`,
`load.pl`, `wasm/client/solve_browser.pl` (+ bench/test harnesses where relevant).
Generated 2026-07-06 by a two-phase multi-agent workflow (see [Method](#method)) on
SWI-Prolog **10.1.10**, against the version-matched manual under
[`reference/swi-manual/`](reference/swi-manual/). Baseline commit `69a3a2d`.

This is the **fourth** audit in the series:

- [`prolog-audit-findings.md`](./prolog-audit-findings.md) — pre-revamp
  `crossword.pl`/`quality.pl` idiom sweep (F001–F026, **closed 2026-06-24**).
- [`revamp-audit-findings.md`](./revamp-audit-findings.md) — revamp correctness /
  spec-conformance (R1–R16, **closed 2026-06-30**).
- [`prolog-idiom-audit-findings.md`](./prolog-idiom-audit-findings.md) — full-codebase
  idiom / predicate-use sweep (P1–P17, **closed 2026-07-02**).
- **this doc** (C1–C47) — cut/impurity census + **controlled experiments**, the first audit
  in the series to *measure* its recommendations before recommending them.

Findings already dispositioned in the earlier three are out of scope and are **not
re-opened** here; where a C-finding touches prior-audit territory it says so explicitly.
In particular the P-audit's [Appendix A](./prolog-idiom-audit-findings.md#appendix-a--verified-correct--do-not-re-flag)
("verified correct — do not re-flag") was honoured throughout, and this doc extends it
(§[Verified-deliberate hand-rolls](#verified-deliberate-hand-rolls--do-not-swap)).

This doc **doubles as the remediation tracker** — see the
[checklist](#remediation-tracker) and append to the [remediation log](#remediation-log).

## How to read each finding

- **severity** — `high` (wrong output / crash / hang on realistic input) > `med`
  (systemic gap, measurement-integrity defect, or a real hot/warm-path cost) > `low`
  (latent / cosmetic-but-real) > `nit` (pure style / recorded-only).
- **behaviour** — `preserving` (cannot change observable output: the 214-plunit suite, the
  8 CLI goldens, the fill identity digests, and the JSON/ipuz/exolve contracts stay
  byte-identical) vs `risk` (a fix *could* change observable output and must be
  golden/diff-checked).
- **flags** — `perf-relevant` (moves a ratchet-gated inference count — arrange
  `search_inf`, fill `search_inf`/`load_inf` — so it needs an explicit
  `bench-record`/`bench-fill-record`), `contract-affecting` (byte-locked output),
  `test-affecting` (plunit asserts on the exact term/signature), `WASM` (matters for the
  browser image).
- **status** — `open` / `applied` / `fixed` / `rejected` / `wontfix`.
- **Line numbers** are as at the audit baseline (`69a3a2d`). The applied Phase-2 diffs
  shift some sites by a few lines.

## Method

**Phase 1 — 8 read-only audit agents** (parallel lanes): two impurity-census agents
(core+arrange; fill+lint+export+stockgrid+metrics+solve_browser), four idiomatic-style
agents (core; arrange; fill+metrics; small modules), one stdlib-reinvention sweep, one
determinism/steadfastness/choicepoint-hygiene agent. Every indexing/determinism claim was
grounded in `docs/reference/swi-manual/jitindex.md` (§2.17 JIT clause indexing, §2.17.1
deep indexing, §2.17.3 body-code indexing) **and empirically verified with live `swipl`
10.1.10 probes** — `call_cleanup/2` determinism probes and `vm_list/1` disassembly — not
taken from the manual alone (the probes *corrected* the manual's suggestion in one
load-bearing place, see §[2b](#b-where-indexing-recovers-determinism-cut-deletion-is-free--and-eliminating-cp-leaks-pays-in-memory)).

**Phase 2 — 6 controlled experiments** (X1–X6), each in an **isolated git worktree** off
`69a3a2d`, one committed step at a time. Every step was gated on: `make unit` (**214
plunit tests, 0 failures**) + **byte-identical goldens** (8 CLI goldens; fill additionally
the 11-rung `check_fill_identity.sh` digests) — and measured on the **deterministic
inference-count benches** (`benchmarks/run_arrange.pl` 5 core + 7 heavy rungs;
`benchmarks/run_fill.pl` 7 core rungs), with double-run baseline confirmation (inference
columns byte-identical across runs before any edit). Wall-clock was treated as
indicative-only (shared machine, see [Methodology appendix](#methodology-appendix)).
X6 additionally ran a designed **negative control** (E8). Steps earned an APPLY verdict
only with all gates green; two census predictions were **dropped by measurement**
(E2, E4 — see §[1f](#f-two-census-predictions-falsified-by-measurement)).

---

## Executive summary

- **Zero high-severity findings, across all eight Phase-1 agents.** The mechanical
  baseline is spotless: zero compile warnings, `check/0` clean, all `list_strings` hits
  intentional. This is the cleanest audit surface in the series (see
  §[What the codebase does well](#what-the-codebase-does-well)).
- **The headline question is now answered empirically** (§[1](#1-the-headline--are-our-cuts-appropriate-and-does-removing-them-help-performance)):
  cut removal is a win **exactly where clause indexing can take over the cut's job**
  (first-argument dispatch, or the `[]`/`[_|_]` two-clause special case) — there it is
  semantics-neutral, inference-neutral *by construction*, and pays in choicepoint-stack
  memory (RSS **−6.3…−6.5%** on the 21×21 heavy rungs) and structural safety. Where
  indexing *cannot* dispatch (body-guard selection over non-discriminating heads), **the
  cut IS the pruning**: the negative control regressed `search_inf` on **7/7 rungs** and
  failed the project's own ratchet. "Delete cuts" is not a direction; *"make clause
  selection indexable, then delete the cut"* is.
- **18 experiment steps earned APPLY and are applied to this worktree** (staged,
  uncommitted; the overlapping X4.3/X5.1 pair merged via git 3-way into 17 distinct
  patches — see the [count note](#applied-step-count-note)). Merged-state gates: full
  `./run_tests.sh` **ALL TESTS PASSED**; `make bench-check` **PASS (no change, 0
  regressions)**; `make bench-fill-check` **PASS** with `search_inf`
  **−0.97% / −2.39% / −2.47% / −0.34% / −4.60% / −3.42% / −2.82%** on
  sq04/g11/g11_seed/sq05/g17/g21/g13 (**6 WINs, 0 regressions**; sq05's −0.34% sits
  inside the 0.5% tolerance). `make bench-fill-record` is owed after commit.
- Net purity movement: **14 of the 42 cuts deleted, 9 of the 12 `once/1` wrappers
  deleted**, one exported predicate (`assign_clue_numbers/2`) made deterministic at
  source, the fill search spine and both counting kernels now cut-free, the greedy
  constructor's grid-copy tax removed (**wall −17…−18%** on heavy greedy workloads).
- **Findings backlog: 47 open items** (C1–C47: **16 med · 24 low · 7 nit**), none
  applied in Phase 2. The top three (C1–C3) are the determinism agent's MEDs — all about
  *measurement/state integrity in persistent processes*, not solver correctness.
- **3 candidate steps dropped by measurement**: E2 (ordered-merge candidates walk,
  **+55.4%** `search_inf` — the census cost model was inverted), E4 (staged compares,
  **+0.21%**, and a GC probe showed there was no GC to save), E8 (the negative control,
  regression by design, reverted).

| open severity | count | ids |
|----------|-------|-----|
| high     | 0     | — |
| med      | 16    | C1–C16 |
| low      | 24    | C17–C40 |
| nit      | 7     | C41–C47 |

---

## 1. The headline — are our cuts appropriate, and does removing them help performance?

The short answer: **the cuts that remain are appropriate; the ones that weren't are now
gone; and "removing cuts helps" is true only under a precise, now probe-and-benchmark-
documented condition.** The long answer, in seven steps:

### (a) Cut changes are invisible to the inference ratchet — by construction

SWI's inference counter counts **call ports only**. Executing a `!` is not an inference,
and — less obviously — **retrying a dead choicepoint's alternative clauses on
backtracking is not an inference either** (the failed head unifications never reach a
call port). So a transformation that only adds/removes cuts or choicepoints cannot move
`search_inf` at all. X2 proved this bit-identically across **12 rungs**: after five
cut-deletion/head-restructure steps on the hot arrange path, every rung's `search_inf`
was *exactly* unchanged —

| rung | search_inf (before = after) |
|---|---:|
| ladder_09x09_08w | 28,212 |
| ladder_15x15_12w | 101,588 |
| ladder_21x21_25w | 302,257 |
| ladder_15x15_28w | 374,521 |
| ladder_15x15_32w | 969,070 |
| ladder_09x09_16w (heavy) | 674,961 |
| ladder_21x21_80w (heavy) | 4,487,855 |
| ladder_15x15_34w (heavy) | 13,650,369 |
| ladder_15x15_36w (heavy) | 38,275,505 |
| ladder_09x09_17w (heavy) | 38,533,195 |
| ladder_15x15_40w (heavy) | 10,411,252 |
| ladder_21x21_82w (heavy) | 6,822,207 |

Δ = 0.00% on all twelve, up to 38.5M inferences. Consequently a cut's real cost/benefit
lives in exactly two places the ratchet cannot see: **wall time** (choicepoint
creation/scan constants) and **choicepoint-stack memory**. Campaigns gating on
`search_inf` must treat CP-elimination steps as wall/RSS candidates, not inference
candidates — inference counts only move when *counted work* changes (see (e)).

### (b) Where indexing recovers determinism, cut deletion is free — and eliminating CP leaks pays in memory

The census's probe-verified boundary (the load-bearing empirical fact of this audit,
`jitindex.md` §2.17/§2.17.3 + live probes): on SWI 10.1.10, cut-free deterministic clause
selection is reliable **only**

1. on the **first** argument (constants and/or functors — e.g. `none` vs `state(_,_)`,
   `unbounded` vs `2`), or
2. via the **two-clause `[]`/`[_|_]` special case** — on *any* argument, including a
   leading body unification hoisted into the head by §2.17.3 body-code indexing.

A NON-first-argument "primary index argument" linear scan did **not** prune the
choicepoint when the matching clause is not last — *even after a JIT hash existed on that
argument* (probe `a3`/`p1`: arg-4 constant dispatch, CP survives after 50 calls with
`jiti_list`-confirmed hash). The manual alone would not have predicted this; the probe
did, and X4's production-predicate results confirmed it.

Inside that boundary, the cut is provably redundant and deleting it is
semantics-AND-inference neutral (X1/X2/X4: every such step had byte-identical goldens
and Δinf = 0, with post-edit `deterministic/1` probes flipping nondet→**DET** at every
site). And where the old shape *leaked* choicepoints, deletion buys real memory:
the `remove_x/3` rewrite (X2 step 5 — old clause order left up to ~|Words| live CPs per
search node stacked above `member/2`'s real alternative) reduced whole-process RSS by
**−6.3% and −6.5%** on the two 21×21 heavy rungs (80/82-word lists), reproduced across
both interleaved A/B pairs; all other rungs ±1% noise. It also unblocks last-call
optimization in the traversal.

### (c) The negative control: where the cut IS the pruning, deletion measurably loses

X6 Part A (census E8) deleted **only** the cut in `remove_slot/4` (fill.pl:336) — a
predicate whose dispatch lives in **body `==` guards over identical `[_|_]` heads**, i.e.
outside the indexable boundary in (b). The census probe had predicted: solution set
unchanged (semantics-safe), but one live CP per matched node survives until the top-level
commit. Both halves materialized:

| rung | search_inf base | E8 | Δ% |
|---|---:|---:|---:|
| sq04_full | 371,343 | 371,492 | +0.04% |
| g11_full | 514,725 | 515,744 | +0.20% |
| g11_full_seed | 419,582 | 420,470 | +0.21% |
| sq05_full | 876,182 | 876,376 | +0.02% |
| g17_full | 1,822,010 | 1,833,328 | **+0.62%** |
| g21_full | 1,853,583 | 1,860,528 | +0.37% |
| g13_full | 2,934,984 | 2,943,213 | +0.28% |

Regression on **7/7 rungs** (+0.02%…+0.62%), deltas scaling with backtracking volume (the
three deepest-backtracking rungs g17/g21/g13 carry the largest absolute deltas:
+11,318 / +6,945 / +8,229). `make bench-fill-check`: **FAIL (1 regression, 0 wins)** —
the project's own gate rejected the "cleanup". Mechanism: every later backtrack past a
matched node re-enters clause 2 and rescans the rest of the counted-slot list to a
failing `[]` before reaching the search's *real* alternative. (Unlike (a)'s dead-CP
retries, these rescans are real calls — hence counted.) Fully reverted; cuts of this
class — `remove_slot/4` is the canonical shape — **must stay**. The control also
demonstrates the audit wasn't confirmation-biased: the campaign design included, ran,
and published its own refutation case.

### (d) The legacy `(A ; B), !` pattern vs if-then-else — VM-verified

`vm_list/1` disassembly: `( Disj1 ; Disj2 ), !` compiles to a **real choicepoint**
(`c_or` + `i_cut`); the equivalent `( Cond -> Then ; Else )` compiles to
`c_fastcond`/`c_fastcut` (no choicepoint at all). On the isolated pattern the ITE form
measured ~11% faster (0.331 s vs 0.371 s over 5M iterations) with **identical inference
counts per call** (2.00/3.00 both forms) — again invisible to the ratchet. X1 converted
the four hot-path sites (`check_letters/6`, `assign_letters/7`, `check_prev_cell/4`,
`check_next_cell/4`); post-edit `vm_list` confirms no `c_or`+`i_cut` remains at any of
them, `search_inf` bit-identical on all 5 rungs, goldens byte-identical. The
`check_prev/next_cell` pair is additionally a robustness win: their branches genuinely
**overlap** (a row-start cell whose numeric predecessor is empty can succeed both ways),
so the old cut was load-bearing — the ITE reproduces the first-solution commit while
structurally removing the duplicate-success hazard instead of papering over it.

X1's interleaved A/B wall measurements leaned the right way (full stack best-vs-best
total **−2.85%**, the two largest rungs −0.29%/−4.19%) but the per-rung spread
(−4.2%…+4.5%) straddled zero under machine load 3–4, and X2's same-code spread reached
14% — so per the wall-noise policy no wall number is quoted from those runs. Idle-machine
serial A/B results for the arrange search + greedy path:

Measured 2026-07-06 on the otherwise-idle bench host (`lana x86_64`), three interleaved
base↔treated pairs (`git apply -R`/`git apply` of the combined 18-step patch between
runs), min-of-3 per rung per state, fresh process per run. Treated = **all** applied
steps, so strict-path deltas conflate X1's ITE conversions with X2's choicepoint work
and X6's `once/1` drops (per-step wall attribution was not re-measured serially).

Strict search (in-process `arrange_best_layout/6`, the `search` bucket of
`benchmarks/run_arrange.pl`):

| rung                  | base ms | treated ms | Δ wall  |
|-----------------------|--------:|-----------:|---------|
| ladder_09x09_08w      |   2.310 |      2.296 | −0.61%  |
| ladder_15x15_12w      |   7.200 |      7.169 | −0.43%  |
| ladder_15x15_28w      |  28.684 |     28.430 | −0.89%  |
| ladder_21x21_25w      |  22.461 |     22.063 | −1.77%  |
| ladder_15x15_32w      |  65.094 |     63.616 | −2.27%  |

Greedy path (fresh-process X3 harness, `statistics(walltime)` around the in-process
greedy APIs):

| workload       | base ms | treated ms | Δ wall   |
|----------------|--------:|-----------:|----------|
| cand17_strict  |    22.0 |       18.0 | −18.2%   |
| be17           |    22.0 |       18.0 | −18.2%   |
| be15_32        |   378.0 |      329.0 | −12.96%  |
| be21_80        |  1747.0 |     1519.0 | −13.05%  |

Reading: the strict-path deltas are individually small, but all five rungs move the
same direction and the two largest search volumes move most (−1.77%/−2.27%) — consistent
with a real, modest win rather than noise (X2's same-code spread on a *loaded* box was
±14%; on the idle box repeat-run spread was well under these deltas). The greedy deltas
are unambiguous and reproduce X3's under-load measurement (−17…−18% there, −13…−18%
here). Between-state inference drift was the expected few counts per strict rung (X6's
`once/1` drops, −2/rung) — confirming the wall delta is not hiding an algorithm change.

### (e) The findall/copy-term tax beats the cut tax

Where inference counts *and* wall both moved, the mechanism was never cuts — it was
**copying and generator machinery**, i.e. counted/C-level work:

- **X3 (greedy double-findall elimination):** the greedy constructor snapshotted the
  whole `gs/2` grid bundle (2×N² args; 882 args each at 21×21) through findall *per legal
  placement per word per greedy step*, then re-copied each word's best grid in the outer
  findall. Rewriting to a small ground descriptor + pure `check_word_fits/5` legality
  probe + one `assign_word/10` for the winner: **wall −17% / −18%** on the heavy greedy
  workloads (be15_32: 477/495 → 396/420 ms; be21_80: 2324/2367 → 1909/2117 ms) with
  inference deltas of only **−2.16% / −1.53%** — findall's copy-term work is C-level and
  never appeared in inference counts. Goldens (including the candidates golden, the one
  golden riding the greedy path) byte-identical on the first attempt; strict rungs zero
  drift.
- **X5-E1 (findall→first-order walk):** replacing
  `findall(P-V, (nth0(P,Vars,V), nonvar(V)), Bound)` in both recount kernels with a
  first-order positional accumulator moved the *gated* metric down on **every** rung:
  −0.48%…−1.12%, core total **−0.78%** (−68,881 inf), heavy g15_full −0.52%. Modest —
  the big findall overheads had already been paid down by the earlier F-H1/F-L1
  campaign — but strictly down everywhere, with no downside.

### (f) Two census predictions falsified by measurement

The census's cost models were wrong twice, and only measurement caught it:

- **E2 (ordered-merge candidate materialization): +55.4% `search_inf`** core total
  (worst rung +72.6%), semantics fine, performance inverted. Two reasons: SWI's
  det-mode `nth0/3` is **3-way unrolled** (~I/3 inferences per index-I lookup, not O(I)),
  and `candidates/4` only ever materializes the **MRV winner** — a small, sparse index
  set inside a length bucket running to ~24k words, so Σ I_k/3 ≪ a full prefix walk.
  The O(Σidx) → O(bucket) framing was exactly backwards for MRV-selected slots.
- **E4 (staged compares replacing `c/3`-key `@=<`): +0.21%** on the gated metric on
  every rung — and the predicted "heap churn/GC" benefit does not exist: a
  `statistics(garbage_collection)` probe over whole searches (g13 + g21, 3 samples each,
  both versions) showed **0 collections, 0 bytes freed, 0 ms GC** under both. The two
  `c/3` cells per comparison are bump-pointer allocations that never trigger a
  collection at this scale, while the staged ITE chain costs more *counted* inferences
  than one term-build + one `@=<` builtin call.

Both reverted (patches preserved as `*_REVERTED.patch` in the campaign archive). The
lesson is the audit's method in one line: **measure, don't assume** — the same discipline
that would have rejected E2/E4 is the one that validated everything in (b)–(e).

### (g) Bottom line — taste vs pragmatism

The purity improvements were free or profitable **exactly** where clause indexing could
take over the cut's job (first-arg dispatch, `[]`/`[_|_]` heads, hoisted body
unifications) — there, deleting 14 cuts and 9 `once/1` wrappers cost nothing on any gate,
made determinism *structural* (provable from head shapes rather than dependent on a
control operator), and delivered the campaign's only reproducible resource wins (−6% RSS,
−17…−18% greedy wall, −0.78…−2.9% fill `search_inf` where counted work shrank). They were
**harmful** exactly where indexing could not dispatch — the E8 class — and the remaining
28 cuts each have a recorded reason to exist (§[2](#2-cut--impurity-census-aggregate)).
The empirical boundary is now documented, probe-verified, and calibrated by a published
negative control. Good taste and pragmatism turn out to agree here, provided "remove the
cut" is always spelled "restructure until indexing carries the determinism, *then*
remove the cut — and keep it where it is doing pruning work indexing cannot".

---

## 2. Cut / impurity census (aggregate)

Full site-by-site tables live in the two Phase-1 census reports; this section is the
stand-alone aggregate. Classification legend — **GREEN**: clauses already mutually
exclusive; removal preserves the solution set. **RED**: removal changes observable
behaviour (the cut is the predicate's contract). **RED-BUT-CONVERTIBLE (RBC)**: red only
because of clause ordering / catch-alls; a guard or head restructure would green it.

### 2.1 Cuts — 42 total at baseline

| class | core+arrange (24) | fill+lint+export+stockgrid+metrics+browser (18) | total |
|---|---:|---:|---:|
| GREEN | 10 | 14 | **24** |
| RED — deliberate contract, keep | 2 | 0 | **2** |
| RED-BUT-CONVERTIBLE | 12 | 4 | **16** |

Notably: **no cut in the fill/lint/export/stockgrid/metrics/browser group changes the
first-solution semantics of its callers** — every red case there is about extra solutions
on backtracking that a guard could eliminate. `export.pl` and `solve_browser.pl` contain
**zero cuts**.

**Removed by the applied Phase-2 steps (14):**

| site | predicate | via |
|---|---|---|
| core.pl:281 | `find_crossword/6` (RBC C5) | X6.B5 — enumerated strategy clauses |
| core.pl:383 | `mrv_count/8` (GREEN) | X2.2 |
| core.pl:532 | `select_inc/10` (RBC C9, hottest predicate) | X2.3 — `[_|_]` head |
| core.pl:561 | `inc_counts/8` (GREEN) | X2.1 |
| core.pl:768 | `check_letters/6` (GREEN) | X1.1 — ITE |
| core.pl:783 | `check_prev_cell/4` (RBC, overlapping branches) | X1.3 — ITE |
| core.pl:795 | `check_next_cell/4` (RBC) | X1.3 — ITE |
| core.pl:821 | `assign_letters/7` (GREEN) | X1.2 — ITE |
| core.pl:1165 | `remove_x/3` (RBC-in-theory) | X2.5 — det rewrite |
| fill.pl:211 | `candidate_count/5` (GREEN, arg-4 → arg-1) | X4.3 |
| fill.pl:302 | `fill_search_inc/5` (GREEN) | X4.1 — `[]`/`[_|_]` heads |
| fill.pl:362 | `shares_cell/2` (GREEN, predicate deleted) | X4.2 — `ord_intersect/2` |
| lint.pl:111 | `lint_atom/2` cl.2 (dead cut) | X6.B2 |
| lint.pl:248 | `connected/2` (redundant under §2.17.3) | X6.B2 |

**Remaining (28) and why each stays:**

- **First-solution / existence contracts — KEEP (2):** `shares_letter/2` (core.pl:604 —
  semidet existence check; the cut *is* the predicate's contract) and `clashing_cell/6`
  (arrange.pl:611 — "report the FIRST clashing cell" lives in this cut; error path only).
- **Perf-load-bearing, non-indexable dispatch — KEEP (the E8 class):**
  `remove_slot/4` (fill.pl:336 — **negative-control-proven**, §1c) and `grab_run/4`
  (stockgrid.pl:109 — catch-all clause; removal yields every run prefix as an extra
  solution → wrong slots; cold load path). Also `reconcile_answer/3` (arrange.pl:617 —
  GREEN in context, semidet early-exit; ≤1 match guaranteed by `check_unique_answers/1`)
  and fill's two KEEP/COLD from the census verdicts: `remove_slot` above plus
  `meta_value_or/4` (fill.pl:670 — cold RBC; disappears entirely if C7's
  `library(option)` swap lands).
- **Cold validate/dispatch RBC — keep as-is or convert opportunistically (→ C47):**
  `require_strategy/1` (core.pl:116), `with_output('',…)` (core.pl:130),
  `load_clues_by_extension/3` ×2 (core.pl:164,167), `select_word/9` mrv-seed clause
  (core.pl:335 — warm but only on the non-default `mrv`/`mrv_capped` strategies; the
  census's probe-verified `[_|_]` restructure is available),
  `reconcile_fragment_size/3` ×2 (arrange.pl:544,545 — couples with C18's steadfastness
  fix), `symmetry_override/4` (lint.pl:164 — removal would duplicate report entries
  inside a findall), `double_unch_end/2` (lint.pl:228 — call-site green).
- **GREEN but harmless/unindexable — keep (cold):** `set_search_seed(-1)`
  (core.pl:486) and `set_check_target(-1)` (arrange.pl:111) — sentinel clauses,
  guard-exclusive, dead CP only; `fragment_atom/2` ×2 (arrange.pl:522,523 —
  type-disjoint guards; cl.2's cut is a no-op); `rows_of/3` (core.pl:982 — cl.2 provably
  fails for `[]` at valid grid sizes); `lint_atom/2` cl.1 (lint.pl:110 — type guards
  can't index); `barred_max_unch/2` ×6 (lint.pl:219–224 — pairwise-disjoint guards;
  P-audit Appendix A records the table as intentional; cosmetic only);
  `word_cells/5` (metrics.pl:61 — the discriminating argument is an *output*, so no head
  restructure can index it).

### 2.2 once/1 — 12 total at baseline, 9 deleted

Deleted (all X6): the five `once(assign_clue_numbers(...))` wraps at
arrange.pl:240,354,666,698,876 (made unnecessary by making the goal deterministic *at
source*, X6.B4) and the four redundant `once/1`s inside if-then-else conditions —
arrange.pl:277,678 (X6.B3, −2 inf on all 5 arrange rungs) and fill.pl:421,437 (X6.B1,
−1 inf on all 7 fill rungs; the arrow was already the commit point).
Remaining (3): `fill.pl:378` `once(select(...))` — **keep** (the off-hot-path reference
selector; its determinism is pinned by the P13 plunit test); `fill.pl:392` and
`stockgrid.pl:79` `once(assign_clue_numbers(...))` — **now redundant post-X6.B4** (the
goal is det); open cleanup item **C23**.

### 2.3 Everything else impure — inventory

- **Soft-cut `*->`: zero occurrences anywhere. `bagof`/`setof`: zero anywhere.**
- **`\+`:** sparse and deliberate. Verified live in the current tree: core.pl has 2 real
  sites — the `\+`-driven loop inside `count_upto2/2` (core.pl:442; the manual's own
  `succeeds_n_times` idiom) and `no_word_merge_bg/2` (core.pl:729, the maximality
  invariant, negation-as-absence over the boundary grid); arrange.pl has 4, all on cold
  filter/report paths (:174, :321, :623, :863); fill 1 (the hot `\+ memberchk(Word,
  Used)` at fill.pl:309 — recorded keep, see F15 in the KEEP list); lint 2 (BFS/parity,
  cold). *(An earlier campaign summary claimed "zero `\+` in core/arrange"; that is
  incorrect — the counts above are grep-verified. Flagged per the no-smoothing rule.)*
- **If-then-else:** ~20 constructs in core, ~34 in arrange, 77 constructs (80 arrows)
  across fill/lint/export/stockgrid/metrics/browser. Dominant patterns: validate-or-throw
  guards (~20 sites, cold, idiomatic), default-on-missing assoc/dict lookups (the hot
  ones are semidet lookups where ITE is the right tool — no head dispatch is possible on
  a runtime key), outcome-tag chains (cold). No bare `(Cond -> Then)` without else. No
  condition with side effects.
- **findall:** core 3 sites / arrange 12 / fill 7 / lint 8 / export 1 / stockgrid 7 /
  metrics 2. The two that were hot (fill.pl:170/222) were rewritten by X5-E1; the greedy
  double-findall by X3. `pair_crossings/3`'s findall is golden-frozen ("ORDER IS
  LOAD-BEARING") — forbidden territory. Two sites deliberately work *around* findall's
  copy-term semantics (`map_list_to_pairs` identity preservation; ground-answer-keyed
  `selectchk`) — preserve in any rewrite.
- **aggregate_all / forall:** `aggregate_all(count, …)` at core ×2, arrange ×2, lint ×1
  — all the right tool; `forall` ×5 (arrange 1, lint 2, stockgrid 2), all proper
  test-iteration uses.
- **Dynamic state — the complete inventory:** 3 config predicates
  (`verbose_mode/0` core.pl:142; `search_seed/1` core.pl:479;
  `check_target_override/1` arrange.pl:101 — all single-writer,
  retractall-before-assertz, `-1`/`false` sentinel-cleared) + 2 tabled pure memos
  (`answer_letters/2` core.pl:595; `pair_crossings/3` core.pl:646 — abolished per
  `mrv_inc` search entry only, → C1/C2) + 1 non-backtrackable local counter
  (`count_upto2/2`'s `nb_setarg/3` on a fresh `counter(0)` holder — reentrant,
  manual-sanctioned). **Zero** dynamics, globals (`nb_setval`/`b_setval`), or
  assert/retract in fill/lint/export/stockgrid/metrics/solve_browser.
  `catch/3`-as-control: one site (fill.pl:641, artifact-read error-translation shim →
  C33). Budget control uses `call_with_inference_limit/3` with reified outcomes and
  explicit why-no-catch comments (the P2 remediation, preserved).
- **var/nonvar as representation:** the grid protocol (empty cell ≡ unbound variable,
  bindings undone by the trail) makes `var/1`/`nonvar/1` load-bearing at ~16 core sites
  plus arrange/fill reads. This is the file's designed-in, documented impurity — not a
  removal candidate; every transformation preserved the "a check never binds a cell"
  invariants (core.pl:731–749, 826–828), and X1/X3 were verified against exactly those.

---

## 3. Applied remediation — the Phase-2 APPLY steps

All steps below are **applied to this worktree (staged, uncommitted)**. Merged-state
gates re-run green: `./run_tests.sh` ALL TESTS PASSED · `make bench-check` PASS (no
change, 0 regressions) · `make bench-fill-check` PASS (6 WINs, 0 regressions; per-rung
deltas in the [Executive summary](#executive-summary)). **Owed after commit:
`make bench-fill-record`** to ratchet the fill baseline down to the new numbers.

<a id="applied-step-count-note"></a>
> **Count note (flagged honestly):** the campaign tallied "17 APPLY steps"; the
> enumeration below contains **18** discrete experiment steps. The X4.3 and X5.1 patches
> overlap in `candidate_count/5` and were merged via a git 3-way merge into a single
> applied patch — 17 distinct patches, 18 semantic changes. Both counts are correct for
> what they count.

| # | Step | Site(s) | Measured effect | Commit-style one-liner |
|---|---|---|---|---|
| 1 | X1.1 | `check_letters/6` core.pl:761–769 | Δinf **0** (5 rungs); VM `c_or`+`i_cut` → `c_fastcond` | perf(core): check_letters disjunction+cut → if-then-else; CP-free inner counting loop |
| 2 | X1.2 | `assign_letters/7` core.pl:807–822 | Δinf **0**; binding kept in else branch | perf(core): assign_letters disjunction+cut → if-then-else (check-never-binds invariant preserved) |
| 3 | X1.3 | `check_prev_cell/4` + `check_next_cell/4` core.pl:774–794 | Δinf **0**; kills latent duplicate-success hazard (overlapping branches) | fix(core): prev/next-cell `(A;B),!` → `(A->true;B)`; first-solution commit made structural |
| 4 | X2.1 | `inc_counts/8` core.pl:560 | Δinf **0**; probe DET (arg-1 `none`/`state(_,_)`) | refactor(core): delete redundant cut — first-arg indexing carries inc_counts dispatch |
| 5 | X2.2 | `mrv_count/8` core.pl:382 | Δinf **0**; probe DET (`unbounded`/`2`) | refactor(core): delete redundant mrv_count cut + factually wrong comment defending it |
| 6 | X2.3 | `select_inc/10` core.pl:531 | Δinf **0**; seed-solution multiset verified unchanged | refactor(core): select_inc cl.2 `[_|_]` head enables []/[_\|_] indexing; RED cut deleted from hottest predicate |
| 7 | X2.4 | `assign_words_inc/9` + `assign_words/9` core.pl:458/304 | Δinf **0**; terminal-`[]` CP probe nondet→DET | refactor(core): cons heads on the search drivers — solutions exit without a live CP |
| 8 | X2.5 | `remove_x/3` core.pl:1162 | Δinf **0**; **RSS −6.3%/−6.5%** on 21×21 heavy rungs (both A/B pairs) | perf(core): deterministic ==-based remove_x — kills up-to-\|Words\| leaked CPs per node, unblocks LCO |
| 9 | X3 | `next_move/5`+`word_best_placement/8` arrange.pl:1008–1046 (+`check_word_fits/5` export, core.pl) | greedy inf −0.97…−2.16%; **wall −17%/−18%** on be15_32/be21_80; strict rungs zero drift | perf(arrange): greedy stops snapshotting grids through the double findall — probe legality, assign winner once |
| 10 | X4.1 | `fill_search_inc/5` fill.pl:302 | Δinf **0**; probe DET; spine cut-free | refactor(fill): []/[_\|_] heads on the search spine; base-case cut deleted |
| 11 | X4.2 | `shares_cell/2` → `ord_intersect/2` fill.pl:355–363 | **search_inf −2.17%** net (−190,783; best g17 −3.84%, worst sq05 +0.14% in-tolerance); ratchet PASS 5 improvements | perf(fill): recount-decision via ord_intersect merge — O(n+m), CP-free; predicate+cut deleted |
| 12 | X4.3 | `candidate_count/5` Masks → arg 1, fill.pl:203–226 (+callers, tests lockstep) | Δinf **0**; probe DET all paths; hottest predicate CP-free | refactor(fill): first-arg none/masks(_) dispatch replaces cut on the count kernel *(3-way-merged with #13)* |
| 13 | X5.1 | `bound_positions/2` replaces findall/nth0 at fill.pl:170+222 | **search_inf −0.78%** core total, −0.48…−1.12% every rung; g15 heavy −0.52%; ratchet PASS 6 improvements | perf(fill): first-order bound-positions walk replaces the hot findall/nth0 kernel |
| 14 | X6.B1 | `fill_attempt/8`+`fill_attempt_masked/9` fill.pl:421/437 | exactly **−1 inf ×7** fill rungs | refactor(fill): drop once/1 inside ITE conditions — the arrow is the commit point |
| 15 | X6.B2 | `lint_atom/2` cl.2 + `connected/2` lint.pl:111/248 | cold; post-edit probes det=true | refactor(lint): delete a dead cut and a §2.17.3-redundant cut |
| 16 | X6.B3 | `construct_one/7`+`construct_fragment_one/6` arrange.pl:277/678 | exactly **−2 inf ×5** arrange rungs | refactor(arrange): drop redundant once/1 inside ITE conditions |
| 17 | X6.B4 | `add_clue_nums/3` restructure core.pl:886–899 + 5 once-wraps deleted (arrange.pl:240,354,666,698,876) | det probe FALSE→true; numbering `variant_sha1` identical; +≤0.02% emit path, both gates PASS | fix(core): clue-number grouping dispatches on []/[_\|_] — assign_clue_numbers det at source, five defensive wrappers deleted |
| 18 | X6.B5 | `find_crossword/6` core.pl:280–299 | **zero** drift ×5 rungs; `all_crossword` counts verified identical for all 4 strategies, fixed + unbound Loc (522/2234) | refactor(core): enumerated strategy clauses replace the var catch-all; RED-convertible cut deleted |

Also resolved *by* these steps (no separate backlog entry): determinism-agent **L2**
(`assign_clue_numbers/2` exported-nondet — fixed at source by #17), **L11**
(`crossword/4` unwrapped call — moot once the goal is det), **L1** (`remove_x` CP leak —
#8), style-core **L8** (`remove_x` layout/self-description — #8's rewrite), and the
core-census once-hygiene items (#14/#16/#17).

**Dropped by measurement (patches preserved, not applied):**

| Candidate | Site | Result | Disposition |
|---|---|---|---|
| X5-E2 ordered-merge `candidates/4` | fill.pl:177–183 | **+55.4%** search_inf (worst rung +72.6%) — census cost model inverted (§1f) | DROP, reverted |
| X5-E4 staged compares `count_le/2` | fill.pl:325/330 | **+0.21%** search_inf; GC probe: 0 collections/0 bytes/0 ms both versions — no GC to save | DROP, reverted |
| X6-E8 delete `remove_slot/4` cut | fill.pl:336 | **negative control**: +0.02…+0.62% on 7/7 rungs, ratchet FAIL (§1c) | reverted by design — datapoint only |

---

## 4. Findings backlog (C1–C47) — found, verified, NOT applied

Deduplicated across all eight Phase-1 reports; severity-ordered. Everything here is
`status: open`. Sources are noted as (census / style-core / style-arrange / style-fill /
style-small / stdlib / determinism).

### Med

#### C1 — Table abolition is strategy-local → run-order-dependent inference counts
`core.pl:280–299` (mrv_inc clause abolishes; `:295-299` generic clause does not) ·
behaviour preserving (measurement integrity, not output) · flags: perf-relevant (bench
validity) · **open** · (determinism M1)

The comment at core.pl:283–289 says the `abolish_table_subgoals` calls exist so "each
search's inference count [is] self-contained and order-independent" — but only the
`mrv_inc` clause abolishes; `baseline`/`mrv`/`mrv_capped` run with whatever tables are
warm. **Probe:** bundled-17 fixture, `find_crossword(baseline, …)` twice in one process:
run 1 **11,454** inferences (cold), run 2 **4,237** (warm) — 2.7×; after an interleaved
mrv_inc run a "first" baseline measures 4,422. Any recorded baseline/mrv/mrv_capped
inference number (the `docs/experiments.md` evidence trail, `run_matrix.pl`'s "portable
metric") is a function of harness run order. **Fix:** hoist the two abolishes into a
shared prologue of `find_crossword/6` (all clauses) or into the benchmark sampler; the
rebuild is ~|Words|² pairs (the code already argues this is cheap). Note X6.B5's
enumerated-clause restructure makes the shared-prologue fix mechanically easy now.

#### C2 — Unbounded table growth across greedy-path solves in a persistent process
`core.pl:646–654` (`:- table pair_crossings/3`), `:595–599` (`answer_letters/2`); greedy
consumers arrange.pl:965–1046 with no abolish on best-effort/candidates/fragment paths ·
behaviour preserving · flags: **WASM** · **open** · (determinism M2)

Only the strict `mrv_inc` path abolishes the tables. A long-lived process serving
repeated **best-effort** solves — exactly what `solve_browser_json` with
`bestEffort:true` exposes — accumulates one table variant per distinct
(Letters, PLetters) pair ever seen, forever. **Probe:** `pair_crossings` variant count
0 → **30** after one best-effort solve (17 words) → **42** after a second solve with 4
new words → drops to 15 only when a strict solve happens to run → climbs again (27).
Also makes greedy-path inference counts call-history-dependent (why X3 had to measure in
fresh processes) — a latent footgun if greedy rungs are ever added to the ratchet.
**Fix:** abolish at every solver entry seam, or add a `reset` predicate the WASM dispatch
calls per request (pairs with C3).

#### C3 — Module-global mutable state persists across in-process solves; browser wrapper has no reset
`core.pl:479–496` (`search_seed/1`), `:142–148` (`verbose_mode/0`), arrange.pl:101–115
(`check_target_override/1`); the owed reset is self-documented at
`wasm/client/solve_browser.pl:71–75` · behaviour preserving today (hazard demonstrated,
not reachable from the current browser input schema) · flags: **WASM** · **open** ·
(determinism M3)

Set-and-forget module globals: any stray write poisons every later in-process solve
until explicitly cleared, and there is no per-request reset in the browser dispatch —
only a comment promising one before seed/checkTarget/verbose keys are wired. **Probe:**
`solve_browser_json/2` twice → byte-identical; after `set_search_seed(42)` → output
CHANGED (layout + `diagnostics.arrange.seed` appear in subsequent "default" solves);
after `set_check_target(1)` → output CHANGED; `set_search_seed(-1)` verified to restore
byte-identical defaults. Writer hygiene is otherwise good (single-writer setters,
retractall-before-assertz, sentinels). **Fix:** as the code itself prescribes —
unconditional `set_search_seed(-1), set_check_target(-1), set_verbose(false)` (plus C2's
abolish) at the top of the browser entry, and a same-worker determinism test.

#### C4 — No PlDoc (`%!`) structured comments or determinism annotations anywhere
whole project (`%!` appears nowhere under `prolog/`; e.g. core.pl:16–66 export list,
arrange.pl:20–30, fill exported API, metrics' 10 exports) · behaviour preserving
(comments only) · **open** · (style-core M1 / style-arrange S2 / style-fill F2 /
style-small M1 — all four lanes independently)

The prose comments are unusually good — often better than typical PlDoc — so this is a
*convention decision*, not sloppiness. What is genuinely lost: machine-readable
mode/type/determinism lines (`%! name(+Arg) is det.`, per
`swi-manual/packages/pldoc.md`) for the cross-module API, IDE/cross-referencer support,
and checkable determinism claims that today live in scattered prose ("fails on anything
else", "always succeeds with a Res tag"). Several findings in this series (L2-class
det-leaks, silent-fail exports) are exactly the class det-annotations surface. **Fix:**
adopt `%!` headers for exported predicates first (content mostly exists in prose
already) — **or** record the plain-comment convention in AGENTS.md so future audits stop
re-raising it.

#### C5 — Explicit-import (P11) discipline is inconsistent; `load.pl` dies under `autoload(false)`
`load.pl:12` (`directory_file_path/3` ← library(filesex), autoload-only — verified: under
`set_prolog_flag(autoload, false)` the alias directive errors and **nothing loads**);
core.pl:84–92 (only http/json, random, aggregate imported — assoc/apply/lists/pairs all
autoload); arrange.pl:32–37 (lists/assoc/ordsets/http-json autoload — while its own pairs
import carries the P11 rationale comment); metrics.pl (lists autoload); lint.pl:20–23
(lists autoload); export.pl (pairs autoload at :142, and its header claims "Zero project
dependencies (it only needs JSON + lists)" — missing pairs) · behaviour preserving ·
flags: WASM-adjacent (qsave deployment) · **open** · (determinism L5 / style-core M3 /
style-arrange S1 / style-small L1+L4)

This is the P11 finding's follow-through: P11 (closed) added `library(pairs)` to
quality.pl; the *rationale* ("survive `qsave_program(..., [autoload(false)])`") is quoted
in arrange.pl/metrics.pl/fill.pl comments but applied to only a few imports. Manual:
`swi-manual/autoload.md:58` — a saved-state build must load autoloaded predicates
explicitly. **Fix:** add the missing `use_module` lines per the determinism agent's
imports probe (core: apply/assoc/lists/pairs; arrange: lists/assoc/ordsets + JSON;
metrics/lint: lists; export: pairs; load.pl: filesex) — or drop the P11 rationale
comments and accept autoload uniformly. Either way, make it consistent.

#### C6 — Legacy `library(http/json)` alias where 10.1.10 has core `library(json)`
solve_browser.pl:40 (strongest case); also lint.pl:20, export.pl:24, stockgrid.pl:31,
core.pl, fill.pl (project-wide pattern) · behaviour preserving (load-time only) · flags:
**WASM** · **open** · (style-small M4; `check/0` also emits "Library was moved:
library(http/json) --> library(json)")

In the **WASM image the `http/json` alias does not resolve at all** — the project's own
`wasm/README.md:165–167` documents the resulting "source_sink library(http/json) does not
exist" load noise, working only because the predicates autoload from `library(json)`;
solve_browser.pl:93–95 re-documents the same gotcha in prose while line 40 keeps the
directive that causes it. Probe-verified: `use_module(library(json))` loads cleanly and
`json_read_dict/3`, `json_write_dict/2`, `atom_json_dict/3` are `imported_from(json)`
(`swi-manual/packages/json.md:56`). **Fix:** `:- use_module(library(json)).` — at minimum
in solve_browser.pl (the file whose environment actually breaks the alias), ideally
project-wide. Coordinate with C5.

#### C7 — `meta_value/3` + `meta_value_or/4` reimplement `library(option)`, via `=..`
fill.pl:669–671 (defs); call sites :651–661; also `memberchk(masks(true), Options)`
at :601 · behaviour preserving · **open** · (stdlib F3 / style-fill F3 — independently
found by both lanes)

The artifact `Meta` is documented as "a keyed list [Key(Value), ...]" (fill.pl:552) —
literally an SWI option list. `meta_value/3` rebuilds `option/2` with a needless
`Probe =.. [Key, Value]` (every call site passes a literal key), `meta_value_or/4`
rebuilds `option/3` (`swi-manual/option.md:24,32`; probe-verified on the exact Meta
shape). **Fix:** `:- use_module(library(option), [option/2, option/3]).`, delete both
helpers and the `=..`; this also deletes the `meta_value_or` RBC cut (§2.1). Cold
artifact-load path; output bytes unaffected.

#### C8 — lint's hand-rolled BFS duplicates `library(ugraphs)`; frontier append is O(n²)
lint.pl:248–264 (`connected/2` + `reach/5`) · behaviour preserving (boolean pass/fail
rule; lint not inference-gated) · **open** · (stdlib F2; subsumes style-small L8's
`append(Q, NewNbs, Q1)` O(n²) frontier and I4's `length/2`-pair vs `Reached == Filled`
note)

~12 lines of queue management (with an O(|Q|) append per node) reimplement reachability.
`vertices_edges_to_ugraph/3` + `reachable/3` (`swi-manual/ugraphs.md:192`,
probe-verified incl. isolated-vertex behaviour) does it in ~6: build the symmetric edge
list from `cell_neighbour/3`, `reachable(Start, G, Reached)`, then `Reached == Filled`
(both ordsets). Needs `:- use_module(library(ugraphs)).` per the C5 discipline.

#### C9 — `solve_browser_str/2` reimplements `atom_json_dict/3`
solve_browser.pl:96–101 · behaviour preserving · flags: WASM · **open** · (stdlib F1)

The `setup_call_cleanup(open_string(...), json_read_dict(..., [default_tag(json)]),
close(...))` dance is exactly what `atom_json_dict/3` does internally
(`swi-manual/packages/json.md:236`; probe-verified with atom and string input including
`default_tag(json)`). The file *already* calls `atom_json_dict/3` itself in
`browser_selftest` (line 145), and its own header confirms it autoloads under the WASM
build. **Fix:** 5 lines → 2, one less stream to leak.

#### C10 — Browser error formals have no `prolog:error_message//1` hooks
solve_browser.pl:105–121 (`browser_missing_size`, `browser_max_mode_unsupported`) ·
behaviour preserving · flags: **WASM** · **open** · (style-small M5)

Unlike every sibling module (lint.pl:321, export.pl:214, stockgrid.pl:165), no message
hook is defined, so an uncaught formal renders as "Unknown error term:
browser_missing_size" (probe-verified) — and the Worker forwards exactly that string to
the page UI (`wasm/client/worker.js` posts `String(err && err.message || err)`). The one
place these errors surface is the place they are least readable. **Fix:** add the two
DCG clauses (pure Prolog, WASM-safe); verify the JS bridge surfaces the printed form.

#### C11 — `stockgrid_report/1` loads and parses the grid file twice
stockgrid.pl:156–161 (calls `stockgrid_load/2`, then `stockgrid_validate_file/2` which
loads again) · behaviour preserving (stderr text identical) · **open** · (style-small M3)

Two opens, two JSON parses; a validation error is thrown from the *second* load. The
single-load composition already exists (`stockgrid_validate/3` takes the loaded term).
**Fix:** load once, destructure, validate the term.

#### C12 — stockgrid re-pairs report entries to source words quadratically and fragilely
stockgrid.pl:136–141 · behaviour preserving (no golden output; Fails order preserved) ·
**open** · (style-small M2)

`member(WR, Ws)` then re-searching the same list with `nth0(I, Ws, WR)` to recover the
index is O(n²) and only correct because no two report dicts can be equal — a future
report-shape change could silently mis-pair cells. **Fix:** enumerate positionally once —
`nth0(I, Ws, WR), nth0(I, Words, SrcW)` inside the findall (`swi-manual/lists.md:68`,
probe-verified in-order generation) — no uniqueness assumption.

#### C13 — Interpreted (non-compile-expanded) yall lambdas on measured paths
arrange.pl:131 (`layout_reward/4`'s `[PW,A,A1]>>(...)` foldl lambda — per placed word per
rescored candidate) and fill.pl:159 (`build_index/2`'s `[K-Is, K-Set]>>list_to_ord_set(Is,
Set)` — per index key group at every dict load) · behaviour preserving · flags:
**perf-relevant** (fill's is on the gated `load_inf` path → needs re-baseline) · **open**
· (determinism L6 / style-arrange S6 / style-fill F4)

Neither module imports `library(yall)`, so `>>` is autoloaded and every call pays a
meta-call + lambda copy — the exact cost mode fill.pl:137–139's own F-L1 comment banishes
from `normalize_word/2` ("~12 inf/char"; `swi-manual/yall.md:52` — lambdas match
hand-written performance only when compile-expanded). Inconsistent with the project's
own measured rationale. **Fix:** first-order named helpers (`add_word_reward/6`,
`pair_ordset/2`), matching the `alpha_chars/2` precedent; fill's needs
`make bench-fill-record`.

#### C14 — `require_strategy/1`'s error message hard-codes the list the registry owns
core.pl:236–237 vs the strategies registry at :98–102 · behaviour preserving
(stderr-only; plunit locks only the `unknown_strategy(S)` term) · **open** ·
(style-core M4)

The registry comment promises "adding a strategy is a one-line entry here", but the
rendered message duplicates the list literally ("expected one of baseline, mrv,
mrv_capped, mrv_inc") — a new strategy silently leaves the diagnostic stale. **Fix:**
render from `strategies/1` at message time via `{}/1` in the DCG
(`swi-manual/printmsg.md`); today's bytes identical.

#### C15 — `assign_word/10` carries a dead `_PlacedWords` parameter in its exported signature
core.pl:665–689 (arg 7 unused in the sole clause body); export at :39; callers
arrange.pl:578,977,1043 + tests all thread a real list into it · behaviour preserving ·
flags: **test-affecting** (arity change; white-box tests call it directly) · **open** ·
(style-core M6)

The boundary-grid refactor made the placed-words argument unnecessary but the arity-10
signature kept it — an exported predicate advertising an input it does not consume is
misleading API documentation (callers may believe legality still depends on the list).
**Fix:** drop to `assign_word/9` (mechanical; 3 arrange call sites + tests + export
list), or if arity is kept deliberately, say so in a comment on the clause head
(minimum fix).

#### C16 — Deprecated `delete/3` in the two letter-normalization helpers
core.pl:596–599 (`answer_letters/2`, tabled) and metrics.pl:56–58 (`word_letters/3`;
+ the entry-shaped head fill must fake at fill.pl:103) · behaviour preserving (probe:
byte-identical on ground char lists) · flags: **perf-relevant** (`word_letters/3` feeds
arrange's gated search/rescore; `answer_letters/2` feeds mrv_inc counting — needs
`make bench-record`) · **open** · (style-core M2 / stdlib F12 / style-fill F7)

`delete/3` is formally deprecated (`swi-manual/lists.md:59–65`) and both helpers do two
passes where one suffices. **Fix:** one pass —
`exclude(separator_char, Cs, Letters)` with `separator_char(' '). separator_char('-').`
(or `subtract(L0, [' ','-'], Letters)`); style-fill F7 additionally suggests giving
metrics a bare-answer `answer_letters/3` core with `word_letters/3` as the entry-shaped
wrapper so fill stops wrapping its atom in a singleton list. Deprecation is the only
real driver; nothing is wrong today. Do opportunistically, with the bench re-record.
(Keep the R1 hyphen semantics — this is the same path revamp R1 fixed.)

### Low

#### C17 — Silent-failure exported APIs (three, batched)
(a) `lint_run/5` unknown profile → plain failure, no diagnosis (lint.pl:128–129;
probe-verified) — CLI pre-validates so no live path hits it, but it is documented
"deliberate API"; fix: `domain_error(lint_profile, Profile)` guard
(`swi-manual/error.md:22`). (b) `export_solve/2` unknown format → silent failure
(export.pl:200–210); **caution:** a catch-all clause would be a *bug* — the existing
clauses also fail by design on malformed cells (the no-partial-output policy); guard
up-front with `must_be(oneof([ipuz, exolve]), Format)` before dispatching
(`swi-manual/error.md:60,87`). (c) `set_search_seed/1` / `set_verbose/1` fail silently
on bad input (core.pl:486–491, 144–145) where the module's error discipline elsewhere is
exemplary; fix: `must_be(nonneg, N)` / `must_be(boolean, V)` — but verify no caller
relies on failure (check tests/arrange.plt seed tests first). · behaviour preserving
(all unreachable from the CLI success path) · **open** · (style-small L3+L5 /
style-core L4)

#### C18 — `reconcile_fragment_size/3` is not steadfast: bound output turns success into a bogus throw
arrange.pl:544–547 · behaviour preserving (misuse path only) · **open** · (determinism L3)

`reconcile_fragment_size(5, none, 6)` **throws**
`fragment_size_mismatch(5, none)` — with the nonsense rendering "gridLength 5 disagrees
with the requested --size none" — where steadfast behaviour is failure. **Fix:** commit
on inputs, unify output after: `reconcile_fragment_size(F, none, Out) :- !, Out = F.`
etc. (also resolves the A4/A5 RBC cuts in §2.1).

#### C19 — `load_clues/2` (.pl path) is not steadfast; JSON path fails cleanly — inconsistent misuse behaviour
core.pl:180–187 (`read_clues_prolog_term/3`: `Term = clues(Words)` unifies the caller's
output while *deciding*) · behaviour preserving · **open** · (determinism L4)

`load_clues('fixtures/bundled_17_clues.pl', bogus)` throws `prolog_no_clues_term(...)`
even though the file contains a clues/1 term. **Fix:** test against a fresh variable,
unify after committing: `( Term = clues(W0) -> Words = W0 ; … )`. The same pattern
exists in `benchmarks/run_matrix.pl:121–125` `read_loop` — fix together (and see C38).

#### C20 — `set_shuffle_seed/1` performs side effects before it can fail
core.pl:504–507 · behaviour preserving (deterministic path never reads the RNG) ·
**open** · (determinism L9)

With a (mis)bound argument, `set_random(seed(random))` reseeds the global RNG from OS
entropy *first*, then `random_between/3` fails — a dirty failure that clobbers RNG
state. **Fix:** `must_be(var, N)` first, or draw into a fresh variable and unify last.

#### C21 — `lint_solve/4` writes the full report to stdout before its output arg can fail
lint.pl:311–317 · behaviour preserving (CLI passes a fresh var) · **open** ·
(determinism L10)

`lint_solve(File, toc, false, 'BOGUS')` fails *after* emitting 8.5 KB of JSON — the
exported seam is not steadfast in the "no output on failure" sense the other solve seams
honour. **Fix:** bind Verdict from the report dict *before* emitting.

#### C22 — Unused / internal-only exports
core.pl:56 `pw_end/2` (zero consumers repo-wide) and core.pl:47 `valid_strategy/1` (used
only by `require_strategy/1` in the same module) · behaviour preserving · **open** ·
(determinism L7)

Contradicts the stated convention ("export list is the verified union of every
consumer's needs"). **Fix:** drop both exports. (`stockgrid_report/1` also has no code
consumer but is documented interactive API — not flagged.)

#### C23 — Two `once(assign_clue_numbers(...))` wraps now redundant post-X6.B4
fill.pl:392 and stockgrid.pl:79 · behaviour preserving · flags: perf-relevant (fill's
deletion shifts its gated rungs by exactly −1 inf → re-record) · **open** · (X6 report
caveat)

X6.B4 made `assign_clue_numbers/2` deterministic at source and deleted arrange's five
wraps; these two sites were out of that step's scope and are now merely defensive.
**Fix:** drop both; run `make bench-fill-record` with the batch.

#### C24 — Stale comments (batch)
core.pl:228–230 ("verified working on the pinned SWI **10.0.2**" — the project runs
10.1.10; hook re-verified working); core.pl:1012–1013 ("keeps the original word *dicts*"
— placed words are `pw/8` compounds now); core.pl:69–81 (header still describes the
working grid as the single `grid/N²` term; the threaded structure is the
`gs(LetterGrid, BoundaryGrid)` bundle documented only at :658–664); fill.pl:14 (header
says "Exports: fill_solve/4 … nothing else" but the module exports four predicates —
`fill_solve/4, fill_solve_index/5, fill_save_index/2,3`) · behaviour preserving
(comments only) · **open** · (style-core L1 / style-fill F1)

#### C25 — solve_browser is not a module and reaches unexported internals (acknowledged spike debt)
solve_browser.pl:50,53–54,82 (`crosswordsmith_core:doc_to_words/2`,
`crosswordsmith_arrange:arrange_best_layout/5`, `…:arrange_diag_layout_dict/5`); also
`benchmarks/subjects.pl:57` reaches `crosswordsmith_arrange:arrange_best_layout/6` ·
behaviour preserving · flags: WASM · **open** · (determinism L8 / style-small L9)

Self-documented at solve_browser.pl:24–31 with a named production fix (exported
`browser.pl` module, plan §9) — recorded here so the promotion isn't lost, plus:
either export `arrange_best_layout/6` for benchmarks or comment subjects.pl as
white-box.

#### C26 — Per-node dynamic `search_seed/1` lookup on the gated hot path
core.pl:519–525 (`order_candidates/2` — one dynamic-predicate call per search node) ·
behaviour preserving (deterministic path provably unchanged: flag=false ≡ lookup-fails) ·
flags: perf-relevant · **open** · (census core, experiment shortlist rank 12 — the one
ranked candidate Phase 2 did not run)

**Fix:** hoist to a once-per-search flag threaded through the existing `select_inc`
state (extend `state/2` → `state/3` or pass a Perturb arg). MEDIUM plumbing; small
expected inference drop; measure like the X-steps.

#### C27 — findall/member used as list projection/filter where maplist/exclude is direct
arrange.pl:171–174 (`dropped_answers/3`), :320–324 (`unplaceable_words/2`), :619–620
(`fragment_answers/2`), :622–625 (`remaining_words/3`), :695–697/:350 (dropped-entry
projections), :870–877 (`arrange_candidates/6`) · behaviour preserving (order preserved
by all) · **open** · (style-arrange S3)

findall copies every collected term — and this very file documents term-copying as a
live footgun for non-ground `[Answer, _{...}]` entries (apply_move comment :992–999);
`remaining_words/3` is the one site whose copies flow onward into search state (safe
today only because removal keys on the ground answer). `maplist/3` / `exclude/3`
(`swi-manual/apply.md`) state the deterministic 1:1 maps directly and copy nothing.
Use named helpers, not lambdas (C13).

#### C28 — Greedy argmax sorts a whole list to take its stable maximum
arrange.pl `best_move` + `word_best_placement` (post-X3 shapes; baseline
:1016–1018/:1040–1046) · behaviour **risk if done wrong** · **open** · (style-arrange S4)

`sort(1, @>=, List, [Best|_])` is correct and documented (stability per
`swi-manual/builtinlist.md`) but O(n log n) per greedy step where a single-pass
`foldl` max keeping the incumbent on ties is O(n). **CAUTION:** the replacement compare
MUST be strict (`@>`) or later ties displace earlier ones and golden bytes change;
`max_member/2` is NOT safe (key ties would fall through to comparing payloads). Note X3
already changed these predicates' value shapes (`best_move/2` → `best_move/5`); re-derive
line numbers from the applied tree. Candidate lists are small — low value; verify with
`make golden`.

#### C29 — `pick_diverse_/6` recomputes `length(Acc)` every iteration
arrange.pl:843–853 (:845) · behaviour preserving · **open** · (style-arrange S5)

O(K) per pool element where a countdown argument is the standard accumulator shape.
Cold (once per `--candidates` run over a ≤20-element pool) — clarity more than cost.

#### C30 — `seed_candidates/2` sorts by negated keys where `sort(1, @>=)` is the file's own idiom
arrange.pl:951–961 (`neg_answer_len/2`) · behaviour preserving (both stable —
`swi-manual/builtinlist.md:93–97`) · **open** · (style-arrange S8)

The file leans on stable-descending `sort/4` everywhere else (:237, :353, :818, :1018,
:1046); sorting positive lengths deletes the negation helper. Both forms legitimately
idiomatic — hence low/possible.

#### C31 — `init_cell_vars/2`: element-wise `put_assoc` where `ord_list_to_assoc` is the direct constructor
fill.pl:75–77 · behaviour preserving · flags: perf-relevant (strict "counts change"
sense; grid-load, small N) · **open** · (style-fill F8)

Keys are already a strictly-ascending ordset (`mask_white_cells/3` ends in
`list_to_ord_set/2`). `pairs_keys_values(Pairs, WhiteCells, _FreshVars),
ord_list_to_assoc(Pairs, Assoc)` (`swi-manual/assoc.md:39`) replaces N rebalances.

#### C32 — `load_dict/3` word collection: findall+member+guard where `convlist/3` is the idiom
fill.pl:121 · behaviour preserving (word set and post-`sort/2` order identical) · flags:
**perf-relevant** (gated `load_inf` → re-baseline) · **open** · (style-fill F9)

`convlist(line_word, Lines, Ws0)` with `line_word(Line, Letters) :-
normalize_word(Line, Letters), Letters \== [].` (`swi-manual/apply.md:67`) — no findall
copy per word.

#### C33 — Artifact-read catch discards the root cause
fill.pl:641–643 · behaviour preserving (normal-case stderr unchanged) · **open** ·
(style-fill F11)

`catch(fill_read_artifact(...), _ReadErr, throw(error(fill_index_unreadable(File), _)))`
swallows whether the cause was a fastrw version mismatch, permissions, or truncation.
**Fix:** carry it — `throw(error(fill_index_unreadable(File), context(_, ReadErr)))`;
the cause then shows under verbose error printing.

#### C34 — `placed_bbox/4`: five passes over the cell set where one fold suffices
metrics.pl:73–82 (`min_list`/`max_list` ×4 after `pairs_keys_values/3`) · behaviour
preserving · flags: perf-relevant (arrange rescore path, not the fill ratchet) · **open**
· (style-fill F12)

A single `foldl/4` over `R-C` pairs computes all four bounds in one pass. Optional —
the current form is arguably more readable.

#### C35 — `slots_to_layout/3`: findall+member where maplist expresses the 1:1 map
fill.pl:393 · behaviour preserving (emit path, order preserved) · **open** ·
(style-fill F13)

`maplist(pw_answer_entry, Placed, InputWords)` with
`pw_answer_entry(PW, [A]) :- pw_answer(PW, A).`

#### C36 — `empty_slots/4` tests emptiness through `candidates/4` instead of the purpose-built counter
fill.pl:524–528 · behaviour preserving (cold failure-report path) · **open** ·
(style-fill F10)

Intent is "count is zero"; the module owns `candidate_count/4` for exactly this. Current
form is not expensive (the `= []` short-circuits) but reads as materialize-and-compare
and its cheapness depends on maplist evaluation order. **Fix:**
`candidate_count(Vars, DictByLen, Index, 0)`.

#### C37 — lint `tally/4` builds three throwaway lists; the file's own counting idiom is `aggregate_all(count, …)`
lint.pl:302–305 · behaviour preserving (values identical) · flags: contract-adjacent
(feeds the golden-locked JSON `summary` counts — rerun `make golden`) · **open** ·
(style-small L2)

`include(==('PASS'), Sevs, Ps), length(Ps, Pass)` ×3 → three
`aggregate_all(count, member(S, Sevs), N)` calls (`swi-manual/aggregate.md:32,66`;
lint.pl:277 already uses the idiom, library already imported).

#### C38 — Five copies of a read-terms-until-`clues/1` loop vs `read_file_to_terms/3`
core.pl:175–187 + verbatim clones in benchmarks/run_arrange.pl:137–143,
run_matrix.pl:126–133, start_sensitivity.pl:66–75, wasm/test/inference_parity.pl:73–80 ·
behaviour preserving for every existing fixture · **open** · (stdlib F5)

`read_file_to_terms(File, Terms, []), memberchk(clues(Words), Terms)`
(`swi-manual/readutil.md:54`; probe-verified against the bundled fixture). Honest
delta: the loop stops parsing at the first `clues/1`; `read_file_to_terms` parses the
whole file, so a fixture with a syntax error *after* the clues term would newly throw —
for bench/test fixtures arguably a feature, but it is a behaviour change. The real win
is deleting the four clones, not the core.pl copy. Combine with C19's steadfastness fix.

#### C39 — `read_file_lines/2` and `fill_file_sha256/2` hand-open files vs `read_file_to_string/3`
fill.pl:129–132 and :681–686 · behaviour preserving (probe: octet read byte-identical;
SHA-256 matches coreutils) · flags: **perf-relevant** (`read_file_lines` is inside the
gated `load_inf` — the swap moves it by a small constant → deliberate re-record) ·
**open** · (stdlib F6)

`read_file_to_string(File, Str, [])` / `read_file_to_string(File, Bytes,
[encoding(octet)])` (`swi-manual/readutil.md:48`). Pure-clarity win otherwise.

#### C40 — Bench harnesses parse JSON strings by hand vs `atom_json_dict/3`
benchmarks/check_baseline.pl:99–101,357–361 and check_fill_baseline.pl:109–111
(+ its parse_history_line twin) · behaviour preserving · **open** · (stdlib F9)

Same replacement as C9, three-plus call sites; harness-only code. The *file*-reading
`json_read_dict` sites are already the correct idiom — only the string-input sites are
reinventions.

### Nit

#### C41 — Adjacent-duplicate detection: `append(_, [D,D|_], Sorted)` → `nextto(D, D, Sorted)`
core.pl:909–918 (`check_unique_answers/1`), arrange.pl:637–643
(`check_fragment_unique/1`) · behaviour preserving (probe: same leftmost pair, same
order, same failure) · **open** · (stdlib F4)

`swi-manual/lists.md:56`. Readability ("two equal neighbours") plus avoiding the prefix
unification; runs once per invocation before the search — not perf-gated.

#### C42 — `bits_max_unch_run/2` hand-rolls run-length encoding vs `clumped/2`
metrics.pl:132–137 · behaviour preserving (probe-verified incl. all-checked case) ·
**open** · (stdlib F7 — honestly marked borderline)

`clumped/2` (`swi-manual/lists.md:157`) + `aggregate_all(max(N), member(0-N, RLE), …)`
is the precise library concept, but trades a transparent 5-line fold for RLE + a guarded
aggregate — reasonable people could call it churn. Recorded; adopt only if touching the
file anyway.

#### C43 — `row_string/2` is a misnamed one-line wrapper
stockgrid.pl:52,56 · behaviour preserving · **open** · (style-small L6)

The name says "string", the output is a char list, and `string_chars/2` is directly
maplist-able. Inline it or rename `row_chars/2`.

#### C44 — Cosmetic batch, core.pl
(a) naming drift: `Word` meaning entry vs answer across neighbouring predicates
(:309 vs :553–554), orphan `Letters2`/`Vals2`, `Length` vs `GridLen` in the geometry
helpers (:1099–1110); (b) whitespace debris: trailing blanks (:889,:892–893,:899,
:1078–1094), double blank lines inside predicate clause groups, 5-space indent at
:1095–1096 (the remove_x tabs were removed by X2.5's rewrite); (c) `==`/`\==` on
freshly-computed integers where `=:=`/`=\=` is the arithmetic idiom (:834,:837,:1089);
(d) `mrv_count_goal/7` builds a goal term where a named `viable_placement/7` helper
would be cross-referencer-visible (:402–407); (e) remaining legacy bare-`(`/`;` layout
sites not converted by X1 (:912–917, :959–964, :972–977 and the adj_is_free chains) ·
behaviour preserving · **open** · (style-core L2/L3/L6/L7 + M5 remainder)

#### C45 — Nit batch, arrange.pl
(a) three emit bodies repeat `current_output(Out), json_write_dict(Out, X), nl(Out)` —
one `emit_payload/1` helper (:383–395, :882–886); (b) the crossings weight `10000` is
commented but unnamed — `crossing_weight(10000).` matching the file's
`candidate_tau_pct/1` pattern (:1049–1054); (c) the two report-failure predicates are
near-duplicates differing only in contract-adjacent stderr wording — leaving them
verbatim is defensible (:301–316 vs :734–749) · behaviour preserving · **open** ·
(style-arrange N1/N2/N3)

#### C46 — Recorded-only stdlib notes (no action recommended)
(a) `dropped_answers/3` could use `subtract/3` — the one-pass findall is arguably as
clear (stdlib F8, borderline-churn); (b) `crypto_file_hash/3` would stream the SHA-256
but drags in the OpenSSL-backed ssl package — a real portability **regression for the
WASM deployment** the fill artifacts target; dictionaries are ~2 MB so slurping is a
non-issue — probably **not** worth doing (stdlib F10); (c) `dict_word_count/2` is
duplicated between fill.pl:674–677 and benchmarks/run_fill.pl:155–157 — the smell is the
duplication, not the fold; run_fill could call the module's copy (stdlib F11); (d)
`verbose_report/2` bypasses `print_message/2` — coherent for a CLI whose stderr contract
is prose lines; change only if log redirection is ever wanted (style-core L5) · **open**
(as a record) 

#### C47 — Cold RED-CONVERTIBLE cut conversions (optional robustness batch)
core.pl:116 (`require_strategy/1` → single ITE clause), :130 (`with_output/2` → guard
`File \== ''`), :164/167 (`load_clues_by_extension/3` → guarded catch-all), :335
(`select_word/9` mrv-seed → give cl.3 a `[_|_]` third-arg head; census probe `sw_*`
verified exact semantics after conversion; warm, non-default strategies only);
lint.pl:164 (`symmetry_override/4` → complement guard), :228 (`double_unch_end/2` →
`Bits \= [0,0|_]` guard) · behaviour preserving after conversion (bare removal is RED at
every one of these — that is the point of the class) · **open** · (census core C1–C4/C6,
census fill/small #5/#8/#15)

Each of these cuts is *correct today*; conversion buys guard-explicit clauses instead of
order-dependent commits. Zero perf case (all cold or off the gated path). Batch
opportunistically; goldens verify. (`reconcile_fragment_size` A4/A5 are handled by C18;
`meta_value_or` by C7; `grab_run`/`remove_slot` are KEEPs — §2.1.)

---

## Verified-deliberate hand-rolls — do NOT swap

Extends the P-audit's Appendix A. Each of these *looks* like a finding; each was checked
against the manual/library source in Phase 1 (and in several cases measured in Phase 2)
and is deliberate. Future audits should not re-flag them without new evidence.

- **`count_upto2/2`** (core.pl:424–434) — saturating solution counter via `nb_setarg/3`
  on a fresh local holder + `\+`-driven loop. This IS the manual's own idiom
  (`succeeds_n_times`, `swi-manual/manipterm.md:145`) plus an early exit; the comment
  names and rejects `aggregate_all/3` and `limit/2` with reasons. Hot, inference-gated,
  reentrant. Keep.
- **`alpha_chars/2`** (fill.pl:144–150) — first-order rewrite of
  `include(char_type_alpha)` with measured numbers in the comment (F-L1: load_inf
  26.6M → 10.7M). Reverting to `include/3` would fail the ratchet. Keep.
- **`select_min_count/3` + `min_count_walk/3` + `remove_slot/4`** (fill.pl:318–338) —
  semantically replaceable by `min_member/2` + `selectchk/3` (probe-verified: standard
  order on `cnt(C, slot(S,D,…))` decides on (C,S,D) and Start+Dir is unique), but this
  is the per-node hot path, byte-documented to reproduce `select_mrv`'s selection, and
  search_inf-gated; `remove_slot`'s cut is **negative-control-proven load-bearing**
  (§1c). A swap is pure churn with ratchet risk. Keep; the equivalence note is the value.
- **`remove_x/3` `==`-semantics** (core.pl, rewritten by X2.5) — the *rewrite* is
  applied, but the `==`-based, absent-tolerant contract stays deliberately NOT
  `selectchk/3`: arrange.pl:992–999 documents the exact footgun (findall-copied
  non-ground `[Answer, _{...}]` entries make `==`-removal drop nothing), and `remove_x`
  succeeds (identity) when the element is absent, where `selectchk` fails. Do not swap.
- **`pair_crossings/3` order** (core.pl:646–654) — the findall deliberately replays a
  former inline conjunction verbatim; "ORDER IS LOAD-BEARING" (golden-visible). Frozen.
- **`rows_of/3`** (core.pl:982–986) — no stdlib chunker exists in 10.1.10;
  `length/2`+`append/3` is the canonical idiom.
- **`enum_tokens`/`enum_clean`** (export.pl:67–90) — must distinguish *which* separator
  sits between runs; split-mode predicates don't apply. Domain logic. (The DCG rewrite
  was evaluated and declined — golden-locked output, marginal clarity.)
- **`split_runs`/`grab_run`** (stockgrid.pl:102–111) — maximal-run scanning over an
  ordset; no stdlib equivalent (per P-audit Appendix A too).
- **`same_text/2`, `signed/2`, `median_of/3`** (benchmarks) — no stdlib equivalents;
  the median deliberately unifies two prior definitions (F025).
- **`other_dir/2` in metrics vs core's unexported `swap_dir/2`** (metrics.pl:88–89) —
  arguably deliberate: metrics' header requires lint→metrics to stay off the solver
  substrate. If ever consolidated, the module-boundary rationale wins, not the DRY urge.
- **`\+ memberchk(Word, Used)`** (fill.pl:309) — at crossword scale (≤ ~40 placed
  words, early-fail atom-list compares) a list scan is the right cost/complexity trade;
  an ordset/assoc would complicate backtrack-restoration and move gated search_inf for
  no measured win. Recorded so a future reader knows it was considered (census E7 —
  deliberately not run).
- **`once(select(Best, Slots, Rest))`** (fill.pl:378) — the retained reference selector;
  its determinism is pinned by a plunit test (P13). Leave.
- **JSON emission** — everything goes through `json_write_dict/2,3`; no hand-rolled JSON
  anywhere. Word representation as lists of one-char atoms is the right choice for the
  cell-variable model (chars unify O(1) as interned atoms) — assessed, no finding.
- **The E2/E4 patch graveyard** (`X5_step2_REVERTED.patch`, `X5_step3_REVERTED.patch`)
  — measured losses (§1f). Do not resurrect without a workload that changes the
  premises (real GC pressure for E4; non-MRV bulk materialization for E2).

---

## What the codebase does well

Across all eight Phase-1 agents: **zero high-severity findings** — stated plainly, this
codebase came through a purity-focused, adversarially-probed audit with nothing worse
than measurement-integrity MEDs and consistency LOWs. The mechanical baseline was clean:
**zero compile warnings** across all modules and test suites, `check/0` **zero findings**
on every pass, all 17 `list_strings` hits intentional. Highlights the agents
independently called out:

- **Data-structure choices are excellent and argued in writing** — the arg/3-compound
  letter grid with trail-based undo, the `gs/2` boundary-grid bundle for O(1) legality,
  assoc-based MRV caches and cell maps, fill's assoc-of-ordsets index — each carrying
  the experiment reference that justified it. No list scans on large data anywhere.
- **The findall-copies-variables trap is dodged and documented** at every place it
  matters (fill's slot specs wired outside findall; `map_list_to_pairs` identity
  preservation; ground-answer-keyed removal) — the classic Prolog grid-model bug,
  handled with comments that teach the reader why.
- **Deliberate, measured first-order code on hot paths** with the library alternative
  named, rejected, and costed in comments (`alpha_chars/2`, `cell_var/3`,
  `count_upto2/2`) — and textbook `library(pairs)`/`sort/4`/`aggregate_all` pipelines
  everywhere else.
- **Determinism engineering throughout**: stable key-only sorts with the stability claim
  stated and correct per the manual; reified outcomes (`none`/`move(...)`,
  `ok/budget/exhausted`) dispatched on indexable heads; `call_with_inference_limit/3`
  handled exactly per `metacall.md` with why-no-catch comments (P2 preserved).
- **Error discipline**: shaped `error/2` formals with per-case DCG messages in every
  module, validation ordered for the clearest diagnosis, `setup_call_cleanup/3` on every
  stream, and a single CLI top-level catch giving curated messages + correct exit codes.
- **Comment discipline is exemplary** — the strongest self-documenting search code the
  style lane reported auditing: nearly every non-obvious choice carries its rationale
  and a cross-reference (spec §, experiment id, prior-audit P/F/R number), including
  honest debt records (solve_browser's spike header and owed state-reset invariant).
- **The dynamic-state surface is tiny and disciplined**: three single-writer config
  predicates with sentinel clears, two tabled pure functions, one local nb_setarg
  counter — and *zero* dynamics/globals in the six library modules.

---

## Methodology appendix

- **Benches.** Arrange: `benchmarks/run_arrange.pl` (5 core rungs ladder_09x09_08w …
  ladder_15x15_32w; 7 heavy rungs up to 38.5M inferences), ratcheted by
  `make bench-check` (`benchmarks/check_baseline.pl`, `search_inf` vs `baseline.json`,
  0.5% tolerance; the gated rungs run the strict `find_crossword(mrv_inc)` path only).
  Fill: `benchmarks/run_fill.pl` (7 core rungs sq04…g13), ratcheted by
  `make bench-fill-check` (`search_inf` + `load_inf`), plus the 11-rung
  `check_fill_identity.sh` byte-digest harness (including the heavy 15×15 ENABLE rung).
  The greedy arrange path is golden-covered but **not** inference-ratcheted — X3 built a
  purpose-made deterministic greedy bench for it.
- **Determinism discipline.** Every experiment ran its baseline **twice** before any
  edit and required byte-identical inference columns (all six did, every time). Greedy
  measurements used **one workload per fresh swipl process**, because `pair_crossings/3`
  is tabled and never abolished on the greedy path, making in-process counts
  call-history-dependent (the C1/C2 findings, discovered *by* this protocol); each
  workload was swept twice per state — counts identical both sweeps.
- **Wall-noise policy.** The host was shared between concurrent experiment agents (load
  0.7 → 4.3 during the campaign); same-code run-to-run wall spread reached **14%**
  (X2) and ±10–30% on dict-load-dominated commands (X6). Policy: wall numbers are quoted
  only where the effect dwarfs the spread (X3's −17…−18%) or as explicitly "indicative,
  interleaved A/B, direction-only" (X1); everything else defers to the deterministic
  inference metric, RSS (reproduced across interleaved pairs), or the pending
  idle-machine serial re-measure (§1d marker).
- **Worktree isolation.** Each experiment ran in its own `git worktree` off `69a3a2d`,
  one commit per step, patch files preserved; the E8 control was applied, measured, and
  `git revert`ed with a verified-empty `git diff` against base. (One X1 A/B script
  briefly checked out the user's branch by mistake; those runs were discarded and the
  A/B redone with explicit SHAs — checkouts only, no user-branch state modified.)
- **Negative-control design.** X6-E8 was chosen *because* the census predicted cut
  removal would lose there (semantics-safe but non-indexable dispatch) — a designed
  falsifiability check on the campaign's central thesis. It failed the ratchet exactly
  as predicted (§1c), calibrating the boundary rather than merely confirming the wins.
- **Probe corrections worth keeping.** `deterministic/1` inside an ITE condition sees
  the arrow's own else-CP (X4 discarded first measurements as artifacts); hand-rolled
  lowercase letters silently produce zero-count probe states because `answer_letters`
  canonicalizes case (X2); `swipl -q load.pl -g Goal` silently ignores `-g` — the
  correct form is `swipl -q -g Goal -t halt load.pl` (determinism agent; an early
  "clean" check/0 run was a no-op because of it).
- **Report inconsistencies found during synthesis** (flagged, not smoothed): the
  campaign's "17 APPLY steps" vs the 18 enumerated steps
  ([count note](#applied-step-count-note)); X6's summary line "7 once/1 wrappers
  deleted" vs the 9 its own steps enumerate (2+2+5 — the 9 is diff-verified); X1's prose
  "all six goldens" where every other report counts 8 (treated as a prose slip; X1's
  full `./run_tests.sh` gate passed); and the campaign summary's "zero `\+` in
  core/arrange" claim, corrected by grep in §2.3.

---

## Remediation tracker

Tick as landed. `→ §Cn` links to the finding; X-steps link to
§[3](#3-applied-remediation--the-phase-2-apply-steps).

**Applied (Phase-2 APPLY steps — committed 2026-07-06)**
- [x] **X1.1–X1.3** hot `(A;B),!` → ITE: check_letters / assign_letters / check_prev+next_cell — 2026-07-06
- [x] **X2.1–X2.5** indexing-carries-determinism batch: inc_counts, mrv_count, select_inc, assign_words(_inc), remove_x — 2026-07-06
- [x] **X3** greedy double-findall grid-snapshot elimination — 2026-07-06
- [x] **X4.1–X4.3** fill spine `[]`/`[_|_]`, shares_cell→ord_intersect, candidate_count arg-1 — 2026-07-06
- [x] **X5.1** bound_positions first-order kernel — 2026-07-06
- [x] **X6.B1–B5** once/cut hygiene + add_clue_nums det-at-source + find_crossword enumerated — 2026-07-06
- [ ] **post-commit**: `make bench-fill-record` (ratchet the 6-WIN fill baseline down)
- [x] serial idle-machine wall A/B for §1d (3 interleaved pairs; tables inline in §1d) — 2026-07-06

**Open — med**
- [ ] **C1** `med` — table abolition strategy-local; run-order-dependent inference counts → `core.pl:280`
- [ ] **C2** `med·WASM` — unbounded table growth on greedy path in persistent process → `core.pl:646,595`
- [ ] **C3** `med·WASM` — browser dynamic-state reset seam missing → `solve_browser.pl:71` · `core.pl:479` · `arrange.pl:101`
- [ ] **C4** `med` — no PlDoc/determinism annotations project-wide (or record the convention)
- [ ] **C5** `med` — P11 explicit-import follow-through; `load.pl` dies under autoload(false) → `load.pl:12` + module headers
- [ ] **C6** `med·WASM` — `library(http/json)` → `library(json)` → `solve_browser.pl:40` +siblings
- [ ] **C7** `med` — `meta_value/3`+`meta_value_or/4` → `library(option)` → `fill.pl:669`
- [ ] **C8** `med` — lint BFS → `library(ugraphs)` → `lint.pl:248`
- [ ] **C9** `med·WASM` — `solve_browser_str/2` → `atom_json_dict/3` → `solve_browser.pl:96`
- [ ] **C10** `med·WASM` — browser error_message hooks for the two unhooked formals → `solve_browser.pl:105`
- [ ] **C11** `med` — stockgrid_report double load → `stockgrid.pl:156`
- [ ] **C12** `med` — stockgrid quadratic report re-pairing → `stockgrid.pl:136`
- [ ] **C13** `med·perf` — interpreted yall pair → `arrange.pl:131` · `fill.pl:159`
- [ ] **C14** `med` — require_strategy message hard-codes the registry → `core.pl:236`
- [ ] **C15** `med·test` — `assign_word/10` dead arg 7 → `core.pl:665`
- [ ] **C16** `med·perf` — deprecated `delete/3` letter-normalization → `core.pl:596` · `metrics.pl:56`

**Open — low**
- [ ] **C17** `low` — silent-failure exports: lint_run profile / export_solve format / set_search_seed+set_verbose
- [ ] **C18** `low` — reconcile_fragment_size steadfastness (bogus throw) → `arrange.pl:544`
- [ ] **C19** `low` — load_clues .pl-path bound-output inconsistency → `core.pl:180`
- [ ] **C20** `low` — set_shuffle_seed side effects before failure → `core.pl:504`
- [ ] **C21** `low` — lint_solve emits before verdict binding → `lint.pl:311`
- [ ] **C22** `low` — unused exports `pw_end/2`, `valid_strategy/1` → `core.pl:47,56`
- [ ] **C23** `low·perf` — redundant once() post-X6.B4 → `fill.pl:392` · `stockgrid.pl:79`
- [ ] **C24** `low` — stale comments (10.0.2 pin, word-dicts, gs/2 header, fill export header)
- [ ] **C25** `low·WASM` — solve_browser module promotion (plan §9) + subjects.pl internals
- [ ] **C26** `low·perf` — per-node search_seed dynamic lookup hoist → `core.pl:519`
- [ ] **C27** `low` — arrange findall-as-map/filter → maplist/exclude (6 sites)
- [ ] **C28** `low·risk` — greedy argmax sort → strict-@> fold (golden-check)
- [ ] **C29** `low` — pick_diverse_ length-per-iteration → `arrange.pl:845`
- [ ] **C30** `low` — seed_candidates negated keys → `sort(1,@>=)` → `arrange.pl:951`
- [ ] **C31** `low·perf` — init_cell_vars → ord_list_to_assoc → `fill.pl:75`
- [ ] **C32** `low·perf` — load_dict findall → convlist → `fill.pl:121`
- [ ] **C33** `low` — artifact-read catch keeps the cause → `fill.pl:641`
- [ ] **C34** `low·perf` — placed_bbox single fold → `metrics.pl:73`
- [ ] **C35** `low` — slots_to_layout maplist → `fill.pl:393`
- [ ] **C36** `low` — empty_slots via candidate_count → `fill.pl:524`
- [ ] **C37** `low·contract-adjacent` — lint tally → aggregate_all → `lint.pl:302`
- [ ] **C38** `low` — read-clues loop ×5 → read_file_to_terms/3
- [ ] **C39** `low·perf` — read_file_lines/sha → read_file_to_string/3 → `fill.pl:129,681`
- [ ] **C40** `low` — bench harness string-JSON → atom_json_dict/3

**Open — nit**
- [ ] **C41** `nit` — duplicate detection via nextto/3 → `core.pl:909` · `arrange.pl:637`
- [ ] **C42** `nit` — bits_max_unch_run via clumped/2 (borderline churn) → `metrics.pl:132`
- [ ] **C43** `nit` — row_string misnamed wrapper → `stockgrid.pl:52`
- [ ] **C44** `nit` — core.pl cosmetic batch (naming / whitespace / `==`→`=:=` / mrv_count_goal / legacy layout)
- [ ] **C45** `nit` — arrange nit batch (emit_payload DRY / crossing_weight constant / report near-dupes)
- [ ] **C46** `nit` — recorded-only stdlib notes (subtract / crypto_file_hash / dict_word_count / verbose_report)
- [ ] **C47** `nit` — cold RED-convertible cut conversions (optional guard restructures)

## Remediation log

Append one entry per landed change (newest last): `YYYY-MM-DD · id · status · commit ·
one-line note`. Update the finding's status line and tick its box above at the same time.

- 2026-07-06 · audit opened · Phase 1 (8 agents) + Phase 2 (6 experiments, X1–X6) complete;
  C1–C47 raised (0 high · 16 med · 24 low · 7 nit), all `open`.
- 2026-07-06 · X1–X6 APPLY steps · **applied (uncommitted)** · — · all 18 APPLY-verdict
  steps applied to this worktree (staged; the overlapping X4.3/X5.1 candidate_count
  patches merged via git 3-way → 17 distinct patches; diff: arrange.pl/core.pl/fill.pl/
  lint.pl/tests/fill.plt, +223/−143). Gates re-run on the merged state: full
  `./run_tests.sh` ALL TESTS PASSED; `make bench-check` PASS (no change, 0 regressions);
  `make bench-fill-check` PASS — search_inf −0.97%/−2.39%/−2.47%/−0.34%/−4.60%/−3.42%/
  −2.82% on sq04/g11/g11_seed/sq05/g17/g21/g13 (6 WINs, 0 regressions).
  **Owed after commit:** `make bench-fill-record`; idle-machine serial wall A/B for §1d.
  Dropped by measurement (not applied): X5-E2 (+55.4% search_inf), X5-E4 (+0.21%, no GC
  to save), X6-E8 (negative control, reverted by design).
