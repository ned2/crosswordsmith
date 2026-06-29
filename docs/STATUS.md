# crosswordsmith — Build status

Progress tracker for the work specified in [`design-spec.md`](./design-spec.md). This is the **system of record** for what is done / in progress / blocked. The spec says *what* and the [`arrange` plan](./arrange-implementation-plan.md) says *how*; this file says *where we are*.

**Update discipline:** move a row's status in the **same commit/PR** as the work. A component is `done` only when its referenced `AC-*` criteria pass (see the spec). Don't add rows for unspecified scope — spec it first (spec §-change discipline).

**Status legend:**

| Mark | Meaning |
|---|---|
| `not started` | Specified, nothing built yet |
| `in progress` | Actively being built |
| `done` | Built; its ACs pass |
| `blocked` ⊘ | Waiting on an open decision (see Open decisions) |
| `deferred` | Spec'd as DEFERRED — not buildable until its decisions resolve |
| `legacy` | Exists in pre-spec code; works but not yet validated/refactored against the spec |

---

## Flavour A — `arrange` (spec §7 · plan phases)

**Architecture reframed by the design-exploration pass + confirmed by the Phase-1.5 gate (2026-06-30):** "deterministic MRV-first layout engine," not "branch-and-bound." The gate measured `best==first` on 4/5 fixtures, `prunes=0` on all 5, and the cap inert (0–1 of 8–26 words bind) — so the **B&B search superstructure, the admissible bound, the incremental delta, and the LNS pass are all dropped.** The engine = **construct + rescore (`layout_reward/4`) + emit**; best-effort + candidates ride the **greedy constructor**. See the plan's "Phase 1.5 — RESULT."

| Phase | Deliverable | Status | ACs |
|---|---|---|---|
| 1 | Scoring infra (`arrange.pl`): per-word `checked` + capped integer reward over a complete placement (the oracle) | **done** (oracle built + sanity-checked; `layout_reward/4`) | AC-ARR-9 |
| **1.5** | **Search-value gate**: measured whether search beats the first MRV incumbent | **done (2026-06-30): DESCOPE** — `best==first` 4/5, `prunes=0` 5/5, cap inert | — (decision gate) |
| 2 | Strict layout (fixed N): **construct + rescore + emit** — best of 4-corner MRV-inc placements, rescored; place-all-or-fail; budget-aware 3-outcome semantics (`placed`/`infeasible`-names-words/`not_proven`) | **done** (`arrange_strict_solve/3`; validated, deterministic) | AC-ARR-1, AC-ARR-3, AC-ARR-4, AC-ARR-9, AC-ARR-10 |
| 3 | Sizing + emit framing: `--size N`; `fixed` (exact N×N) vs `max` (tight square crop, side max(H,W)), default `max` | **done** (`emit_arrange/4`; both framings validated) | AC-ARR-5 |
| 4 | Best-effort (drop): served by the **greedy constructor path** (drops naturally), not a drop-branch on the strict DFS; lexicographic most-placed → reward across seeds; report dropped | **done** (`arrange_best_effort/6`; +3 plunit tests) | AC-ARR-2 |
| 5 | Fragment seeding: parse emit-schema fragment, reconcile by answer, pre-place + validate, search remainder (words-only v1) | **done** (`seed_from_fragment/6` + `arrange_fragment_strict/6` / `arrange_fragment_best_effort/7`; +10 plunit, +1 golden) | AC-FRAG-1, AC-FRAG-2, AC-FRAG-3, AC-EMIT-2 (AC-FRAG-4 thin-form deferred) |
| 6 | Candidates: distinctness from **constructor breadth + greedy diversity**, τ-filtered (not top-K B&B leaves) | not started | AC-ARR-7 |
| 7 | CLI + migration: subcommand dispatch; `--enumerate`; "did you mean `arrange`?" shim; README/`run_tests.sh`/golden updates | not started | AC-CLI-1…3, AC-ARR-6, AC-ARR-8 |
| — | ~~LNS polish pass~~ — **dropped** by the Phase-1.5 gate (no reward headroom; bound never pruned) | dropped | — |

Reachability calibration (`--check-target`, ε, τ) is a **required pre-weighting step** against `toc_demo` (cap inert) + `quality_22` (cap active), not an afterthought (OD-9).

---

## Shared substrate (spec §6) & CLI (spec §5)

| Component | Status | Notes / ACs |
|---|---|---|
| Input & `meta` passthrough (§6.1) | legacy | Exists in `crossword.pl`; validate against AC-IN-1/2. |
| Square grid model (§6.2) | legacy | `init_grid` etc.; reused as-is. |
| Clue numbering & enumeration (§6.3) | legacy | Numbering exists (AC-NUM-1); enumeration string derivation to confirm (AC-ENUM-1). |
| Metric predicates (§6.4) | legacy | Split across `crossword.pl`/`quality.pl`; lift `crossing_count`/`placed_bbox`/`word_meets_half` to shared layer (arrange Phase 1/7). |
| Emit / canonical JSON (§6.5) | legacy | Stable sorted-key JSON exists; confirm round-trip AC-EMIT-1/2. |
| Fragment-grid primitive (§6.6) | **done** (words-only v1) | Realized in arrange Phase 5 (`arrange.pl`): emit-schema parse + reconcile + pin-via-legality-core + remainder search. AC-FRAG-1/2/3 + AC-EMIT-2 pass; AC-FRAG-4 (thin form) deferred. |
| CLI contract + migration (§5) | legacy → not started | Old flag CLI exists; subcommand cutover is arrange Phase 7; AC-CLI-1…3. |

---

## Flavour B (spec §8)

| Component | Spec | Status | Notes / ACs |
|---|---|---|---|
| `lint` (validator/profiles) | §8.1 | not started | blocked-uk/toc/american profiles first; AC-LINT-1…4. Barred profile ⊘ on OD-7. |
| `export` (ipuz v2 / Exolve) | §8.2 | not started | Transformations of canonical JSON; AC-EXP-1…3. |
| Stock-grid library / profiles | §8.3 | ⊘ blocked | PARTIAL — schema + grid set open (OD-5, OD-6). |
| `fill` engine (grid-first, open-dict) | §8.4 | deferred | Not buildable until OD-1…4 resolved. |
| Backlog features | §8.5 | — | Unspecified; each needs its own decision pass before it gets a row here. |

---

## Open decisions (spec §10) — these gate `blocked`/`deferred` rows above

| # | Gates | Status |
|---|---|---|
| OD-1 | `fill`: blocked-only v1 vs barred-compatible from start | open |
| OD-2 | `fill`: dictionary integration + default lexicon | open |
| OD-3 | `fill`: seed/fragment semantics into open-dict fill | open |
| OD-4 | `fill`: v1 profiles + no-fill failure contract | open |
| OD-5 | Stock-grid: template schema | open |
| OD-6 | Stock-grid: which grids seed the library | open |
| OD-7 | `lint`: barred-Ximenean unch table + barred symmetry codes | open |
| OD-8 | Backlog: per-feature decision pass | open |
| OD-9 | `arrange` impl detail: ε/target/τ calibration, thin-form syntax, dup-answer disambiguation | tracked in plan |

---

## At a glance

- **Done:** `arrange` Phase 1 (oracle) + 1.5 (gate → DESCOPE) + **2 (strict)** + **3 (size framing)** + **4 (best-effort via greedy)** + **5 (fragment seeding)**. All in `arrange.pl`; **30 plunit tests (`tests/arrange.plt`) + golden regression (fixed + max + fragment) wired into `run_tests.sh`/`make test`** — full suite 111/111 + 4 goldens green.
- **Next buildable, unblocked:** `arrange` **Phase 6** (candidates) → **7** (CLI + migration).
- **Blocked:** stock-grid library (OD-5/6), `lint` barred profile (OD-7).
- **Deferred:** `fill` engine (OD-1…4).
- **Dropped (by the gate):** `arrange` B&B search loop, admissible bound, incremental delta, LNS polish.
- **Nothing in progress.**
