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
| 6 | Candidates: distinctness from **constructor breadth + greedy diversity**, τ-filtered (not top-K B&B leaves) | **done** (`arrange_candidates/6`; greedy seed×corner pool, translation-invariant placement distance, τ=0.30; +7 plunit, +1 golden) | AC-ARR-7 |
| 7 | CLI + migration: subcommand dispatch; `--enumerate`; "did you mean `arrange`?" shim; README/`run_tests.sh`/golden updates | **done** (`crosswordsmith` script; `crossword.pl` library-ized; goldens run through the CLI; README rewritten; +1 plunit) | AC-CLI-1, AC-CLI-2, AC-CLI-3, AC-ARR-6, AC-ARR-8 |
| — | ~~LNS polish pass~~ — **dropped** by the Phase-1.5 gate (no reward headroom; bound never pruned) | dropped | — |

Reachability calibration (`--check-target`, ε, τ) is a **required pre-weighting step** against `toc_demo` (cap inert) + `quality_22` (cap active), not an afterthought (OD-9).

---

## Shared substrate (spec §6) & CLI (spec §5)

| Component | Status | Notes / ACs |
|---|---|---|
| Input & `meta` passthrough (§6.1) | legacy | Exists in `crossword.pl`; validate against AC-IN-1/2. |
| Square grid model (§6.2) | legacy | `init_grid` etc.; reused as-is. |
| Clue numbering & enumeration (§6.3) | legacy | Numbering exists (AC-NUM-1); enumeration string derivation to confirm (AC-ENUM-1). |
| Metric predicates (§6.4) | legacy | Live in `quality.pl` (shared metric layer); per spec §4 they stay a separate module — target `prolog/crosswordsmith/metrics.pl` via [`source-structure-migration-plan.md`](./source-structure-migration-plan.md) Phases 2–3, not lifted into `crossword.pl`. |
| Emit / canonical JSON (§6.5) | legacy | Stable sorted-key JSON exists; confirm round-trip AC-EMIT-1/2. |
| Fragment-grid primitive (§6.6) | **done** (words-only v1) | Realized in arrange Phase 5 (`arrange.pl`): emit-schema parse + reconcile + pin-via-legality-core + remainder search. AC-FRAG-1/2/3 + AC-EMIT-2 pass; AC-FRAG-4 (thin form) deferred. |
| CLI contract + migration (§5) | **done** (`arrange` verb) | `crosswordsmith` script: subcommand dispatch, bare→usage, old-style→migration hint, `arrange` flags incl. `--enumerate`/`--candidates`/`--fragment`; the substrate is a library (`prolog/crosswordsmith/core.pl`) and root `crossword.pl` a message-only migration shim; `--shuffle`/`--strategy` removed. AC-CLI-1/2/3, AC-ARR-6/8. `lint`/`export`/`fill` verbs recognised but report not-built/deferred. |

---

## Flavour B (spec §8)

| Component | Spec | Status | Notes / ACs |
|---|---|---|---|
| `lint` (validator/profiles) | §8.1 | **done** (toc / blocked-uk / american / barred-ximenean) | `lint.pl` + `crosswordsmith lint` verb: per-word/per-rule PASS/WARN/FAIL report + verdict over the canonical layout JSON; reuses the shared metric predicates. AC-LINT-1/2/3/4. 18 plunit + 1 golden. **barred-ximenean built (OD-7 resolved, DP-3)**: primary-sourced Ximenean per-length unch band; symmetry relaxed to advisory. |
| `export` (ipuz v2 / Exolve) | §8.2 | **done** | `export.pl` + `crosswordsmith export --to ipuz\|exolve`: transformations of the canonical JSON. ipuz v2 (puzzle/solution/clues, enumerations derived from spaces/hyphens); Exolve plain text. AC-EXP-1/3 (structure + enumeration/clue preservation); AC-EXP-2 structure. 11 plunit + 2 goldens. Real kotwords/Exet ingestion is a manual step. |
| Stock-grid library / profiles | §8.3 | **done** | `stockgrid.pl` + `grids/` (mask schema OD-5, grid set OD-6). Ships 3 lint-validated 180°-symmetric blocked grids (`blocked_13a`/`13b`/`15a`); legality is a CI regression (7 plunit). LOCKED. |
| `fill` engine (grid-first, open-dict) | §8.4 | **done** | `fill.pl` + `crosswordsmith fill --grid <mask> [--seeds <frag>] [--dict <words>]`: each white cell is a shared logical variable, MRV backtracking over an in-memory pattern index, fragment seeds as hard pins; deterministic; reports unfillable slots + fails. 7 plunit + 1 golden. AC-FILL-1…4. Bundles a small sample wordlist; real fills via `--dict UKACD18`. |
| Backlog features | §8.5 | — | Unspecified; each needs its own decision pass before it gets a row here. |

---

## Open decisions (spec §10) — these gate `blocked`/`deferred` rows above

| # | Gates | Status |
|---|---|---|
| OD-1 | `fill`: blocked-only v1 vs barred-compatible from start | **resolved (DP-1): blocked-only** |
| OD-2 | `fill`: dictionary integration + default lexicon | **resolved (DP-2): in-memory pattern index; UKACD18 default via --dict** |
| OD-3 | `fill`: seed/fragment semantics into open-dict fill | **resolved (DP-1): pin-and-fill (fragment primitive)** |
| OD-4 | `fill`: v1 profiles + no-fill failure contract | **resolved (DP-2): stock-grids-as-profiles; report unfillable slots + fail** |
| OD-5 | Stock-grid: template schema | **resolved (DP-1): black-square mask (slots derived)** |
| OD-6 | Stock-grid: which grids seed the library | **resolved: ships blocked_13a/13b/15a (lint-validated)** |
| OD-7 | `lint`: barred-Ximenean unch table + barred symmetry codes | **resolved (DP-3): Ximenean band primary-sourced + built; symmetry relaxed** |
| OD-8 | Backlog: per-feature decision pass | open |
| OD-9 | `arrange` impl detail: ε/target/τ calibration, thin-form syntax, dup-answer disambiguation | **resolved (DP-1): 5:1 / ceil(L/2) / τ=0.30; thin-form deferred; unique answers** |

---

## At a glance

- **Done — Flavour A `arrange` is feature-complete:** Phase 1 (oracle) + 1.5 (gate → DESCOPE) + **2 (strict)** + **3 (size framing)** + **4 (best-effort via greedy)** + **5 (fragment seeding)** + **6 (candidates)** + **7 (CLI + migration)**. Engine in `arrange.pl`; CLI in `crosswordsmith` (`crossword.pl` is now a library); **38 plunit tests (`tests/arrange.plt`) + 4 CLI goldens (fixed + max + fragment + candidates)** wired into `run_tests.sh`/`make test` — full suite **168 plunit + 8 goldens + 3 CLI exit-code checks, all green** (counts current as of the 2026-06-30 revamp audit + remediation, all 16 findings resolved — see [`revamp-audit-findings.md`](./revamp-audit-findings.md)).
- **Flavour B `lint` (§8.1): done** — `lint.pl` + `crosswordsmith lint` verb (toc / blocked-uk / american / barred-ximenean). Consumes the canonical layout JSON, reports per-word/per-rule PASS/WARN/FAIL + verdict; 18 plunit + 1 golden. The barred-ximenean band is primary-sourced (DP-3).
- **Flavour B `export` (§8.2): done** — `export.pl` + `crosswordsmith export --to ipuz|exolve`: ipuz v2 JSON + Exolve plain text, enumerations derived from the answer; 11 plunit + 2 goldens. (Real kotwords/Exet round-trip is a manual verification step — checklist in [`exet-verification.md`](./exet-verification.md).)
- **Flavour B stock-grid library (§8.3): done** — DP-1 fixed OD-5 (mask schema); the build resolved OD-6 (ships 3 lint-validated grids). `stockgrid.pl` + `grids/`; 7 plunit.
- **Flavour B `fill` (§8.4): done** — DP-2 resolved OD-2/OD-4 (completing OD-1…4); `fill.pl` + `crosswordsmith fill`. Grid-first MRV backtracking over an in-memory pattern index, fragment seeds as hard pins, deterministic; 7 plunit + 1 golden. (Ships a sample wordlist; real fills via `--dict UKACD18`.)
- **Every spec'd component is now built, and OD-1…7 are all resolved.** The only thing still open is **OD-8** (backlog features in §8.5, each needs its own decision pass + spec section before implementation). The CLI does `arrange` / `lint` / `export` / `fill`.
- **Dropped (by the gate):** `arrange` B&B search loop, admissible bound, incremental delta, LNS polish.
- **Nothing in progress; nothing deferred or blocked.** Only OD-8 (unspec'd backlog) remains open.
- **Audited + remediated (2026-06-30): done.** A full multi-agent code review found 16 findings; **all 16 are resolved** and all four coverage gaps are closed. See *Audit & remediation* below and [`revamp-audit-findings.md`](./revamp-audit-findings.md).

### Audit & remediation (2026-06-30) — done

A full 7-lane multi-agent code review (spec-conformance AC-by-AC, per-module correctness for `arrange`/`fill`/`lint`/`stockgrid`/`export`, CLI + cross-cutting invariants, tests + docs integrity) with adversarial per-finding verification produced **16 findings** (2 high, 2 med, 11 low, 1 nit); **all 16 are fixed**. Highlights:

- **The two high defects were real correctness bugs, now fixed + regression-tested:** hyphenated answers placed a literal `-` grid cell (R1); a `fill` seed absent from the dictionary rejected a grid that has a legal fill (R2).
- **`--check-target` (§7.2 MUST) was implemented** during remediation (R8), resolving a LOCKED spec self-contradiction.
- **Doc/spec drift swept:** stale test counts, `AGENTS.md`, §6.5 word-object fields, the stock-grid `symmetry` annotation, and the mislabelled UKACD18 license.

The audit's **four coverage gaps are all closed:**

- **INV-4 license/provenance** — audited every bundled/vendored asset; **no AC-X-4 violation**. Corrected UKACD18's license (redistributable freeware, ship its notice verbatim — **not** BSD-3) and flagged the vendored SWI manual as CC BY-SA 3.0.
- **§7.3 worst-case latency** — the strict 4-corner sweep now shares one inference budget (R7; `toc_demo`@15 ≈100 s → ≈28 s).
- **INV-2 determinism** — `tests/determinism_fuzz.sh` (`make fuzz`): a 54-case verb × flag × degenerate-input fuzz, each run as 3 processes for byte-identity; **INV-2 holds (0 nondeterministic cases, 0 hangs)**.
- **AC-EXP-2 Exet round-trip** — un-automatable in-repo; a step-by-step manual checklist + audit log ships at [`exet-verification.md`](./exet-verification.md) (the one remaining human-in-the-loop step).

Full per-finding record + remediation log: [`revamp-audit-findings.md`](./revamp-audit-findings.md). Post-remediation suite: **168 plunit + 8 goldens + 3 CLI exit-code checks** (`make test`), plus the on-demand `make fuzz`.

### SWI-Prolog idiom audit (2026-07-01) — remediation done

A 5-lane parallel review swept the *whole* current core (`crossword.pl`, `arrange.pl`, `fill.pl`, `lint.pl`, `quality.pl`, `export.pl`, `stockgrid.pl` + harness) specifically for **predicate-use correctness, stdlib reuse, and idiom**, grounding every claim in the version-matched SWI manual under [`reference/swi-manual/`](./reference/swi-manual/). Verdict: high-quality, idiomatic Prolog — **no deprecated predicates, no `format/2` mismatches, no state leaks**. Produced **17 findings (0 high · 4 med · 7 low · 5 nit)**; **all 17 are now fixed** under the adopted SWI 10.1.10 runtime:

- **P1 (med, the only behaviour-risk item): fixed** — a one-cell `is_end_cell(down,…)` off-by-one (`crossword.pl`, `>=` → `>`) let a down word merge collinearly at cell `(L-1)*L`; fixed + 4 regression plunit, all goldens byte-identical.
- **P2 (med): fixed** — dropped the broad `catch/3` in both engines (`call_with_inference_limit/3` handles the budget itself, so it only ever swallowed real errors); a genuine error now surfaces via `main/0` instead of as "infeasible". +2 plunit.
- **P3/P4 (med): fixed** — two hot-path efficiency wins, benchmarked in INFERENCES: `fill` MRV counting no longer materializes candidate lists (**−56%** on the counting map), and the checked-bitmap metric is hoisted into `quality.pl` with `dir_cells` computed once per lint run (**−31%** `lint_run(toc)`).
- **P8 (low·B): fixed after SWI 10.1.10 adoption** — the existing `entry_letters/2` now replaces the three duplicated inline letter-normalization blocks. This cleanup was originally rejected on SWI 10.0.2 because the `assign_words_inc/9` production search loop triggered a GC-related stack overflow; SWI 10.1.10 no longer reproduces the failure, and the full suite passes after regenerating JSON goldens for 10.1.10 pretty-printer whitespace. **P17** was resolved as **doc** (invariant comments, not throws).

Post-remediation suite: **179 plunit + 8 goldens + 3 CLI exit-code checks** (`make test`), plus `make fuzz` (54 cases). Per-finding record + remediation log (checklist + per-commit hashes): [`prolog-idiom-audit-findings.md`](./prolog-idiom-audit-findings.md). **Status: done as of 2026-07-02.**

### De-accretion / retirement roadmap

The new `arrange` engine grew on top of the old machinery's primitives and orphaned its drivers. Tracking the cleanup so it doesn't just accrete:

- **In progress:** migrate source structure toward an SWI-Prolog `prolog/`
  library layout, then introduce modules and explicit exports from leaves
  inward. Tracked in
  [`source-structure-migration-plan.md`](./source-structure-migration-plan.md).
  Landed 2026-07-02: Phase 0 (spec §4 alignment) and Phase 1 (implementation
  files moved to `prolog/crosswordsmith/`, `crossword.pl` split into
  `core.pl` + a root message-only shim that loads nothing, root `load.pl` as
  the single owner of load order — driver/tests/benchmarks all load through
  it).
- **Done:** removed the dead Phase-1.5 `gate_*` measurement harness + the orphaned `arrange_*_run` convenience runners (the `crosswordsmith` CLI is the entry point); `arrange.pl` 903 → 729 lines. The CLI fragment path now checks input uniqueness like the other modes.
- **Done (lint phase, opening move):** **deleted the dead `--quality` engine** from `quality.pl` (`quality_solve`/`quality_layout`/`grid_candidates`/`layout_score`/`quality_weights`/the floor subsystem) + its 9 tests; `quality.pl` 318 → 213 lines, re-framed as "shared metrics + the greedy density constructor." The lint-rule metrics (`word_meets_half`/`word_max_unch_run`/`checked_cells`/`dir_cells`/`word_checked_count`) now live in a file with no dead weight, ready for `lint` to consume.
- **Superseded (2026-07-02):** the deferred "relocate the shared metric predicates from `quality.pl` into `crossword.pl`" tidy-up is dropped — spec §4 now keeps metrics as a separate module (`prolog/crosswordsmith/metrics.pl`), preserving lint's metrics-only dependency boundary. Instead, `quality.pl` is renamed to `metrics.pl` and sheds the greedy constructor to `arrange.pl`: [`source-structure-migration-plan.md`](./source-structure-migration-plan.md) Phases 2–3.
- **Decided — keep:** the legacy `crossword/3,4` top-level and the alternate strategies `baseline`/`mrv`/`mrv_capped` stay as a **benchmark-only research surface** (the evidence base for choosing `mrv_inc`, the only production strategy); to be relabelled/test-migrated off the production path during the dissolution above. The `legacy_main` migration shim retires once old muscle-memory fades.
