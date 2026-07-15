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
| 5 | Fragment seeding: parse emit-schema fragment (canonical + thin `[{answer,row,col,dir}]` forms, desugared at the parse boundary), reconcile by answer, pre-place + validate, search remainder (words-only) | **done** (`seed_from_fragment/6` + `arrange_fragment_strict/6` / `arrange_fragment_best_effort/7` + `load_fragment/4`; +17 plunit, 2 golden checks over 1 golden file; thin seeds reach `fill --seeds` too — see the fill row) | AC-FRAG-1, AC-FRAG-2, AC-FRAG-3, AC-FRAG-4, AC-EMIT-2 |
| 6 | Candidates: distinctness from **constructor breadth + greedy diversity**, τ-filtered (not top-K B&B leaves) | **done** (`arrange_candidates/6`; greedy seed×corner pool, translation-invariant placement distance, τ=0.30; +7 plunit, +1 golden) | AC-ARR-7 |
| 7 | CLI + migration: subcommand dispatch; `--enumerate`; README/`run_tests.sh`/golden updates | **done** (`crosswordsmith` script; `crossword.pl` library-ized then removed with the old-style migration hint once no pre-migration consumers remained; goldens run through the CLI; README rewritten; +1 plunit) | AC-CLI-1, AC-CLI-2, AC-CLI-3, AC-ARR-6, AC-ARR-8 |
| — | ~~LNS polish pass~~ — **dropped** by the Phase-1.5 gate (no reward headroom; bound never pruned) | dropped | — |

Reachability calibration (`--check-target`, ε, τ) is a **required pre-weighting step** against `toc_demo` (cap inert) + `quality_22` (cap active), not an afterthought (OD-9).

---

## Shared substrate (spec §6) & CLI (spec §5)

| Component | Status | Notes / ACs |
|---|---|---|
| Input & `meta` passthrough (§6.1) | legacy | Exists in `prolog/crosswordsmith/core.pl`; validate against AC-IN-1/2. |
| Square grid model (§6.2) | legacy | `init_grid` etc.; reused as-is. |
| Clue numbering & enumeration (§6.3) | legacy | Numbering exists (AC-NUM-1); enumeration string derivation to confirm (AC-ENUM-1). |
| Metric predicates (§6.4) | legacy | Live in `prolog/crosswordsmith/metrics.pl` (the shared metric layer, formerly `quality.pl`); per spec §4 they stay a separate module, not lifted into core. The greedy constructor still departs to `arrange.pl` ([`source-structure-migration-plan.md`](./source-structure-migration-plan.md) Phase 3). |
| Emit / canonical JSON (§6.5) | legacy | Stable sorted-key JSON exists; confirm round-trip AC-EMIT-1/2. |
| Fragment-grid primitive (§6.6) | **done** (words-only; canonical + thin forms) | Realized in arrange Phase 5 (`arrange.pl`): emit-schema parse + reconcile + pin-via-legality-core + remainder search. AC-FRAG-1/2/3/4 + AC-EMIT-2 pass. The thin `[{answer,row,col,dir}]` form desugars at the parse boundary (`load_fragment/4`) for both consumers — `arrange --fragment` (framed by `--size`/`--max-size`, default 15) and `fill --seeds` (framed by the fill grid's own side). Letter-level pins (ninas) remain future scope. |
| CLI contract + migration (§5) | **done** (`arrange` verb) | `crosswordsmith` script: subcommand dispatch, bare→usage, old-style→unknown-verb usage error, `arrange` flags incl. `--enumerate`/`--candidates`/`--fragment`; the substrate is a library (`prolog/crosswordsmith/core.pl`); `--shuffle`/`--strategy` removed. AC-CLI-1/2/3, AC-ARR-6/8. The migration shim (root `crossword.pl`) and the dedicated old-style migration hint were retired once there were no pre-migration consumers. `lint`/`export`/`fill` verbs recognised but report not-built/deferred. |

---

## Flavour B (spec §8)

| Component | Spec | Status | Notes / ACs |
|---|---|---|---|
| `lint` (validator/profiles) | §8.1 | **done** (toc / blocked-uk / american / barred-ximenean) | `lint.pl` + `crosswordsmith lint` verb: per-word/per-rule PASS/WARN/FAIL report + verdict over the canonical layout JSON; reuses the shared metric predicates. AC-LINT-1/2/3/4. 18 plunit + 1 golden. **barred-ximenean built (OD-7 resolved, DP-3)**: primary-sourced Ximenean per-length unch band; symmetry relaxed to advisory. |
| `export` (ipuz v2 / Exolve) | §8.2 | **done** | `export.pl` + `crosswordsmith export --to ipuz\|exolve`: transformations of the canonical JSON. ipuz v2 (puzzle/solution/clues, enumerations derived from spaces/hyphens); Exolve plain text. Optional top-level `title`/`author` (json-output-spec §6.5) pass through when present — ipuz `title`/`author`, Exolve `exolve-title`/`exolve-setter`; nothing invented **except** Exolve's `exolve-title: Untitled` when a title is absent (Exet Save defence). The constant `exolve-id` is retired (Exolve auto-derives it) and the invented ipuz `"Untitled"` is dropped (P1 native-schema uplift, 2026-07-07). AC-EXP-1/3 (structure + enumeration/clue preservation); AC-EXP-2 structure. 24 plunit + 4 goldens. Real kotwords/Exet ingestion is a manual step. |
| Stock-grid library / profiles | §8.3 | **done** | `stockgrid.pl` + `grids/` (mask schema OD-5, grid set OD-6). Ships 3 lint-validated 180°-symmetric blocked grids (`blocked_13a`/`13b`/`15a`); legality is a CI regression (7 plunit). LOCKED. |
| `fill` engine (grid-first, open-dict) | §8.4 | **done** | `fill.pl` + `crosswordsmith fill --grid <mask> [--seeds <frag>] [--dict <words>]`: each white cell is a shared logical variable, MRV backtracking over an in-memory pattern index, fragment seeds as hard pins (canonical or thin §6.6 form, framed by the grid's own side); deterministic; reports unfillable slots + fails. 7 plunit + 1 golden at delivery (thin seeds later added +2 plunit + a seeded golden, 2 checks). AC-FILL-1…4. Bundles a small sample wordlist; real fills via `--dict UKACD18`. |
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
| OD-9 | `arrange` impl detail: ε/target/τ calibration, thin-form syntax, dup-answer disambiguation | **resolved (DP-1): 5:1 / ceil(L/2) / τ=0.30; thin-form deferred *(since landed 2026-07-06: AC-FRAG-4, arrange + fill)*; unique answers** |

---

## At a glance

- **Done — Flavour A `arrange` is feature-complete:** Phase 1 (oracle) + 1.5 (gate → DESCOPE) + **2 (strict)** + **3 (size framing)** + **4 (best-effort via greedy)** + **5 (fragment seeding)** + **6 (candidates)** + **7 (CLI + migration)**. Engine in `arrange.pl`; CLI in `crosswordsmith` (`crossword.pl` is now a library); **65 plunit tests (`tests/arrange.plt`) + 5 CLI golden checks (fixed + max + fragment canonical/thin + candidates)** wired into `run_tests.sh`/`make test` — full suite **306 plunit + 12 golden checks + 4 CLI exit-code checks + 2 fail-report checks, all green** (counts as of 2026-07-06, the thin-fragment-form change incl. its fill --seeds extension; the 2026-06-30 revamp audit + remediation resolved all 16 findings — see [`revamp-audit-findings.md`](./revamp-audit-findings.md)).
- **Flavour B `lint` (§8.1): done** — `lint.pl` + `crosswordsmith lint` verb (toc / blocked-uk / american / barred-ximenean). Consumes the canonical layout JSON, reports per-word/per-rule PASS/WARN/FAIL + verdict; 18 plunit + 1 golden. The barred-ximenean band is primary-sourced (DP-3).
- **Flavour B `export` (§8.2): done** — `export.pl` + `crosswordsmith export --to ipuz|exolve`: ipuz v2 JSON + Exolve plain text, enumerations derived from the answer; optional `title`/`author` pass-through (json-output-spec §6.5), the invented ipuz `"Untitled"` retired and the constant `exolve-id` dropped (P1 native-schema uplift, 2026-07-07); 24 plunit + 4 goldens. (Real kotwords/Exet round-trip is a manual verification step — checklist in [`exet-verification.md`](./exet-verification.md).)
- **Flavour B stock-grid library (§8.3): done** — DP-1 fixed OD-5 (mask schema); the build resolved OD-6 (ships 3 lint-validated grids). `stockgrid.pl` + `grids/`; 7 plunit.
- **Flavour B `fill` (§8.4): done** — DP-2 resolved OD-2/OD-4 (completing OD-1…4); `fill.pl` + `crosswordsmith fill`. Grid-first MRV backtracking over an in-memory pattern index, fragment seeds as hard pins, deterministic; 7 plunit + 1 golden. (Ships a sample wordlist; real fills via `--dict UKACD18`.)
- **Every spec'd component is now built, and OD-1…7 are all resolved.** The only thing still open is **OD-8** (backlog features in §8.5, each needs its own decision pass + spec section before implementation). The CLI does `arrange` / `lint` / `export` / `fill`.
- **Dropped (by the gate):** `arrange` B&B search loop, admissible bound, incremental delta, LNS polish.
- **Nothing in progress; nothing deferred or blocked.** Only OD-8 (unspec'd backlog) remains open.
- **Audited + remediated (2026-06-30): done.** A full multi-agent code review found 16 findings; **all 16 are resolved** and all four coverage gaps are closed. See *Audit & remediation* below and [`revamp-audit-findings.md`](./revamp-audit-findings.md).

### Research-derived backlog priorities (2026-07-15) — need spec + plan

A 2026 competitive re-scan ([`research/setter-tool-landscape-2026.md`](./research/setter-tool-landscape-2026.md), three verification passes) re-baselined crosswordsmith against the current field and identified the three highest-leverage **OD-8 backlog** items. Per the OD-8 discipline, each needs its own §10 decision pass + spec section + implementation plan before work starts. **Item 1 (scored fill) is now BUILT (DP-4 → §8.4a, amended by DP-5; FS-1 complete 2026-07-15);** the other two are recorded, not yet built, and still need their decision passes. Tracked in their home specs:

| # | Item | Home | Why now | Status |
|---|---|---|---|---|
| 1 | **Scored fill** (`fill --min-score` + fill-quality report) | design-spec **§8.4a** (DP-4 + DP-5) | The #1 competitive gap — every serious filler scores; crosswordsmith's scoreless MRV fills junk (measured, [`../benchmarks/fill_quality/`](../benchmarks/fill_quality/README.md)). License-clean via Spread the Wordlist. | **built (2026-07-15)** per [`plans/scored-fill-implementation.md`](./plans/scored-fill-implementation.md) — spec'd DP-4 → §8.4a (LOCKED), refined by DP-5 (native units, default prune `score ≥ 1`, uniform score 1); scored ingestion + score-desc ordering + `--min-score` prune + `--report-json` fill-quality report; FS-3(c) CI subset in `make test`/`make fuzz` (16 plunit + 3 goldens + 4 fuzz cases, no external deps); FS-3(a) STW/ingrid harness gate **PASSED** (native `--min-score 50` = mean 50.0 / min 50 / 0 below-clean on all four completable grids, `--report-json` agrees with post-hoc scoring); fill ratchet incl. `--heavy` PASS at +0.00% (search untouched). FS-3(b) mask spread done (2026-07-15, `matrix.sh` completion × min-score frontier over 8 masks — the ceiling bites only on blocked UK-style grids; `blocked_13a` at min ≤ 30 is the sole measured vs-ingrid gap); FS-4 **decided** (DP-6, 2026-07-15: `--budget N` + `fill --seed`/`--shuffle` spec'd as §8.4b [LOCKED] and **built the same day** per [`plans/fill-search-levers-implementation.md`](./plans/fill-search-levers-implementation.md) — budget option + ties-only load/slot seeds on the engine-internal PRNG, deterministic path byte-identical [ratchet +0.00%], seeded golden + 5 fuzz cases + shuffle-provenance roundtrip in CI; both measured as control/variety, NOT completion fixes — ×20 budget and 0/8 seeded reorders fail the gap row; crossing-aware forward-checking elevated to the named necessary path, its own future decision pass per the §8.5 row) |
| 2 | **`stats`/`inspect` + `diff` verbs** on `xword` | [`xword-status.md`](./xword-status.md) ★ · [`plans/xword-breadth-expansion.md`](./plans/xword-breadth-expansion.md) | The differentiated "breadth" first cut (read-only, deterministic, hub-amplifying) — research says compete on render + hub, not convert coverage. | needs decision pass + spec + plan |
| 3 | **Clue stockpile keyed by answer** (`meta`) | design-spec §8.5 ★ | Unserved across the entire surveyed field (both passes) — a genuine differentiator, low-effort, pure `meta` plumbing. | needs decision pass + spec + plan |

These do not change any `done` row above; they sit under the still-open **OD-8** (unspec'd backlog). The full ranked, surface-split roadmap is in the research doc §F.

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

### SWI-Prolog purity & cut audit (2026-07-06) — remediation done

A whole-core cut/impurity census plus six controlled, benchmark-gated experiments (the first audit in the series to *measure* its recommendations), grounded in the version-matched SWI 10.1.10 manual. Verdict: **zero high-severity findings** — the cleanest surface in the series. Applied 18 measured cut/purity steps (14 cuts + 9 `once/1` deleted; −6% RSS on 21×21 heavy rungs, −17…−18% greedy wall, −0.78…−2.9% fill `search_inf` where counted work shrank), with a published negative control proving where cuts must stay. Backlog: **63 findings** (C1–C47 sweep + C48–C63 browser.pl addendum), now **63/63 dispositioned** (fixed, with C42/C46 wont-fix and C44(c) reverted on measurement). One export cleanup is intentionally deferred — see the de-accretion roadmap. Full per-finding record + remediation log: [`prolog-purity-audit-findings.md`](./prolog-purity-audit-findings.md). **Status: done as of 2026-07-06.**

### De-accretion / retirement roadmap

The new `arrange` engine grew on top of the old machinery's primitives and orphaned its drivers. Tracking the cleanup so it doesn't just accrete:

- **Done (2026-07-03):** source-structure + module migration, per
  [`source-structure-migration-plan.md`](./source-structure-migration-plan.md)
  (all phases complete; per-phase record + deviations in the plan's
  checklists). End state: implementation under `prolog/crosswordsmith/`, one
  module per file (`crosswordsmith_{core,metrics,arrange,lint,export,
  stockgrid,fill}`) with explicit export lists; root `load.pl` (alias + the
  seven `use_module`s) is the single loader for driver/tests/benchmarks;
  the old root `crossword.pl` migration shim has been retired; `quality.pl`
  retired by rename to `metrics.pl`, the greedy constructor now lives in
  `arrange.pl`, and lint's metrics-only dependency boundary is enforced by
  imports. White-box tests reach internals as `Module:Pred(...)`; exports
  carry only the real inter-module/CLI/benchmark API. Notable deviations
  from the plan's verified map, found at module-ization: four extra exports
  with real consumers (`add_word_cells/3`, `emit_arrange/4`, `valid_loc/1`,
  `strategies/1` — closure references and latent paths that call-site greps
  missed; `list_undefined` is the gate that catches them).
- **Done:** removed the dead Phase-1.5 `gate_*` measurement harness + the orphaned `arrange_*_run` convenience runners (the `crosswordsmith` CLI is the entry point); `arrange.pl` 903 → 729 lines. The CLI fragment path now checks input uniqueness like the other modes.
- **Done (lint phase, opening move):** **deleted the dead `--quality` engine** from `quality.pl` (`quality_solve`/`quality_layout`/`grid_candidates`/`layout_score`/`quality_weights`/the floor subsystem) + its 9 tests; `quality.pl` 318 → 213 lines, re-framed as "shared metrics + the greedy density constructor." The lint-rule metrics (`word_meets_half`/`word_max_unch_run`/`checked_cells`/`dir_cells`/`word_checked_count`) now live in a file with no dead weight, ready for `lint` to consume.
- **Superseded (2026-07-02):** the deferred "relocate the shared metric predicates from `quality.pl` into `crossword.pl`" tidy-up is dropped — spec §4 now keeps metrics as a separate module (`prolog/crosswordsmith/metrics.pl`), preserving lint's metrics-only dependency boundary. Instead, `quality.pl` is renamed to `metrics.pl` and sheds the greedy constructor to `arrange.pl`: [`source-structure-migration-plan.md`](./source-structure-migration-plan.md) Phases 2–3.
- **Done (2026-07-07 — dissolution executed):** the legacy `crossword/3,4` top-level and the alternate strategies `baseline`/`mrv`/`mrv_capped` stay (the recorded keep decision stands — they are the evidence base for choosing `mrv_inc`, the only production strategy), now **relabelled and test-migrated off the production path**: `crossword/3,4` carries a RESEARCH/BENCHMARK-ONLY banner in `core.pl` (internal, routes through `find_crossword/6`'s C1/C48 seam); `strategies/1`'s PlDoc names all three live research-strategy consumers (`bench-matrix`, `start_sensitivity.pl`, and the `baseline` enumeration oracle guarding the **production** `all_crossword/5` / `arrange --enumerate` seam — do not "clean up" `all_crossword/5`, it is production); `tests/crossword.plt` renamed to `tests/core.plt` (suite-per-module naming) with its 8 production-machinery emit/CLI tests re-pointed at the exported seams (`solve_emit_json` helper = `crossword/4`'s exact composition) and 2 re-purposed tests kept under a research-surface banner (driver glue: default-strategy resolution + uniqueness-check-before-search ordering). Test-count delta: 0. The `legacy_main` migration shim was verified already absent from code (2026-07-07; only doc mentions remain). Ref: [`plans/legacy-surface-dissolution.md`](./plans/legacy-surface-dissolution.md).
- **Done (2026-07-06 — trigger fired):** export `crosswordsmith_arrange:arrange_best_layout/6`. The thin-form work (AC-FRAG-4) changed `arrange.pl`'s export surface, which was the recorded trigger. `/6` is now on the export list with a PlDoc contract (`is det`, probe-verified); the benchmark search sampler (`benchmarks/subjects.pl`) plain-imports it and the sanctioned white-box annotation (C25, commit `085a44c`) is retired. The access point itself is unchanged — every recorded `search_inf` baseline is defined against `/6` and the ratchet verified every gated count unmoved (+0.00%, core + `--heavy`). Ref: [`prolog-purity-audit-findings.md`](./prolog-purity-audit-findings.md) C25.
- **Done (2026-07-06):** export the five `crosswordsmith_fill` benchmark seams — `load_dict/3`, `fill_grid/4`, `fill_attempt/8`, `apply_seeds/4`, `seeded_slot/2`. Fresh observation during the C25 closure (no pre-existing deferred row): `benchmarks/fill_subjects.pl` reached these as un-annotated `crosswordsmith_fill:Pred(...)` white-box — the exact pattern C25 retired for arrange. Now on `fill.pl`'s export list with PlDoc contracts (det tags probe-verified on the benchmark call shapes, incl. both `fill_attempt/8` outcome paths); `fill_subjects.pl` plain-imports them. `fill_attempt/7` stays internal (default-budget delegate, mirroring `arrange_best_layout/5`). The seams themselves are unchanged — every `fill_baseline.json` / `fill_identity.sha256` count is defined against these predicates and the fill ratchet verified every gated count unmoved (+0.00%, core + `--heavy`; identity digests all match). The point-in-time campaign probes (`benchmarks/probe_f1/`, `probe_fh2/`) keep their historical qualified calls.
- **Done (2026-07-07 — invented-metadata retirement, P1 native-schema uplift):** `export.pl` stops fabricating data it never had. The ipuz transform no longer injects a `"title": "Untitled"` (both title/author are optional in ipuz, so kotwords is unaffected), and the Exolve transform drops the constant `exolve-id: crosswordsmith-export` (Exolve documents id as optional and auto-derives it from a content signature — a constant carried no information and did not preserve solving state across edits, per xword-spec ~:242). The one invented value that STAYS is Exolve's `exolve-title: Untitled` when the layout is title-less — a format-ecosystem requirement, not parity-chasing (Exet's Save crashes on a null title, [`exet-verification.md`](./exet-verification.md)), now emitted only when no title is present. In its place the engine passes through the optional top-level `title`/`author` when a layout carries them (json-output-spec §6.5; `arrange` still emits neither). Two goldens shifted intentionally (`export_bundled_17.ipuz` loses the title line; `.exolve` loses the id line); one new titled+authored fixture (`fixtures/titled_layout.json`) byte-pins pass-through. Export outside the gated paths — arrange + fill `search_inf` ratchets +0.00%. Ref: [`plans/native-schema-uplift-plan.md`](./plans/native-schema-uplift-plan.md) P1.
