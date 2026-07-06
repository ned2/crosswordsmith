# Legacy-surface dissolution — relabel + test-migrate the research surface

Status: REVIEWED 2026-07-07 — adversarially reviewed against the repo; all
inventory claims confirmed, corrections folded in (marked ⊳). Closes the last
open line of `docs/STATUS.md`'s de-accretion roadmap ("Decided — keep: the
legacy `crossword/3,4` top-level and the alternate strategies ... to be
relabelled/test-migrated off the production path ... The `legacy_main`
migration shim retires once old muscle-memory fades.", STATUS.md:149).

Branch: `deaccrete/legacy-surface-dissolution` from current main. This plan doc
lands at `docs/plans/legacy-surface-dissolution.md` as the first commit.

## 1. Verified inventory (2026-07-07) — the roadmap sentence has drifted

| Item | Verified reality |
|---|---|
| `crossword/3,4` (core.pl:336–344) | Already INTERNAL — not on the export list; core.pl:10–12 header already labels it "legacy ... benchmark-only research surface". Sole callers anywhere: `tests/crossword.plt` — 10 qualified sites (:202, :213, :224, :234, :380, :431, :437, :488, :490, :503), ⊳ all of them `crossword/3`; nothing exercises `/4`'s explicit-strategy arm. |
| `all_crossword/5` | EXPORTED with a PRODUCTION consumer: `arrange.pl:1188–1193`, the `--enumerate` seam (design-spec §7.1, AC-ARR-8 at design-spec.md:267; guarded by tests/arrange.plt:518). **Not legacy** — must NOT be relabelled research. |
| `find_crossword/6` | Exported. The shared driver: `crossword/3,4` and `all_crossword/5` route through it; `benchmarks/run_matrix.pl:114` and tests call it directly. Carries the C1/C48 seam invariant (comment core.pl:369–415). ⊳ `find_crossword/5` is INTERNAL — re-points must not quietly substitute it. |
| Strategies `baseline`/`mrv`/`mrv_capped` | Research-only as *production* strategies, but ⊳ live consumers are THREE, not one: (1) `make bench-matrix` → `run_matrix.pl` (strategies/1's only out-of-module caller, :56); (2) `benchmarks/start_sensitivity.pl:40,66` sweeps `[baseline, mrv_inc]` (run per docs/experiments.md:58; on no make target); (3) `tests/crossword.plt:295–306` uses `baseline` as the enumeration-count **correctness oracle** for `all_crossword/5` (role recorded at docs/benchmark-rework-plan.md:274) — that is a guard of the production `--enumerate` seam, NOT a research-surface test. |
| `default_strategy/1` (= `mrv_inc`) | Production: imported and called by `arrange.pl` (:91, :1189). |
| `legacy_main` | ALREADY RETIRED — zero references in code; only docs mention it (STATUS.md:149, source-structure-migration-plan.md:60,261). |
| `tests/crossword.plt` (599 lines, 8 sub-suites at :31,:91,:133,:174,:316,:351,:458,:517) | **core.pl's suite under its pre-migration name**. Most sub-suites are legitimate white-box guards of the shared substrate; the :517 `quality` block covers metrics.pl. Legacy-driver usage is confined to the 10 `crossword/3` sites in `solver`/`json_input`/`cli`. |
| ⊳ Suite discovery | `run_tests.sh` just runs `swipl -q tests/run_tests.pl`; the explicit consult list at **`tests/run_tests.pl:17`** (`consult('tests/crossword.plt')`) is the single point of failure for the rename. No CI config exists; Makefile and the wasm harness reference no .plt filename. |
| ⊳ Other filename/phrasing references | `core.pl:289–290` comment and `docs/json-output-spec.md:306,:354` name `crossword.plt` (the :354 one cites a line number). `fill.pl:84`, `fill.pl:690`, `core.pl:1340` describe `crossword/4` as "the solve path" — misleading once it is bannered research-only. Historical audit/plan docs stay point-in-time (untouched). |

## 2. Work items (one commit after the plan commit; tracker rows move with it)

### W1 — relabel in place (core.pl + comment sweep)
- No code motion. Give `crossword/3,4` a plain-prose banner (⊳ AGENTS.md: `%!`
  PlDoc + verified tags are the convention for exported predicates; internals
  use `%` prose — so prose, and if any determinism note is written it must say
  the search is `nondet` through `find_crossword/6`, core.pl:427–434, not det).
  The banner states: RESEARCH/BENCHMARK-ONLY; sanctioned consumers are this
  repo's tests; ⊳ it ROUTES THROUGH the production driver `find_crossword/6`
  and its C1/C48 reset seam (do not describe it as isolated from production
  machinery — the seam-invariant comment at :369–371 names it as a router and
  must stay true).
- Sharpen `strategies/1`'s PlDoc: the registry exists for the research matrix;
  production resolves through `default_strategy/1` only. ⊳ Name all three
  live strategy consumers (bench-matrix, start_sensitivity, the baseline
  enumeration oracle) — do NOT claim bench-matrix is the only reach.
- ⊳ Comment-consistency sweep in the same commit: rephrase `fill.pl:84`,
  `fill.pl:690`, `core.pl:1340` so they reference the production pipeline
  (e.g. "the legacy crossword/4 driver" or the underlying seam) rather than
  implying crossword/4 is the production solve path.
- Explicitly REJECTED alternative: extracting the legacy surface to a
  `research.pl` module — the strategy clauses interleave with the shared
  driver and the C1/C48 memo-hygiene seam; extraction churns a seam every
  gated baseline depends on, for zero behavioral gain.

### W2 — test triage (tests/crossword.plt → tests/core.plt)
- **Rename `tests/crossword.plt` → `tests/core.plt`** (aligns with the
  one-suite-per-module convention, AGENTS.md:42; the other six suites all
  match module names — ⊳ noting core.plt still hosts the metrics.pl `quality`
  block, so "one per module" stays approximate). ⊳ Update the consult at
  `tests/run_tests.pl:17`. Hard acceptance: the plunit pass/test count after
  the rename is IDENTICAL to before (record both numbers) — a silently
  un-run suite is this change's biggest real risk.
- Re-point the legacy call sites at the exported seams. `crossword/4`'s body
  (core.pl:340–344) is exactly `check_unique_answers/1 → find_crossword/6 →
  assign_clue_numbers/2 → emit_json/3`, all exported; ⊳ the sites are
  `crossword/3`, so re-points also call `default_strategy/1` for the Strategy
  arg (use exported `find_crossword/6`; `find_crossword/5` is internal).
  - `solver` (:202, :213, :224): JSON emit shape → the composition directly.
  - ⊳ `solver` (:234, `duplicate_answer_rejected`): guards DRIVER GLUE — that
    the uniqueness check runs BEFORE the search. Re-pointing it to
    `check_unique_answers/1` alone would silently un-guard that ordering:
    this is one of the KEPT bannered research-surface tests (below).
  - `json_input` (:380, :431, :437): loading round-trips → the composition.
  - `cli` (:488, :490, :503): `with_output/2` contract — goal-agnostic
    meta-predicate (core.pl:183–189), safe to re-point.
  - ⊳ `count_upto2`/oracle tests (:295–306): NOT legacy — they guard the
    production `--enumerate` seam using `baseline` as the oracle. Leave them
    where they are; do not banner them as research.
- KEEP a small bannered research-surface section: 1–2 ⊳ RE-PURPOSED EXISTING
  tests (not additions — keeping the count-delta promise exact) that still
  call `crosswordsmith_core:crossword/3`: :234 (driver glue/ordering) and one
  full-emit test (e.g. :202's shape), under a banner saying they guard the
  research surface + default-strategy resolution so `bench-matrix`'s entry
  doesn't rot untested between campaigns.
- No test deletions without naming the equivalent coverage in the commit
  message. Expected count delta: exactly 0.

### W3 — docs (same commit as W1+W2)
- STATUS.md:149 roadmap row → **Done (dated)**, rewritten to verified
  reality: crossword/3,4 internal + bannered; strategies relabelled with all
  three consumers named; `all_crossword/5` is production (`--enumerate`) —
  state explicitly so a future sweep doesn't "clean it up"; `legacy_main`
  verified absent 2026-07-07 (retired earlier, undated).
- ⊳ Filename sweep for the rename: `core.pl:289–290` comment,
  `docs/json-output-spec.md:306,:354` (living spec — update; re-verify the
  :354 cited line number against the renamed file). README does not name the
  file (verified). Historical audit/experiment docs stay untouched.

## 3. Verification (report real results; machine assumed idle)
- `make unit` + `make test` — all green; plunit test count identical
  before/after the rename (record both numbers).
- `make bench-matrix` — CSV renders; all four strategies still resolve.
- ⊳ `swipl -q benchmarks/start_sensitivity.pl` — the other live `baseline`
  consumer, on no make target. If its full sweep is long, run it once and
  report the runtime honestly; do not silently skip.
- No-change sentinels: `swipl -q benchmarks/check_baseline.pl` and
  `swipl -q benchmarks/check_fill_baseline.pl` — every gated rung exactly
  +0.00%; do not re-record anything.
- `list_undefined` gate over load.pl + benchmark files.

## 4. Non-goals
- No deletion of `crossword/3,4`, the strategies, `run_matrix.pl`, or
  `start_sensitivity.pl` — the recorded keep decision stands.
- No golden, baseline, or history changes. No new exports (AGENTS.md: never
  add an export for a test).
- No module extraction / code motion in core.pl.
