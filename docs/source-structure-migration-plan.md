# Source Structure + Module Migration Plan

Status: planned. Stress-tested against the codebase on 2026-07-02: full
cross-file coupling map, white-box test audit, docs audit, and SWI-Prolog
10.1.10 mechanism checks. The dependency map and export lists below are
verified against actual clause heads and call sites, not aspirational.

This plan tracks a behaviour-preserving migration from the current top-level
plain `.pl` implementation files toward a conventional SWI-Prolog library
layout, followed by a second pass that introduces modules and explicit exports.

The migration is intentionally staged:

1. Align the design spec with this plan's end state (the spec is normative).
2. Move files and loading paths, with no predicate visibility changes.
3. Two small code-motion steps (rename `quality.pl`; extract the greedy
   constructor into `arrange.pl`) while everything is still module-free.
4. Add `:- module(...)` declarations and explicit `use_module(...)` imports,
   one file per commit, leaves inward.

Every checked-off implementation step must preserve the public CLI contract and
keep `make test` green. JSON golden output should not change unless a separate,
intentional output-contract change is made.

## Resolved Decisions

These were open questions; they are now settled (2026-07-02 review):

1. **Metrics stay a separate module; `quality.pl` is retired by rename, not
   absorption.** `docs/design-spec.md` §4 currently says the metric helpers
   get lifted into `crossword.pl` and `quality.pl` retired; this plan
   supersedes that (Phase 0 amends the spec). Rationale: `lint.pl` calls zero
   `crossword.pl` predicates today — its only project dependency is the metric
   layer, and folding metrics into core would force the validator to depend on
   the whole solver substrate; core's module-ization is already the riskiest
   step and shouldn't grow by ~15 exports; the spec's language predates this
   plan and was aimed at retiring the stale `--quality`-era *name*, which the
   rename to `metrics.pl` achieves.
2. **The greedy constructor moves into `arrange.pl`.** `quality.pl` holds two
   unrelated things: pure metric predicates and the greedy constructor
   (`seed_candidates/2`, `greedy_construct/6`, `greedy_loop/6`). The
   constructor's only consumer is arrange, and its project dependencies are
   core's *search* primitives, not the metrics. Moving it (Phase 3) keeps the
   metrics module honestly named and avoids exporting constructor internals
   across a module boundary for a single caller.
3. **White-box tests use module qualification, not exports.** The six `.plt`
   suites call ~50 internal predicates by unqualified name today (see the
   audit summary below). After module-ization they call internals as
   `Module:Pred(...)`. Verified on SWI 10.1.10: qualified calls reach
   non-exported predicates and module-local dynamics. Export lists therefore
   document only the real inter-module/CLI API, never test scaffolding.
4. **Root `crossword.pl` becomes a message-only shim.** Shebang +
   `legacy_main` migration message + `halt(1)`, roughly 15 lines; it no longer
   loads the library. Nothing in-repo consults it after Phase 1 (tests and
   benchmarks switch to `load.pl`). Side benefit: the
   `initialization(_, main)` last-registration-wins race between
   `crossword.pl` and the `crosswordsmith` script disappears, because
   `core.pl` carries no initialization directive.

## Research Basis

- SWI-Prolog project guidance says source organisation depends on project size:
  small projects may use one directory of linked files, while larger projects
  are normally split into sub-project directories with search paths and load
  files. See [`docs/reference/swi-manual/projectfiles.md`](reference/swi-manual/projectfiles.md).
- SWI-Prolog pack guidance defines a `prolog/` directory as the conventional
  library source root for packs; when a pack is attached, that directory is
  added to the `library` file search path. See
  [`docs/reference/swi-manual/pack-structure.md`](reference/swi-manual/pack-structure.md).
- SWI-Prolog modules are normally created by a source file whose first term is
  `:- module(Name, PublicList).`, making the exported predicate surface explicit.
  See [`docs/reference/swi-manual/defmodule.md`](reference/swi-manual/defmodule.md).
- PlUnit supports separate `.plt` files that load the files under test, so the
  existing `tests/` shape is already compatible with the target layout. See
  [`docs/reference/swi-manual/packages/plunit.md`](reference/swi-manual/packages/plunit.md).
- **Verified by experiment on SWI-Prolog 10.1.10** (the adopted baseline):
  (a) module code resolves predicates *defined in* `user` (modules inherit
  from `user`); (b) predicates *imported into* `user` via `use_module` are
  likewise visible to modules; (c) `M:pred(...)` reaches non-exported and
  dynamic predicates. (a)+(b) are the bridge that lets Phase 4 module-ize one
  file at a time while the rest stay plain; (c) is the white-box test
  mechanism.

## Verified Dependency Map (2026-07-02)

Load graph today: the `crosswordsmith` driver `ensure_loaded`s, in order,
`arrange.pl` (which chain-loads `crossword.pl`, which chain-loads
`quality.pl`), then `lint.pl`, `export.pl`, `stockgrid.pl`, `fill.pl`.
`lint`/`export`/`stockgrid`/`fill` perform **no project loads of their own** —
they rely entirely on the driver's (or `tests/run_tests.pl`'s) load order.
`load.pl` becomes the single owner of that order.

Cross-file runtime calls (caller → provider), the basis for every export list
below:

- `quality` → `crossword`: `init_grid/2`, `start_loc/4`, `remove_x/3`,
  `fits_on_grid/4`, `assign_word/10`, `find_intersecting_word/6`,
  `next_cell/4`. All but `next_cell/4` sit on the constructor side and move to
  arrange's import list in Phase 3 (re-verify the exact split at extraction).
- `arrange` → `crossword`: `start_locs/1`, `assign_clue_numbers/2`,
  `find_crossword/6`, `check_unique_answers/1`, `shares_letter/2`,
  `build_grid_rows/3`, `build_words/4`, `answer_meta_assoc/2`, `cell_coord/3`,
  `init_grid/2`, `assign_word/10`, `fits_on_grid/4`, `assign_words_inc/9`,
  `default_strategy/1`, `all_crossword/5`.
- `arrange` → `quality`: `word_checked_count/3`, `word_letters/3`,
  `placed_bbox/4`, `word_cells/5`, plus the constructor entry points
  (`seed_candidates/2`, `greedy_construct/6`, `greedy_loop/6`) until Phase 3
  internalizes them.
- `lint` → `quality`: `layout_dir_cells/2`, `word_checked_bitmap/3`,
  `bits_checked_count/2`, `bits_max_unch_run/2`, `word_half_threshold/2`,
  `cell_rc/4`. `lint` → `crossword`: **none** — the validator deliberately
  depends only on the JSON contract plus metrics. Preserve this boundary.
- `export` → nothing. Fully self-contained apart from library modules.
- `stockgrid` → `crossword`: `assign_clue_numbers/2`. `stockgrid` → `lint`:
  `lint_run/5` (this edge dictates module-izing lint *before* stockgrid).
- `fill` → `crossword`: `assign_clue_numbers/2`, `emit_json/3`. `fill` →
  `quality`: `word_letters/3`. `fill` → `arrange`: `load_fragment/3`. `fill` →
  `stockgrid`: `stockgrid_load/2`, `mask_white_cells/3`, `grid_run/4`.
- CLI driver → core: `with_output/2`, `load_clues/2`; → arrange: the four
  `*_solve` entry points, `load_fragment/3`, `reconcile_fragment_size/3`, and
  the check-target override (see below); → lint: `lint_solve/4`,
  `lint_known_profile/1`; → export: `export_solve/2`; → fill: `fill_solve/4`.
- Benchmarks (today consulting only `crossword.pl`): `find_crossword/6`,
  `default_strategy/1`, `valid_strategy/1`, `require_strategy/1`,
  `start_locs/1`. They read fixtures with their own `clues/1` readers, not
  `load_clues/2`.

Global-state and directive inventory (everything module-ization can trip on):

- `check_target_override/1` — the **only** `:- dynamic` in the seven files
  (arrange.pl). Written by the driver's `set_check_target/1`
  (retractall+assertz), read by arrange's `check_target/2`, and asserted as
  `user:check_target_override(_)` by `arrange.plt`. Phase 4.6 moves
  `set_check_target/1` into arrange and exports it; the dynamic itself stays
  module-private; driver and tests go through the setter (tests may also
  qualify).
- `script_directory/1` — asserted by the driver only; stays in the driver.
- `prolog:error_message//1` — multifile hook with clauses in six files. It is
  explicitly module-qualified and declared multifile, which is the intended
  cross-module extension pattern: **module-safe as-is, no change needed.**
- No `op`, `table`, `discontiguous`, or `nb_setval`/`b_setval` declarations
  anywhere; no ordinary predicate has clauses split across files.
- `:- meta_predicate` declarations exist in crossword.pl (`with_output/2`,
  `capped/2`) and export.pl (`collect_by_number/4`); meta arguments get module
  sensitivity once inside modules — keep the declarations with the moved code.

White-box test audit summary (what Phase 4 must convert to qualified calls):
`crossword.plt` calls ~20 crossword.pl internals (geometry, solver, loader
helpers) plus 12 quality.pl metric helpers; `arrange.plt` ~19 arrange.pl
internals plus `user:check_target_override/1`; `lint.plt` 5 lint internals plus
2 quality helpers; `export.plt` 5 export predicates; `stockgrid.plt` 6
stockgrid predicates plus `lint_run/5`; `fill.plt` ~11 fill internals plus
`lint_load/3`/`lint_run/5`. None of the `.plt` files loads anything itself —
`tests/run_tests.pl` is the single loader to update.

## Target Shape

The intended end state is:

```text
crosswordsmith
crossword.pl                  # message-only migration shim (~15 lines)
load.pl                       # project load/search-path setup (single owner of load order)
prolog/
  crosswordsmith/
    core.pl                   # current core substrate from crossword.pl
    arrange.pl                # Flavour-A engine + the greedy constructor (from quality.pl)
    metrics.pl                # current quality.pl minus the constructor, renamed
    lint.pl
    export.pl
    stockgrid.pl
    fill.pl
tests/
fixtures/
grids/
docs/
```

The exact module names use a project prefix because SWI module names live in a
flat namespace:

- `crosswordsmith_core`
- `crosswordsmith_arrange`
- `crosswordsmith_metrics`
- `crosswordsmith_lint`
- `crosswordsmith_export`
- `crosswordsmith_stockgrid`
- `crosswordsmith_fill`

`load.pl` defines a search alias (e.g.
`user:file_search_path(crosswordsmith, <load.pl dir>/prolog/crosswordsmith)`)
so loads read `crosswordsmith(core)` rather than long relative paths. The
layout is pack-*shaped* but this is not a pack; no `pack.pl` is planned.

## Current File Classification

| File | Current role | Target disposition |
|---|---|---|
| `crosswordsmith` | Active CLI script | Keep at repo root; load project via `load.pl`. |
| `crossword.pl` | Active shared substrate plus legacy direct-CLI shim | Split: root file becomes a message-only shim; implementation moves to `prolog/crosswordsmith/core.pl`. |
| `arrange.pl` | Active Flavour-A engine | Move, absorb the greedy constructor (Phase 3), then module-ize. |
| `quality.pl` | Active shared metrics + greedy constructor; legacy name | Move first; rename to `metrics.pl` (Phase 2); shed the constructor (Phase 3); module-ize. |
| `lint.pl` | Active validator | Move, then module-ize early (before stockgrid — stockgrid calls `lint_run/5`). |
| `export.pl` | Active ipuz/Exolve transformer | Move, then module-ize first (zero project dependencies). |
| `stockgrid.pl` | Active stock-grid library | Move, then module-ize after lint. |
| `fill.pl` | Active fill engine | Move, then module-ize after stockgrid (arrange/metrics/core resolve via `user` inheritance until their own steps). |

## Phase 0: Align The Design Spec

The spec is normative (spec-driven development); it must describe this end
state *before* the code moves.

- [x] Amend `docs/design-spec.md` §4 "Module layout (target)": replace the
      "metric helpers lifted into `crossword.pl` / `quality.pl` retired"
      language with the `prolog/crosswordsmith/{core,metrics,arrange,lint,
      export,stockgrid,fill}` module layout, the constructor-lives-in-arrange
      decision, and the lint-depends-only-on-metrics boundary.
- [x] Update `docs/STATUS.md` line ~131 (the deferred metric-relocation tidy-up)
      to point at this plan's Phases 2–3 instead (also aligned the §6.4
      metric-predicates row, which carried the same superseded language).

## Phase 1: Layout Migration Only

Goal: move source files under `prolog/crosswordsmith/` without changing
predicate visibility, module semantics, CLI behaviour, or emitted output.

- [x] Add root `load.pl`.
  - [x] Define the `crosswordsmith` file-search alias rooted at
        `prolog/crosswordsmith/`, relative to `load.pl`'s own directory
        (`prolog_load_context/2`), so it is cwd-independent. (Landed pointing
        at the repo root in the no-moves commit; repointed at
        `prolog/crosswordsmith/` in the move commit.)
  - [x] Load implementation files in the driver's current known-good order:
        arrange (chain-loads core → quality), lint, export, stockgrid, fill.
  - [x] Keep the loader side-effect free (no initialization directives, no
        flag changes beyond what the files themselves do) so consulting it
        from the test harness and benchmarks is safe.
- [x] Create `prolog/crosswordsmith/`.
- [x] Move active implementation files into `prolog/crosswordsmith/`.
  - [x] Move `arrange.pl`.
  - [x] Move `quality.pl` without renaming yet.
  - [x] Move `lint.pl`.
  - [x] Move `export.pl`.
  - [x] Move `stockgrid.pl`.
  - [x] Move `fill.pl`.
  - [x] Move the implementation portion of `crossword.pl` (~875 of 886 lines)
        into `prolog/crosswordsmith/core.pl`. `core.pl` gets **no**
        `initialization` directive and no shebang.
  - [x] Keep the intra-package sibling chain-loads during the transition
        (`arrange.pl` `ensure_loaded`s `core.pl`; `core.pl` `ensure_loaded`s
        `quality.pl` — same `prolog_load_context` idiom, paths unchanged since
        the files move together), so individual files stay consultable.
- [x] Reduce root `crossword.pl` to the message-only shim: shebang,
      `:- initialization(legacy_main, main)`, `legacy_main/0` migration
      message, `halt(1)`. It does **not** load the library.
- [x] Update `crosswordsmith` to load `load.pl` (resolved relative to the
      script directory) instead of the five sibling files. Keep the
      `script_directory/1` assert and the foot-of-file
      `:- initialization(main, main)`.
- [x] Update `tests/run_tests.pl`: replace the six implementation consults
      with `consult('load.pl')`; keep the six `.plt` consults unchanged.
- [x] Update benchmark harnesses (`run_benchmarks.pl`, `run_matrix.pl`,
      `start_sensitivity.pl`) to load `load.pl` via their existing
      repo-root-relative resolution, replacing their direct `crossword.pl`
      consults. (Heavier load than they need; acceptable and temporary until
      they can `use_module(crosswordsmith(core))` after Phase 4.7.)
- [x] Update documentation references that assume implementation files live at
      repo root: `README.md` source-map table (~lines 300–307), `AGENTS.md`
      "Project Shape" (~lines 25–34; `CONTEXT.md` is a symlink to it).
- [x] Run `make test`.
- [x] Run `make fuzz` — path resolution and the `--out` contract are exactly
      what this phase touches.
- [x] Confirm no golden output changes.

## Phase 2: Rename `quality.pl` To `metrics.pl`

Goal: remove the stale file name after relocation, before module exports make
the name more permanent.

- [x] Rename `prolog/crosswordsmith/quality.pl` to
      `prolog/crosswordsmith/metrics.pl`.
- [x] Update the `core.pl` chain-load and comments.
- [x] Update docs that describe `quality.pl` as the current shared metric layer.
- [x] Run `make test`.
- [x] Confirm no golden output changes.

## Phase 3: Extract The Greedy Constructor Into `arrange.pl`

Goal: pure code motion while everything is still module-free, so the metrics
module is honestly named and constructor internals never cross a module
boundary.

- [x] Move `seed_candidates/2`, `greedy_construct/6`, `greedy_loop/6` and their
      private helpers from `metrics.pl` into `arrange.pl`. (Helpers moved:
      `neg_answer_len/2`, `seed_word/7`, `apply_move/6`, `next_move/5`,
      `best_move/2`, `word_best_placement/8`, `placement_key/8`,
      `crossing_count/6`+`cc_/7`, `bbox_growth/4`, `extend_cell/4`.
      `word_letters/3` and `word_cells/5` stay in metrics — they have direct
      arrange/fill consumers and are on the 4.5 export list.)
- [x] Re-verify the metrics→core dependency split: the constructor takes
      `init_grid/2`, `start_loc/4`, `remove_x/3`, `fits_on_grid/4`,
      `assign_word/10`, `find_intersecting_word/6` with it; `metrics.pl`
      should retain little beyond `next_cell/4` from core (confirm at
      extraction). **Confirmed at extraction (2026-07-02): metrics retains
      exactly `next_cell/4` from core (via `word_cells/5`), nothing else.
      One new arrange→metrics edge: the moved `extend_cell/4` calls
      `cell_rc/4`, which stays in metrics and is already on the 4.5 export
      list — so 4.6's import of metrics must include `cell_rc/4`.**
- [x] Update `tests/crossword.plt`'s quality block: `seed_candidates/2` tests
      now target arrange (test moved to `arrange.plt`; `word_letters/3` test
      stays with metrics).
- [x] Run `make test`.
- [x] Confirm no golden output changes.

## Phase 4: Introduce Modules From Leaves Inward

Goal: replace load-order coupling with explicit module boundaries and imports.
Module-ize low-coupling leaves first, then the shared substrate last.

**Mechanism (verified on SWI 10.1.10):** at each step, `load.pl` swaps that
file's `consult`-style load for a `use_module` from `user` context. The new
module's exports land in `user`, so not-yet-module-ized plain files keep
resolving them; the module's own calls to still-plain predicates resolve via
inheritance from `user`. Already-module-ized consumers add an explicit
`use_module` of the new module in the same commit. Each step also converts
that module's white-box test calls to `Module:Pred(...)` qualified form.

Export-list policy: export only what another implementation file, the CLI
driver, the benchmarks, or a deliberate library API needs. Test-only
predicates are never exported — tests qualify.

### 4.1 `crosswordsmith_export`

Zero project dependencies; the safest first target.

- [ ] Add `:- module(crosswordsmith_export, [...]).`.
- [ ] Export (driver uses `export_solve/2`; the rest are deliberate library
      API, currently test-covered):
  - [ ] `export_solve/2`
  - [ ] `layout_to_ipuz/2`
  - [ ] `layout_to_exolve/2`
  - [ ] `answer_enumeration/2`
- [ ] Qualify `export.plt` internals (`export_load/2`).
- [ ] Run `make test`.

### 4.2 `crosswordsmith_lint`

Before stockgrid, which calls `lint_run/5`. Lint's metrics/core calls resolve
via `user` inheritance until 4.5/4.7.

- [ ] Add `:- module(crosswordsmith_lint, [...]).`.
- [ ] Export:
  - [ ] `lint_solve/4` (driver)
  - [ ] `lint_run/5` (stockgrid)
  - [ ] `lint_known_profile/1` (driver)
  - [ ] `lint_load/3` (deliberate API; today only tests use it)
- [ ] Qualify `lint.plt` internals (`lint_dict_layout/3`, `barred_max_unch/2`,
      `eval_word_rule/5`) and update `stockgrid.plt`/`fill.plt` cross-calls to
      import or qualify `lint_run/5`/`lint_load/3`.
- [ ] Run `make test`.

### 4.3 `crosswordsmith_stockgrid`

- [ ] Add `:- module(crosswordsmith_stockgrid, [...]).`.
- [ ] Import `crosswordsmith_lint` explicitly (for `lint_run/5`). Core's
      `assign_clue_numbers/2` resolves via `user` inheritance until 4.7.
- [ ] Export (fill uses the first three; the validate/report trio is
      deliberate library API):
  - [ ] `stockgrid_load/2`
  - [ ] `mask_white_cells/3`
  - [ ] `grid_run/4`
  - [ ] `stockgrid_validate/3`
  - [ ] `stockgrid_validate_file/3`
  - [ ] `stockgrid_report/1`
- [ ] Qualify `stockgrid.plt` internals (`grid_lights/3`).
- [ ] Run `make test`.

### 4.4 `crosswordsmith_fill`

- [ ] Add `:- module(crosswordsmith_fill, [...]).`.
- [ ] Export:
  - [ ] `fill_solve/4` — nothing else; every other fill predicate is internal.
- [ ] Import `crosswordsmith_stockgrid` explicitly. Arrange's
      `load_fragment/3`, metrics' `word_letters/3`, and core's
      `assign_clue_numbers/2`/`emit_json/3` resolve via `user` inheritance
      until 4.5–4.7.
- [ ] Qualify `fill.plt` internals (`fill_grid/4`, `apply_seed/3`,
      `seeded_slot/2`, `load_dict/3`, `fill_attempt/7,8`, `select_mrv/6`,
      `candidate_count/4`, `candidates/4`, ...).
- [ ] Run `make test`.

### 4.5 `crosswordsmith_metrics`

- [ ] Add `:- module(crosswordsmith_metrics, [...]).`.
- [ ] Export the metric predicates with real cross-module consumers:
  - [ ] `layout_dir_cells/2` (lint)
  - [ ] `word_checked_bitmap/3` (lint)
  - [ ] `bits_checked_count/2` (lint)
  - [ ] `bits_max_unch_run/2` (lint)
  - [ ] `word_half_threshold/2` (lint)
  - [ ] `cell_rc/4` (lint)
  - [ ] `word_checked_count/3` (arrange)
  - [ ] `word_letters/3` (arrange, fill)
  - [ ] `placed_bbox/4` (arrange)
  - [ ] `word_cells/5` (arrange)
- [ ] Do **not** export: `crossing_count/6` (internal-only — note: the earlier
      draft of this plan listed it as `/5`, which was wrong on both counts),
      `checked_cells/2`, `word_meets_half/2`, `word_max_unch_run/3`,
      `dir_cells/3` (test-only — tests qualify).
- [ ] Update `crosswordsmith_lint` and `crosswordsmith_fill` to import
      metrics explicitly; `load.pl`'s `use_module` keeps plain `arrange.pl`
      resolving until 4.6.
- [ ] Qualify the `crossword.plt` quality-block calls.
- [ ] Run `make test`.

### 4.6 `crosswordsmith_arrange`

- [ ] Add `:- module(crosswordsmith_arrange, [...]).`.
- [ ] Export CLI-facing predicates:
  - [ ] `arrange_solve/4`
  - [ ] `arrange_fragment_solve/5`
  - [ ] `arrange_candidates_solve/5`
  - [ ] `arrange_enumerate_solve/2`
  - [ ] `load_fragment/3` (driver and fill)
  - [ ] `reconcile_fragment_size/3` (driver — missing from the earlier draft)
  - [ ] `set_check_target/1` (new home; see next item)
- [ ] Move `set_check_target/1` from the `crosswordsmith` driver into this
      module; keep `:- dynamic check_target_override/1` module-private. The
      driver calls the exported setter; `arrange.plt` drops its
      `user:check_target_override/1` asserts in favour of the setter (or
      qualified asserts).
- [ ] Import `crosswordsmith_metrics` and (once 4.7 lands) core explicitly;
      until 4.7, core predicates resolve via `user` inheritance.
- [ ] White-box test helpers (`arrange_best_layout/5,6`, `emit_arrange/4`,
      `fragment_dict_words/3`, `construct_one/7`, ...) are **not** exported —
      qualify in `arrange.plt` (per Resolved Decision 3).
- [ ] Update CLI/`fill`/tests.
- [ ] Run `make test`.

### 4.7 `crosswordsmith_core`

The largest export surface; last so every consumer's needs are already
explicit. The list below is the verified union of all consumers (metrics,
arrange, stockgrid, fill, driver, benchmarks):

- [ ] Add `:- module(crosswordsmith_core, [...]).`.
- [ ] Export:
  - [ ] I/O + emit: `with_output/2`, `load_clues/2`, `emit_json/3`
  - [ ] Clue numbering + layout build: `assign_clue_numbers/2`,
        `build_grid_rows/3`, `build_words/4`, `answer_meta_assoc/2`
  - [ ] Grid geometry: `cell_coord/3`, `init_grid/2`, `start_loc/4`,
        `start_locs/1`, `next_cell/4`, `fits_on_grid/4`
  - [ ] Search primitives: `assign_word/10`, `find_intersecting_word/6`,
        `assign_words_inc/9`, `find_crossword/6`, `all_crossword/5`
  - [ ] Strategy registry: `default_strategy/1`, `valid_strategy/1`,
        `require_strategy/1`
  - [ ] Utilities: `remove_x/3`, `shares_letter/2`, `check_unique_answers/1`
- [ ] Keep the `:- meta_predicate` declarations (`with_output/2`, `capped/2`)
      with the code.
- [ ] Update all module imports (metrics, arrange, stockgrid, fill) and the
      driver; switch benchmarks to `use_module(crosswordsmith(core))` if the
      lean load is wanted, else leave them on `load.pl`.
- [ ] Qualify `crossword.plt` internals (`swap_dir/2`, `calc_num/5`,
      `calc_start/5`, `prev_cell/4`, `is_start_cell/3`, `is_end_cell/3`,
      `doc_to_words/2`, `read_clues_json/2`, `valid_loc/1`, ...).
- [ ] Run `make test`.
- [ ] Run `make fuzz` as an additional determinism check.

## Phase 5: Cleanup + Documentation

- [ ] Remove the temporary intra-package `ensure_loaded` chain-loads and any
      remaining `use_module`-into-`user` compatibility loads in `load.pl`
      (its job shrinks to alias + `use_module` of the top-level modules).
- [ ] Remove comments that describe "consult AFTER" ordering once imports make
      the dependency graph explicit.
- [ ] Update `README.md` project-shape section.
- [ ] Update `AGENTS.md` source map (again, if 4.x changed anything beyond
      Phase 1's update; `CONTEXT.md` follows via symlink).
- [ ] Confirm `docs/design-spec.md` §4 matches the landed reality (Phase 0
      wrote the intent; this confirms it).
- [ ] Update `docs/STATUS.md` and this checklist in the same commits as work
      lands.
- [ ] Run `make test`.
- [ ] Run `make fuzz`.

## Suggested Commit Boundaries

1. Phase 0: design-spec §4 + STATUS.md alignment.
2. Add `load.pl` and search-path setup, no file moves.
3. Move files under `prolog/crosswordsmith/`, split the root `crossword.pl`
   shim from `core.pl`, update loaders/tests/benchmarks/docs.
4. Rename `quality.pl` to `metrics.pl`.
5. Extract the greedy constructor into `arrange.pl`.
6. Module-ize `export`.
7. Module-ize `lint`.
8. Module-ize `stockgrid`.
9. Module-ize `fill`.
10. Module-ize `metrics`.
11. Module-ize `arrange` (including the `set_check_target/1` relocation).
12. Module-ize `core`.
13. Final docs/comment cleanup.

## Verification Gate For Each Commit

- [ ] `make test`
- [ ] No unintended golden diff
- [ ] Root `crossword.pl` still prints the migration message and exits
      non-zero when run directly
- [ ] No change to the canonical JSON input/output contract unless specs and
      goldens are intentionally updated in the same commit
- [ ] For each Phase-4 step: that module's white-box test calls converted to
      qualified form in the same commit, and no test-only predicate appears in
      an export list
