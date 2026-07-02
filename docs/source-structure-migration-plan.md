# Source Structure + Module Migration Plan

Status: planned.

This plan tracks a behaviour-preserving migration from the current top-level
plain `.pl` implementation files toward a conventional SWI-Prolog library
layout, followed by a second pass that introduces modules and explicit exports.

The migration is intentionally split in two:

1. Move files and loading paths first, with no predicate visibility changes.
2. Add `:- module(...)` declarations and explicit `use_module(...)` imports
   after the path churn has settled.

Every checked-off implementation step must preserve the public CLI contract and
keep `make test` green. JSON golden output should not change unless a separate,
intentional output-contract change is made.

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

## Target Shape

The intended end state is:

```text
crosswordsmith
crossword.pl                  # compatibility/migration shim only
load.pl                       # project load/search-path setup
prolog/
  crosswordsmith/
    core.pl                   # current core substrate from crossword.pl
    arrange.pl
    metrics.pl                # current quality.pl, renamed after relocation
    lint.pl
    export.pl
    stockgrid.pl
    fill.pl
tests/
fixtures/
grids/
docs/
```

The exact module names should use a project prefix because SWI module names live
in a flat namespace. Proposed names:

- `crosswordsmith_core`
- `crosswordsmith_arrange`
- `crosswordsmith_metrics`
- `crosswordsmith_lint`
- `crosswordsmith_export`
- `crosswordsmith_stockgrid`
- `crosswordsmith_fill`

## Current File Classification

| File | Current role | Target disposition |
|---|---|---|
| `crosswordsmith` | Active CLI script | Keep at repo root; load project via `load.pl` or search path. |
| `crossword.pl` | Active shared substrate plus legacy direct-CLI shim | Split: root file remains shim; implementation moves to `prolog/crosswordsmith/core.pl`. |
| `arrange.pl` | Active Flavour-A engine | Move to `prolog/crosswordsmith/arrange.pl`, then module-ize. |
| `quality.pl` | Active shared metrics + greedy constructor; legacy name | Move first, then rename to `metrics.pl` when imports are explicit. |
| `lint.pl` | Active validator | Move, then module-ize. |
| `export.pl` | Active ipuz/Exolve transformer | Move, then module-ize early. |
| `stockgrid.pl` | Active stock-grid library | Move, then module-ize early. |
| `fill.pl` | Active fill engine | Move, then module-ize after stockgrid/arrange dependencies are explicit. |

## Phase 1: Layout Migration Only

Goal: move source files under `prolog/crosswordsmith/` without changing predicate
visibility, module semantics, CLI behaviour, or emitted output.

- [ ] Add root `load.pl`.
  - [ ] Define a project file-search alias rooted at `prolog/`.
  - [ ] Load implementation files in the current known-good order.
  - [ ] Keep the loader side-effect profile compatible with the test harness.
- [ ] Create `prolog/crosswordsmith/`.
- [ ] Move active implementation files into `prolog/crosswordsmith/`.
  - [ ] Move `arrange.pl`.
  - [ ] Move `quality.pl` without renaming yet.
  - [ ] Move `lint.pl`.
  - [ ] Move `export.pl`.
  - [ ] Move `stockgrid.pl`.
  - [ ] Move `fill.pl`.
  - [ ] Move the implementation portion of `crossword.pl` into
        `prolog/crosswordsmith/core.pl`.
- [ ] Keep root `crossword.pl` as the compatibility/migration shim.
- [ ] Update `crosswordsmith` to load through `load.pl` or the project search
      path rather than top-level sibling files.
- [ ] Update `tests/run_tests.pl` to load `load.pl`.
- [ ] Update benchmark harnesses to load `load.pl`.
- [ ] Update documentation references that assume implementation files live at
      repo root.
- [ ] Run `make test`.
- [ ] Confirm no golden output changes.

## Phase 2: Rename `quality.pl` To `metrics.pl`

Goal: remove the stale file name after relocation, before module exports make
the name more permanent.

- [ ] Rename `prolog/crosswordsmith/quality.pl` to
      `prolog/crosswordsmith/metrics.pl`.
- [ ] Update load paths and comments.
- [ ] Update docs that describe `quality.pl` as the current shared metric layer.
- [ ] Run `make test`.
- [ ] Confirm no golden output changes.

## Phase 3: Introduce Modules From Leaves Inward

Goal: replace load-order coupling with explicit module boundaries and imports.
Module-ize low-coupling leaves first, then the shared substrate last.

### 3.1 `crosswordsmith_export`

- [ ] Add `:- module(crosswordsmith_export, [...]).`.
- [ ] Export:
  - [ ] `export_solve/2`
  - [ ] `layout_to_ipuz/2`
  - [ ] `layout_to_exolve/2`
  - [ ] `answer_enumeration/2`
- [ ] Update tests/imports.
- [ ] Run `make test`.

### 3.2 `crosswordsmith_stockgrid`

- [ ] Add `:- module(crosswordsmith_stockgrid, [...]).`.
- [ ] Export:
  - [ ] `stockgrid_load/2`
  - [ ] `stockgrid_validate/3`
  - [ ] `stockgrid_validate_file/3`
  - [ ] `stockgrid_report/1`
  - [ ] `mask_white_cells/3`
  - [ ] `grid_run/4`
- [ ] Import the lint/core APIs explicitly once available; use temporary
      compatibility loading if this step precedes `lint` module-ization.
- [ ] Update tests/imports.
- [ ] Run `make test`.

### 3.3 `crosswordsmith_lint`

- [ ] Add `:- module(crosswordsmith_lint, [...]).`.
- [ ] Export:
  - [ ] `lint_solve/4`
  - [ ] `lint_run/5`
  - [ ] `lint_load/3`
  - [ ] `lint_known_profile/1`
- [ ] Import metrics/core APIs explicitly once available.
- [ ] Update tests/imports.
- [ ] Run `make test`.

### 3.4 `crosswordsmith_fill`

- [ ] Add `:- module(crosswordsmith_fill, [...]).`.
- [ ] Export:
  - [ ] `fill_solve/4`
  - [ ] narrowly selected helper predicates only where tests or other modules
        establish a real API need.
- [ ] Import stockgrid, arrange fragment-loading, and core helpers explicitly.
- [ ] Update tests/imports.
- [ ] Run `make test`.

### 3.5 `crosswordsmith_metrics`

- [ ] Add `:- module(crosswordsmith_metrics, [...]).`.
- [ ] Export shared metric predicates used by arrange/lint:
  - [ ] `layout_dir_cells/2`
  - [ ] `word_checked_bitmap/3`
  - [ ] `word_checked_count/3`
  - [ ] `word_max_unch_run/3`
  - [ ] `word_meets_half/2`
  - [ ] `word_half_threshold/2`
  - [ ] `checked_cells/2`
  - [ ] `crossing_count/5`
  - [ ] `placed_bbox/4`
  - [ ] `word_cells/5`
- [ ] Export greedy-constructor predicates used by arrange:
  - [ ] `seed_candidates/2`
  - [ ] `greedy_construct/6`
- [ ] Import core helpers explicitly.
- [ ] Update arrange/lint/fill imports.
- [ ] Run `make test`.

### 3.6 `crosswordsmith_arrange`

- [ ] Add `:- module(crosswordsmith_arrange, [...]).`.
- [ ] Export CLI-facing predicates:
  - [ ] `arrange_solve/4`
  - [ ] `arrange_fragment_solve/5`
  - [ ] `arrange_candidates_solve/5`
  - [ ] `arrange_enumerate_solve/2`
  - [ ] `load_fragment/3` if `fill` continues to reuse it.
- [ ] Decide whether white-box test helpers are exported or accessed through
      module qualification in tests.
- [ ] Import core and metrics explicitly.
- [ ] Update CLI/tests/imports.
- [ ] Run `make test`.

### 3.7 `crosswordsmith_core`

- [ ] Add `:- module(crosswordsmith_core, [...]).`.
- [ ] Export only deliberate shared substrate:
  - [ ] clue/input helpers used by CLI/arrange
  - [ ] grid/search primitives used by arrange/metrics/benchmarks
  - [ ] clue numbering and canonical emit helpers
  - [ ] `with_output/2`
- [ ] Keep direct execution behaviour in the root `crossword.pl` shim, not in
      the implementation module.
- [ ] Update all imports.
- [ ] Update tests to import or qualify core predicates explicitly.
- [ ] Run `make test`.
- [ ] Run `make fuzz` as an additional determinism check.

## Phase 4: Cleanup + Documentation

- [ ] Remove temporary compatibility loads.
- [ ] Remove comments that describe "consult AFTER" ordering once imports make
      the dependency graph explicit.
- [ ] Update `README.md` project-shape section.
- [ ] Update `AGENTS.md` if the top-level source map changes.
- [ ] Update `docs/design-spec.md` target module-layout language.
- [ ] Update `docs/STATUS.md` and this checklist in the same commits as work
      lands.
- [ ] Run `make test`.
- [ ] Run `make fuzz`.

## Suggested Commit Boundaries

1. Add `load.pl` and search-path setup, no file moves.
2. Move files under `prolog/crosswordsmith/`, update loaders/tests/benchmarks.
3. Split root `crossword.pl` shim from implementation `core.pl`.
4. Rename `quality.pl` to `metrics.pl`.
5. Module-ize `export`.
6. Module-ize `stockgrid` and `lint`.
7. Module-ize `fill`.
8. Module-ize `metrics` and `arrange`.
9. Module-ize `core`.
10. Final docs/comment cleanup.

## Verification Gate For Each Commit

- [ ] `make test`
- [ ] No unintended golden diff
- [ ] No CLI migration-regression for root `crossword.pl`
- [ ] No change to the canonical JSON input/output contract unless specs and
      goldens are intentionally updated in the same commit

