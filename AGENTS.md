# AGENTS.md

This file is a lightweight context shim for coding agents. Project context
belongs in [README.md](README.md) so it stays useful to both humans and agents.
Add agent-only instructions here over time when they are not appropriate for
the README.

## Start Here

- Read [README.md](README.md) before making changes. It is the canonical source
  for the project overview, requirements, CLI usage, data structures, solver
  behavior, and implementation limitations.
- Read the relevant specs in [docs/](docs/) before changing JSON input or
  output behavior:
  - [docs/json-input-spec.md](docs/json-input-spec.md)
  - [docs/json-output-spec.md](docs/json-output-spec.md)
- Need an SWI-Prolog predicate's signature, mode, or semantics? Look it up in the
  local, version-matched manual under
  [docs/reference/swi-manual/](docs/reference/swi-manual/) (grep it, or see its
  `INDEX.md`) instead of guessing or going online — see
  [docs/reference/README.md](docs/reference/README.md).

## Project Shape

- The implementation is a SWI-Prolog crossword layout-and-fill toolkit.
- `crosswordsmith` is the CLI entry point (verbs: `arrange`, `lint`, `export`,
  `fill`). The implementation lives under `prolog/crosswordsmith/`, one module
  per file (`crosswordsmith_core` etc.) with explicit export lists; root
  `load.pl` loads all seven modules (the CLI, tests, and benchmarks all go
  through it) and defines the `crosswordsmith` file-search alias.
- Export lists carry only the real inter-module/CLI/benchmark API. White-box
  tests call module internals as `Module:Pred(...)` — never add an export
  for a test.
- `prolog/crosswordsmith/core.pl` is the shared substrate (grid model,
  free-canvas legality core, clue numbering, emit, input).
- The engines are `arrange.pl` (Flavour A, including the greedy constructor)
  and `fill.pl` / `lint.pl` / `export.pl` / `stockgrid.pl` (Flavour B);
  `metrics.pl` holds the shared metric predicates (all under
  `prolog/crosswordsmith/`).
- `fixtures/bundled_17_clues.pl` is the bundled clue dataset used by examples,
  benchmarks, and the golden regression.
- `tests/*.plt` hold the plunit coverage (one suite per module); the
  deterministic CLI goldens live under `tests/golden/`. Run all three layers
  (plunit + goldens + CLI exit-code checks) via `./run_tests.sh` or `make test`.

## Working Instructions

- Keep broad project explanations in `README.md`, not here.
- Keep behavioral changes covered by focused plunit tests and, when JSON output
  intentionally changes, update and review the golden fixture.
- Prefer the existing Prolog style and SWI-Prolog libraries already in use.
- Do not change the JSON input or output contract without updating the
  corresponding spec, README, tests, and golden output as needed.

## Verification

- Full test suite: `make test` or `./run_tests.sh`
- Plunit only: `make unit`
- Golden output only: `make golden`
- Regenerate golden output only after an intentional output change:
  `make update-golden`
