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
- `xword/` is the Python conversion companion (terminal viewer + format
  multitool: native ⇄ ipuz/exolve, rendering). Its tests are NOT run by
  `make test` — use `make test-xword`. The engine↔xword cross-check is
  specified in [docs/xword-spec.md](docs/xword-spec.md) §11/§14 and has a
  one-command harness: `make xword-parity` (see Verification). Engine↔xword
  byte-parity is best-effort policy, NOT a contract (§14): engine JSON goes
  through SWI's `json_write_dict` (type-dependent colon spacing) while xword
  uses `json.dumps` — don't hand-roll probes to chase the residual gap.

## Browser / WASM

- The browser build lives under `wasm/` — read [wasm/README.md](wasm/README.md)
  first; it is canonical for the build, the wire contract, and the test battery.
- Build: `wasm/build/build-wasm.sh` (pins in `wasm/build/pins.sh`). It runs for
  MINUTES (emsdk/cmake/ninja) — run it as a background task, never under the
  default 2-minute Bash timeout. Same for any cmake/ninja step in
  `~/src/swipl-devel`.
- Test: `make test-wasm`. It stages swipl-web artifacts, qcompiles the app qlf,
  and self-manages a static server on a free port — do not hand-stage
  `swipl-web.*`/qlf files and do not hand-roll `http.server`/`pkill`.
- Never `git checkout` or dirty the shared `~/src/swipl-devel` tree; the wasm
  build uses its own isolated `build.wasm/`.
- Gotcha: the `USE_GMP=OFF` wasm build seeds a DIFFERENT RNG than native SWI,
  so `set_random(seed(N))` diverges CLI-vs-browser. The seeded arrange path
  therefore owns an engine-internal xorshift PRNG — read the PRNG section of
  [wasm/README.md](wasm/README.md) before touching seeding.
- Build/bench logs can exceed the 256KB Read cap — tail or grep them, never
  read them whole. To wait on a long-running step, use a background task or a
  Monitor-style until-loop, not `sleep N; tail`.

## Working Instructions

- Keep broad project explanations in `README.md`, not here.
- Keep behavioral changes covered by focused plunit tests and, when JSON output
  intentionally changes, update and review the golden fixture.
- Prefer the existing Prolog style and SWI-Prolog libraries already in use.
- Every exported predicate carries a PlDoc `%!` header with mode annotations
  (`+`/`-`/`?`) and a determinism tag (`det`/`semidet`/`nondet`/`multi`).
  The tags are load-bearing, verified claims (probed with `deterministic/1`
  and failure-path checks), not aspirations — when you add or change an
  export, verify its determinism before writing the tag, and keep the header
  in sync. Internal predicates use plain `%` prose as needed.
- Do not change the JSON input or output contract without updating the
  corresponding spec, README, tests, and golden output as needed.
- After a context compaction, Read a file before your first Edit to it — the
  "already read" state does not survive the summary, even for files you
  authored earlier in the session.
- On long build/spike sessions, land durable state EARLY: extract the working
  recipe into a script plus a ≤1-screen README before the first compaction,
  rather than accumulating it in a growing prose plan. Post-compaction
  recovery then reloads a 40-line script instead of re-reading a 777-line
  document (which one session did 34 times).
- Items marked best-effort / optional / Stretch carry an implicit STOP
  condition: if the happy path fails, document the gap and move on — that IS
  the completed outcome. Do not research or build probes to close the gap,
  and do not escalate it as a decision unless the plan is genuinely ambiguous.
- Ask nuanced tradeoff questions in prose (or state "proceeding with X unless
  you object"); reserve AskUserQuestion menus for clean either/or forks with
  stable framing.

## Verification

- Full native suite: `make test` or `./run_tests.sh` (plunit + goldens + CLI
  exit-code checks).
- Fast inner loop: `make unit` (plunit only) — prefer this while iterating on
  Prolog; run `make test` before committing.
- Golden output only: `make golden`. Regenerate only after an intentional
  output change (`make update-golden`) and review the diff.
- WASM battery: `make test-wasm` (see Browser / WASM above).
- xword companion: `make test-xword`. Engine↔xword output comparison:
  `make xword-parity [XWORD_PARITY_FIXTURE=<layout.json>]` — one command;
  don't hand-roll export/convert/diff probes.
